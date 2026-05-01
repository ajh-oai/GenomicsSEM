use crate::{
    munge_qc, sumstats_qc, KernelError, KernelResult, MungeQcOutput, SumstatsQcOutput,
    MUNGE_QC_COUNT_LEN, SUMSTATS_QC_COUNT_LEN,
};
use flate2::{read::GzDecoder, write::GzEncoder, Compression};
use rayon::prelude::*;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};

pub const MUNGE_FUSED_COLS: usize = 8;
pub const SUMSTATS_FUSED_COLS: usize = 9;

const SNP_COL: usize = 0;
const A1_COL: usize = 1;
const A2_COL: usize = 2;
const EFFECT_COL: usize = 3;
const P_COL: usize = 4;
const N_COL: usize = 5;
const INFO_COL: usize = 6;
const MAF_COL: usize = 7;
const SE_COL: usize = 8;

#[derive(Clone, Debug)]
struct MungeCandidate {
    snp: String,
    n: f64,
    a1_ref: i32,
    a2_ref: i32,
    a1_file: i32,
    a2_file: i32,
    effect: f64,
    p: f64,
    info: f64,
    maf: f64,
}

#[derive(Clone, Debug)]
struct SumstatsCandidate {
    a1_ref: i32,
    a2_ref: i32,
    a1_file: i32,
    a2_file: i32,
    effect: f64,
    se: f64,
    p: f64,
    n: f64,
    maf_ref: f64,
    maf_file: f64,
    info: f64,
}

#[derive(Clone, Debug)]
pub struct MungeFusedReport {
    pub unsupported: bool,
    pub rows_total: usize,
    pub rows_joined: usize,
    pub rows_written: usize,
    pub counts: [i32; MUNGE_QC_COUNT_LEN],
}

impl Default for MungeFusedReport {
    fn default() -> Self {
        Self {
            unsupported: false,
            rows_total: 0,
            rows_joined: 0,
            rows_written: 0,
            counts: [0; MUNGE_QC_COUNT_LEN],
        }
    }
}

#[derive(Clone, Debug)]
pub struct SumstatsFusedReport {
    pub unsupported: bool,
    pub rows_total: usize,
    pub rows_duplicate_removed: usize,
    pub rows_joined: usize,
    pub rows_written: usize,
    pub counts: [i32; SUMSTATS_QC_COUNT_LEN],
}

impl Default for SumstatsFusedReport {
    fn default() -> Self {
        Self {
            unsupported: false,
            rows_total: 0,
            rows_duplicate_removed: 0,
            rows_joined: 0,
            rows_written: 0,
            counts: [0; SUMSTATS_QC_COUNT_LEN],
        }
    }
}

pub struct SumstatsFusedOutput<'a> {
    pub keep: &'a mut [i32],
    pub beta: &'a mut [f64],
    pub se: &'a mut [f64],
}

#[derive(Clone, Debug)]
pub struct MungeFusedInput<'a> {
    pub filename: &'a str,
    pub output_path: &'a str,
    pub col_indices: &'a [i32],
    pub provided_n: Option<f64>,
    pub n_multiplier: f64,
}

fn open_text_reader(path: &str) -> KernelResult<Box<dyn BufRead>> {
    let file = File::open(path).map_err(|_| KernelError::BadDimensions)?;
    if path.to_ascii_lowercase().ends_with(".gz") {
        return Ok(Box::new(BufReader::new(GzDecoder::new(file))));
    }

    Ok(Box::new(BufReader::new(file)))
}

#[inline]
fn col_index(col_indices: &[i32], col: usize) -> Option<usize> {
    col_indices
        .get(col)
        .and_then(|idx| (*idx >= 0).then_some(*idx as usize))
}

#[inline]
fn field<'a>(fields: &'a [&str], idx: Option<usize>) -> Option<&'a str> {
    idx.and_then(|i| fields.get(i).copied())
}

#[inline]
fn parse_f64(value: Option<&str>) -> f64 {
    let Some(value) = value else {
        return f64::NAN;
    };
    let value = value.trim();
    if value.is_empty() || value == "." || value.eq_ignore_ascii_case("NA") {
        return f64::NAN;
    }

    value.parse::<f64>().unwrap_or(f64::NAN)
}

#[inline]
fn allele_code(value: Option<&str>) -> i32 {
    let Some(value) = value else {
        return 0;
    };
    let value = value.trim();
    if value.len() != 1 {
        return 0;
    }

    match value.as_bytes()[0].to_ascii_uppercase() {
        b'A' => 1,
        b'C' => 2,
        b'G' => 3,
        b'T' => 4,
        _ => 0,
    }
}

#[inline]
fn allele_label(code: i32) -> &'static str {
    match code {
        1 => "A",
        2 => "C",
        3 => "G",
        4 => "T",
        _ => "NA",
    }
}

fn required_indices_present(col_indices: &[i32], required: &[usize]) -> bool {
    required
        .iter()
        .all(|col| col_index(col_indices, *col).is_some())
}

fn build_ref_map<'a>(ref_snps: &'a [&'a str]) -> Option<HashMap<&'a str, usize>> {
    let mut ref_map = HashMap::with_capacity(ref_snps.len());
    for (idx, snp) in ref_snps.iter().enumerate() {
        if snp.is_empty() || ref_map.insert(*snp, idx).is_some() {
            return None;
        }
    }
    Some(ref_map)
}

fn write_munge_rows<W: Write>(
    writer: &mut W,
    candidates: &[MungeCandidate],
    keep: &[i32],
    z: &[f64],
    out_n: usize,
) -> KernelResult<()> {
    writeln!(writer, "SNP\tN\tZ\tA1\tA2").map_err(|_| KernelError::BadDimensions)?;
    for out_idx in 0..out_n {
        let row_idx = (keep[out_idx] - 1) as usize;
        let row = &candidates[row_idx];
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}",
            row.snp,
            row.n,
            z[out_idx],
            allele_label(row.a1_ref),
            allele_label(row.a2_ref)
        )
        .map_err(|_| KernelError::BadDimensions)?;
    }
    Ok(())
}

fn write_munge_output(
    output_path: &str,
    candidates: &[MungeCandidate],
    keep: &[i32],
    z: &[f64],
    out_n: usize,
) -> KernelResult<()> {
    let out_file = File::create(output_path).map_err(|_| KernelError::BadDimensions)?;
    if output_path.to_ascii_lowercase().ends_with(".gz") {
        let writer = BufWriter::new(out_file);
        let mut encoder = GzEncoder::new(writer, Compression::fast());
        write_munge_rows(&mut encoder, candidates, keep, z, out_n)?;
        let mut writer = encoder.finish().map_err(|_| KernelError::BadDimensions)?;
        writer.flush().map_err(|_| KernelError::BadDimensions)?;
        return Ok(());
    }

    let mut writer = BufWriter::new(out_file);
    write_munge_rows(&mut writer, candidates, keep, z, out_n)?;
    writer.flush().map_err(|_| KernelError::BadDimensions)?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn munge_fused_with_ref_map(
    filename: &str,
    output_path: &str,
    ref_snps: &[&str],
    ref_a1: &[i32],
    ref_a2: &[i32],
    ref_map: &HashMap<&str, usize>,
    col_indices: &[i32],
    provided_n: Option<f64>,
    n_multiplier: f64,
    info_filter: f64,
    maf_filter: f64,
) -> KernelResult<MungeFusedReport> {
    if ref_snps.len() != ref_a1.len()
        || ref_snps.len() != ref_a2.len()
        || col_indices.len() < MUNGE_FUSED_COLS
        || !required_indices_present(col_indices, &[SNP_COL, A1_COL, A2_COL, EFFECT_COL, P_COL])
        || (!provided_n.is_some_and(|x| x.is_finite()) && col_index(col_indices, N_COL).is_none())
    {
        return Ok(MungeFusedReport {
            unsupported: true,
            ..MungeFusedReport::default()
        });
    }

    let snp_idx = col_index(col_indices, SNP_COL);
    let a1_idx = col_index(col_indices, A1_COL);
    let a2_idx = col_index(col_indices, A2_COL);
    let effect_idx = col_index(col_indices, EFFECT_COL);
    let p_idx = col_index(col_indices, P_COL);
    let n_idx = col_index(col_indices, N_COL);
    let info_idx = col_index(col_indices, INFO_COL);
    let maf_idx = col_index(col_indices, MAF_COL);

    let mut reader = open_text_reader(filename)?;
    let mut line = String::new();
    reader
        .read_line(&mut line)
        .map_err(|_| KernelError::BadDimensions)?;

    let mut report = MungeFusedReport::default();
    let mut seen_joined = HashSet::new();
    let mut candidates = Vec::new();

    line.clear();
    while reader
        .read_line(&mut line)
        .map_err(|_| KernelError::BadDimensions)?
        != 0
    {
        if line.trim().is_empty() {
            line.clear();
            continue;
        }

        report.rows_total += 1;
        let fields: Vec<&str> = line.split_ascii_whitespace().collect();
        let Some(snp) = field(&fields, snp_idx) else {
            line.clear();
            continue;
        };
        let Some(&ref_idx) = ref_map.get(snp) else {
            line.clear();
            continue;
        };

        if !seen_joined.insert(snp.to_owned()) {
            report.unsupported = true;
            return Ok(report);
        }

        let mut n = provided_n.unwrap_or_else(|| parse_f64(field(&fields, n_idx)) * n_multiplier);
        if provided_n.is_some() {
            n = provided_n.unwrap();
        }

        let maf = parse_f64(field(&fields, maf_idx));
        let maf = if maf > 0.5 { 1.0 - maf } else { maf };

        candidates.push(MungeCandidate {
            snp: ref_snps[ref_idx].to_owned(),
            n,
            a1_ref: ref_a1[ref_idx],
            a2_ref: ref_a2[ref_idx],
            a1_file: allele_code(field(&fields, a1_idx)),
            a2_file: allele_code(field(&fields, a2_idx)),
            effect: parse_f64(field(&fields, effect_idx)),
            p: parse_f64(field(&fields, p_idx)),
            info: parse_f64(field(&fields, info_idx)),
            maf,
        });

        line.clear();
    }

    if candidates
        .iter()
        .any(|row| row.p.is_finite() && (row.p < 0.0 || row.p > 1.0))
    {
        report.unsupported = true;
        return Ok(report);
    }

    candidates.sort_by(|a, b| a.snp.cmp(&b.snp));
    report.rows_joined = candidates.len();

    let n = candidates.len();
    let a1_ref: Vec<i32> = candidates.iter().map(|row| row.a1_ref).collect();
    let a2_ref: Vec<i32> = candidates.iter().map(|row| row.a2_ref).collect();
    let a1_file: Vec<i32> = candidates.iter().map(|row| row.a1_file).collect();
    let a2_file: Vec<i32> = candidates.iter().map(|row| row.a2_file).collect();
    let effect: Vec<f64> = candidates.iter().map(|row| row.effect).collect();
    let p: Vec<f64> = candidates.iter().map(|row| row.p).collect();
    let info: Vec<f64> = candidates.iter().map(|row| row.info).collect();
    let maf: Vec<f64> = candidates.iter().map(|row| row.maf).collect();

    let mut keep = vec![0; n];
    let mut z = vec![0.0; n];
    let mut counts = [0; MUNGE_QC_COUNT_LEN];
    let info_slice = info_idx.map(|_| info.as_slice());
    let maf_slice = maf_idx.map(|_| maf.as_slice());
    let out_n = munge_qc(
        &a1_ref,
        &a2_ref,
        &a1_file,
        &a2_file,
        &effect,
        &p,
        info_slice,
        maf_slice,
        info_filter,
        maf_filter,
        &mut MungeQcOutput {
            keep: &mut keep,
            z: &mut z,
            counts: &mut counts,
        },
    )?;

    write_munge_output(output_path, &candidates, &keep, &z, out_n)?;
    report.rows_written = out_n;
    report.counts = counts;
    Ok(report)
}

#[allow(clippy::too_many_arguments)]
pub fn munge_fused(
    filename: &str,
    output_path: &str,
    ref_snps: &[&str],
    ref_a1: &[i32],
    ref_a2: &[i32],
    col_indices: &[i32],
    provided_n: Option<f64>,
    n_multiplier: f64,
    info_filter: f64,
    maf_filter: f64,
) -> KernelResult<MungeFusedReport> {
    let Some(ref_map) = build_ref_map(ref_snps) else {
        return Ok(MungeFusedReport {
            unsupported: true,
            ..MungeFusedReport::default()
        });
    };

    munge_fused_with_ref_map(
        filename,
        output_path,
        ref_snps,
        ref_a1,
        ref_a2,
        &ref_map,
        col_indices,
        provided_n,
        n_multiplier,
        info_filter,
        maf_filter,
    )
}

#[allow(clippy::too_many_arguments)]
pub fn munge_fused_batch(
    inputs: &[MungeFusedInput<'_>],
    ref_snps: &[&str],
    ref_a1: &[i32],
    ref_a2: &[i32],
    info_filter: f64,
    maf_filter: f64,
    n_threads: usize,
) -> KernelResult<Vec<MungeFusedReport>> {
    if n_threads == 0 {
        return Err(KernelError::BadDimensions);
    }
    let Some(ref_map) = build_ref_map(ref_snps) else {
        return Ok(inputs
            .iter()
            .map(|_| MungeFusedReport {
                unsupported: true,
                ..MungeFusedReport::default()
            })
            .collect());
    };

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .map_err(|_| KernelError::ThreadPoolBuild)?;
    pool.install(|| {
        inputs
            .par_iter()
            .map(|input| {
                munge_fused_with_ref_map(
                    input.filename,
                    input.output_path,
                    ref_snps,
                    ref_a1,
                    ref_a2,
                    &ref_map,
                    input.col_indices,
                    input.provided_n,
                    input.n_multiplier,
                    info_filter,
                    maf_filter,
                )
            })
            .collect()
    })
}

#[allow(clippy::too_many_arguments)]
pub fn sumstats_fused(
    filename: &str,
    ref_snps: &[&str],
    ref_a1: &[i32],
    ref_a2: &[i32],
    ref_maf: &[f64],
    col_indices: &[i32],
    provided_n: Option<f64>,
    info_filter: f64,
    ols: bool,
    beta_is_character: bool,
    linprob: bool,
    se_logit: bool,
    out: &mut SumstatsFusedOutput<'_>,
) -> KernelResult<SumstatsFusedReport> {
    if ref_snps.len() != ref_a1.len()
        || ref_snps.len() != ref_a2.len()
        || ref_snps.len() != ref_maf.len()
        || out.keep.len() < ref_snps.len()
        || out.beta.len() < ref_snps.len()
        || out.se.len() < ref_snps.len()
        || col_indices.len() < SUMSTATS_FUSED_COLS
        || !required_indices_present(col_indices, &[SNP_COL, A1_COL, A2_COL, EFFECT_COL, P_COL])
        || (!provided_n.is_some_and(|x| x.is_finite())
            && (ols || linprob)
            && col_index(col_indices, N_COL).is_none())
        || (!(ols && beta_is_character) && col_index(col_indices, SE_COL).is_none())
    {
        return Ok(SumstatsFusedReport {
            unsupported: true,
            ..SumstatsFusedReport::default()
        });
    }

    let mut ref_map = HashMap::with_capacity(ref_snps.len());
    for (idx, snp) in ref_snps.iter().enumerate() {
        if snp.is_empty() || ref_map.insert(*snp, idx).is_some() {
            return Ok(SumstatsFusedReport {
                unsupported: true,
                ..SumstatsFusedReport::default()
            });
        }
    }

    let snp_idx = col_index(col_indices, SNP_COL);
    let a1_idx = col_index(col_indices, A1_COL);
    let a2_idx = col_index(col_indices, A2_COL);
    let effect_idx = col_index(col_indices, EFFECT_COL);
    let p_idx = col_index(col_indices, P_COL);
    let n_idx = col_index(col_indices, N_COL);
    let info_idx = col_index(col_indices, INFO_COL);
    let maf_idx = col_index(col_indices, MAF_COL);
    let se_idx = col_index(col_indices, SE_COL);

    let mut reader = open_text_reader(filename)?;
    let mut line = String::new();
    reader
        .read_line(&mut line)
        .map_err(|_| KernelError::BadDimensions)?;

    let mut report = SumstatsFusedReport::default();
    let mut seen_counts: HashMap<String, usize> = HashMap::new();
    let mut ref_rows: Vec<Option<SumstatsCandidate>> = std::iter::repeat_with(|| None)
        .take(ref_snps.len())
        .collect();

    line.clear();
    while reader
        .read_line(&mut line)
        .map_err(|_| KernelError::BadDimensions)?
        != 0
    {
        if line.trim().is_empty() {
            line.clear();
            continue;
        }

        report.rows_total += 1;
        let fields: Vec<&str> = line.split_ascii_whitespace().collect();
        let Some(snp) = field(&fields, snp_idx) else {
            line.clear();
            continue;
        };
        let snp = snp.to_owned();
        let count = seen_counts.entry(snp.clone()).or_insert(0);
        *count += 1;
        if *count == 2 {
            report.rows_duplicate_removed += 2;
            if let Some(&ref_idx) = ref_map.get(snp.as_str()) {
                ref_rows[ref_idx] = None;
            }
            line.clear();
            continue;
        }
        if *count > 2 {
            report.rows_duplicate_removed += 1;
            line.clear();
            continue;
        }

        let Some(&ref_idx) = ref_map.get(snp.as_str()) else {
            line.clear();
            continue;
        };
        let n = provided_n.unwrap_or_else(|| parse_f64(field(&fields, n_idx)));

        ref_rows[ref_idx] = Some(SumstatsCandidate {
            a1_ref: ref_a1[ref_idx],
            a2_ref: ref_a2[ref_idx],
            a1_file: allele_code(field(&fields, a1_idx)),
            a2_file: allele_code(field(&fields, a2_idx)),
            effect: parse_f64(field(&fields, effect_idx)),
            se: parse_f64(field(&fields, se_idx)),
            p: parse_f64(field(&fields, p_idx)),
            n,
            maf_ref: ref_maf[ref_idx],
            maf_file: parse_f64(field(&fields, maf_idx)),
            info: parse_f64(field(&fields, info_idx)),
        });

        line.clear();
    }

    let mut candidates = Vec::with_capacity(ref_rows.len());
    let mut ref_indices = Vec::with_capacity(ref_rows.len());
    for (idx, row) in ref_rows.into_iter().enumerate() {
        if let Some(row) = row {
            candidates.push(row);
            ref_indices.push(idx + 1);
        }
    }
    report.rows_joined = candidates.len();

    if candidates
        .iter()
        .any(|row| row.p.is_finite() && (row.p < 0.0 || row.p > 1.0))
    {
        report.unsupported = true;
        return Ok(report);
    }

    let n = candidates.len();
    let a1_ref: Vec<i32> = candidates.iter().map(|row| row.a1_ref).collect();
    let a2_ref: Vec<i32> = candidates.iter().map(|row| row.a2_ref).collect();
    let a1_file: Vec<i32> = candidates.iter().map(|row| row.a1_file).collect();
    let a2_file: Vec<i32> = candidates.iter().map(|row| row.a2_file).collect();
    let effect: Vec<f64> = candidates.iter().map(|row| row.effect).collect();
    let se: Vec<f64> = candidates.iter().map(|row| row.se).collect();
    let p: Vec<f64> = candidates.iter().map(|row| row.p).collect();
    let n_values: Vec<f64> = candidates.iter().map(|row| row.n).collect();
    let maf_ref: Vec<f64> = candidates.iter().map(|row| row.maf_ref).collect();
    let maf_file: Vec<f64> = candidates.iter().map(|row| row.maf_file).collect();
    let info: Vec<f64> = candidates.iter().map(|row| row.info).collect();

    let mut keep = vec![0; n];
    let mut beta = vec![0.0; n];
    let mut se_out = vec![0.0; n];
    let mut counts = [0; SUMSTATS_QC_COUNT_LEN];
    let maf_file_slice = maf_idx.map(|_| maf_file.as_slice());
    let info_slice = info_idx.map(|_| info.as_slice());
    let out_n = sumstats_qc(
        &a1_ref,
        &a2_ref,
        &a1_file,
        &a2_file,
        &effect,
        &se,
        &p,
        &n_values,
        &maf_ref,
        maf_file_slice,
        info_slice,
        info_filter,
        ols,
        beta_is_character,
        linprob,
        se_logit,
        &mut SumstatsQcOutput {
            keep: &mut keep,
            beta: &mut beta,
            se_out: &mut se_out,
            counts: &mut counts,
        },
    )?;

    for out_idx in 0..out_n {
        let candidate_idx = (keep[out_idx] - 1) as usize;
        out.keep[out_idx] = ref_indices[candidate_idx] as i32;
        out.beta[out_idx] = beta[out_idx];
        out.se[out_idx] = se_out[out_idx];
    }

    report.rows_written = out_n;
    report.counts = counts;
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{allele_flip, median_finite, p_to_z_abs, valid_allele};

    #[test]
    fn allele_labels_match_r_factor_codes() {
        assert_eq!(allele_label(1), "A");
        assert_eq!(allele_label(2), "C");
        assert_eq!(allele_label(3), "G");
        assert_eq!(allele_label(4), "T");
        assert_eq!(allele_label(0), "NA");
    }

    #[test]
    fn native_sumstats_duplicate_counter_matches_drop_all_copies_rule() {
        let mut seen_counts: HashMap<String, usize> = HashMap::new();
        let mut removed = 0usize;
        for snp in ["rs1", "rs2", "rs1", "rs1"] {
            let count = seen_counts.entry(snp.to_owned()).or_insert(0);
            *count += 1;
            if *count == 2 {
                removed += 2;
            } else if *count > 2 {
                removed += 1;
            }
        }

        assert_eq!(removed, 3);
    }

    #[test]
    fn p_to_z_is_finite_for_regular_p_values() {
        let z = p_to_z_abs(0.05);
        assert!(z.is_finite());
        assert!(z > 0.0);
    }

    #[test]
    fn imported_qc_helpers_remain_reachable() {
        assert!(valid_allele(1));
        assert!(allele_flip(1, 2, 1));
        let mut values = vec![0.0, 1.0, 2.0];
        assert_eq!(median_finite(&mut values), 1.0);
    }
}
