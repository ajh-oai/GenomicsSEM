use genomicssem_core::{fill_v_snp, fill_v_snp_batch, GenomicControl};
use std::env;
use std::time::Instant;

fn next_unit(state: &mut u64) -> f64 {
    *state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
    ((*state >> 11) as f64) / ((1_u64 << 53) as f64)
}

fn runif(state: &mut u64, lo: f64, hi: f64) -> f64 {
    lo + (hi - lo) * next_unit(state)
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let n_snp = args
        .get(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(2_000_000);
    let k = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(8);

    let mut rng = 1_u64;
    let se_snp: Vec<f64> = (0..(n_snp * k))
        .map(|_| runif(&mut rng, 0.01, 0.5))
        .collect();
    let var_snp: Vec<f64> = (0..n_snp).map(|_| runif(&mut rng, 0.05, 0.95)).collect();

    let mut i_ld = vec![0.0; k * k];
    for col in 0..k {
        for row in 0..k {
            i_ld[row + col * k] = runif(&mut rng, 0.05, 0.4);
        }
    }
    for col in 0..k {
        for row in 0..k {
            let avg = 0.5 * (i_ld[row + col * k] + i_ld[col + row * k]);
            i_ld[row + col * k] = avg;
            i_ld[col + row * k] = avg;
        }
    }
    for d in 0..k {
        i_ld[d + d * k] = runif(&mut rng, 1.0, 1.8);
    }

    let coords_nrow = k * k;
    let mut coords = vec![0_i32; coords_nrow * 2];
    for col in 0..k {
        for row in 0..k {
            let p = row + col * k;
            coords[p] = (row + 1) as i32;
            coords[p + coords_nrow] = (col + 1) as i32;
        }
    }

    let mut one = vec![0.0; k * k];
    let start = Instant::now();
    let mut checksum = 0.0;
    for i_zero in 0..n_snp {
        fill_v_snp(
            &se_snp,
            n_snp,
            k,
            i_zero,
            &i_ld,
            k,
            k,
            &var_snp,
            &coords,
            coords_nrow,
            2,
            k,
            GenomicControl::Standard,
            &mut one,
        )
        .unwrap();
        checksum += one[0] + one[k * k - 1];
    }
    println!(
        "backend,threads,n_snp,k,elapsed_sec,checksum\nrust_loop,1,{n_snp},{k},{:.6},{:.12}",
        start.elapsed().as_secs_f64(),
        checksum
    );

    for n_threads in [1_usize, 2, 4, 8, 16] {
        let mut out = vec![0.0; n_snp * k * k];
        let start = Instant::now();
        fill_v_snp_batch(
            &se_snp,
            n_snp,
            k,
            &i_ld,
            k,
            k,
            &var_snp,
            &coords,
            coords_nrow,
            2,
            k,
            GenomicControl::Standard,
            n_threads,
            &mut out,
        )
        .unwrap();
        let checksum: f64 = out.chunks_exact(k * k).map(|m| m[0] + m[k * k - 1]).sum();
        println!(
            "rust_batch,{n_threads},{n_snp},{k},{:.6},{:.12}",
            start.elapsed().as_secs_f64(),
            checksum
        );
    }
}
