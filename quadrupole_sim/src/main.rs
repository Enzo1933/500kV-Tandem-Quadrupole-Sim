#![allow(unused, non_snake_case)]
mod physics;

use std::io::{read_to_string, stdin};

use crate::physics::{Beam, Tracker};

/// Constants in terms of Natural Units
const PROTON_MASS: f64 = 938.7; // The mass of a proton in MeV/c^2

/// SI Units and Conversions
const MU0: f64 = 1.256_637_061_4e-6; // Magnetic permeability T*m / A
const C_TM: f64 = 299.792_458; // Speed of light c = 1 in Natural Units; conversion factor might be needed: MeV/T*m

const IN_TO_M: f64 = 0.0254; // Conversion factor 
const MM_TO_M: f64 = 1e-3; // Conversion factor

fn main() -> std::io::Result<()> {
    let mut buffer: String = String::new();

    // let split = buffer.split_whitespace();
    // L_mag_m: f64
    // gap_m: f64
    // drift_m: f64
    // energy_MeV: f64
    // x0: f64
    // xp0: f64

    println!(
        "Please input the entire geometry, currents, and material properties with spaces in between each new criteria

                L_mag_m: f64
                gap_m: f64
                drift_m: f64
                energy_MeV: f64
                x0: f64
                xp0: f64
                n_turns: usize
                mu_r: f64
                bore: f64
                i1: f64
                i2: f64
    "
    );
    println!("Enter here: ");
    let buffer_2 = stdin().read_line(&mut buffer)?;
    let split = buffer
        .split_whitespace()
        .collect::<Vec<&str>>()
        .iter()
        .map(|x| x.parse::<f64>().unwrap())
        .collect::<Vec<f64>>();

    let L_mag_m: f64 = split[0] * IN_TO_M;
    let gap_m: f64 = split[1] * IN_TO_M;
    let drift_m: f64 = split[2] * IN_TO_M;
    let r = split[8] * IN_TO_M;
    let x0: f64 = split[4] * IN_TO_M;
    let xp0: f64 = split[5] * IN_TO_M;
    
    let energy_MeV: f64 = split[3];
    let n_turns = split[6] as usize;
    let mu_r = split[7];
    let i1 = split[9];
    let i2 = split[10];

    let beam = Beam::new(L_mag_m, gap_m, drift_m, energy_MeV, x0, xp0);
    let export_femm = Tracker::export_femm_lookup(&beam, i1, i2, n_turns, mu_r, r);
    let export_ibsimu = Tracker::export_to_ibsimu(&beam, i1, i2, n_turns, mu_r, r);

    Ok(())
}
