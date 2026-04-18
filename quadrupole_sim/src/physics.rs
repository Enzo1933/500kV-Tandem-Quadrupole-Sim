#![allow(non_snake_case)]
use ndarray::{Array2, array};

use crate::PROTON_MASS;

/// Calculates the beam rigidity (B_rho)
/// Dimensions: MV/c
/// Parameters: ke_mev [the kinetic energy in MeV]
pub fn beam_rigidity(ke_mev: f32) -> f32 {
    let p = ((ke_mev + PROTON_MASS).powi(2) - PROTON_MASS.powi(2)).sqrt(); // Momentum

    p / 299.792458
}

/// Calculates the quadrupole transfer matrix
pub fn quad_transfer_matrix(
    g: f32,     // Field gradient
    L: f32,     // Effective length
    B_rho: f32, // The beam rigidity
) -> (Array2<f32>, Array2<f32>) {
    let mut k2 = g / B_rho; // The magnetic gradient
    let mut M_f: Array2<f32> = Array2::default((2, 2)); // Focusing matrix
    let mut M_d: Array2<f32> = Array2::default((2, 2)); // Defocusing matrix
    let k = k2.abs().sqrt();

    if k2.abs() < 1e-9 {
        // Near-zero gradient: both planes are drifts
        let drift = drift_matrix(L);
        return (drift.clone(), drift);
    }

    if k2 > 0.0 {
        // Focusing in x, Defocusing in y
        M_f = array![
            [(L * k.sqrt()).cos(), (L * k.sqrt()).sin() / k.sqrt()],
            [-1.0 * (L * k.sqrt()).sin() * k.sqrt(), (L * k.sqrt()).cos()]
        ];
        M_d = array![
            [(L * k.sqrt()).cosh(), (L * k.sqrt()).sinh() / k.sqrt()],
            [
                -1.0 * (L * k.sqrt()).sinh() * k.sqrt(),
                (L * k.sqrt()).cosh()
            ]
        ];
    } else {
        // Defocusing in x, Focusing in y
        M_f = array![
            [(L * k.sqrt()).cosh(), (L * k.sqrt()).sinh() / k.sqrt()],
            [
                -1.0 * (L * k.sqrt()).sinh() * k.sqrt(),
                (L * k.sqrt()).cosh()
            ]
        ];
        M_d = array![
            [(L * k.sqrt()).cos(), (L * k.sqrt()).sin() / k.sqrt()],
            [-1.0 * (L * k.sqrt()).sin() * k.sqrt(), (L * k.sqrt()).cos()]
        ];
    }

    (M_f, M_d)
}

/// Returns a drift matrix
pub fn drift_matrix(L: f32) -> Array2<f32> {
    array![[1.0,L], [0.0,1.0]]
}


