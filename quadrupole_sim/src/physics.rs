#![allow(non_snake_case)]
use ndarray::{Array2, array};

use crate::PROTON_MASS;

/// Calculates the beam rigidity (B_rho)
/// Dimensions: MV/c
/// Parameters: ke_mev [the kinetic energy in MeV]
pub fn beam_rigidity(ke_mev: f32) -> f32 {
    let p = ((ke_mev + PROTON_MASS).powi(2) - PROTON_MASS.powi(2)).sqrt(); // Momentum

    p
}

/// Calculates the quadrupole transfer matrix
pub fn quad_transfer_matrix(
    g: f32,     // Field gradient
    L: f32,     // Effective length
    B_rho: f32, // The beam rigidity
) -> (Array2<f32>, Array2<f32>) {
    let mut k = g / B_rho; // The magnetic gradient
    let mut M_f: Array2<f32> = Array2::default((2, 2)); // Focusing matrix
    let mut M_d: Array2<f32> = Array2::default((2, 2)); // Defocusing matrix

    if k > 0.0 {
        // Focusing
        M_f = array![
            [(L * k.sqrt()).cos(), (L * k.sqrt()).sin() / k.sqrt()],
            [-1.0 * (L * k.sqrt()).sin() * k.sqrt(), (L * k.sqrt()).cos()]
        ]
    } else {
        // Defocusing
        k = k.abs();
        M_d = array![
            [(L * k.sqrt()).cosh(), (L * k.sqrt()).sinh() / k.sqrt()],
            [
                -1.0 * (L * k.sqrt()).sinh() * k.sqrt(),
                (L * k.sqrt()).cosh()
            ]
        ]
    }

    (M_f, M_d)
}
