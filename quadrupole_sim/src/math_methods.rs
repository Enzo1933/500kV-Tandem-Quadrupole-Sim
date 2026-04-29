use nalgebra::{Matrix2xX, SMatrix, Vector1, Vector2, matrix, vector};

use crate::{
    beam_and_tracker::{Beam, Tracker},
    magnet::MagnetGeometry,
};

pub fn x_prime(state: Vector2<f64>, k: f64) -> Vector2<f64> {
    vector![
        state[1],      // x'
        -k * state[0], // x''
    ]
}

pub fn y_prime(state: Vector2<f64>, k: f64) -> Vector2<f64> {
    vector![
        state[1],     // y'
        k * state[0], // y''
    ]
}

pub fn rk4_step<F>(state: Vector2<f64>, z: f64, dz: f64, f: F) -> Vector2<f64>
where
    F: Fn(f64, Vector2<f64>) -> Vector2<f64>,
{
    let k1 = f(z, state);
    let k2 = f(z + 0.5 * dz, state + 0.5 * dz * k1);
    let k3 = f(z + 0.5 * dz, state + 0.5 * dz * k2);
    let k4 = f(z + dz, state + dz * k3);

    state + (dz / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
}

pub fn get_residuals_from_mmf(mmf1: f64, mmf2: f64, beam: &Beam, geo: &MagnetGeometry) -> Vector2<f64> {
    let g1 = geo.field_gradient(mmf1);
    let g2 = geo.field_gradient(mmf2);
    let target_spot = beam.x0 * 0.1;

    let t = Tracker::new(beam, geo, g1, g2, 75).unwrap();
    // residual 0: x/y asymmetry        → drives mmf1/mmf2 ratio
    // residual 1: average spot size    → drives overall MMF scale
    let avg = (t.x_f.abs() + t.y_f.abs()) / 2.0;

    vector![(t.x_f - t.y_f), (avg - target_spot)]
}

pub fn jacobian(mmf1: f64, mmf2: f64, beam: &Beam, geo: &MagnetGeometry) -> SMatrix<f64, 2, 2> {
    let eps = 50.0;

    let r = get_residuals_from_mmf(mmf1, mmf2, beam, geo);
    let r1 = get_residuals_from_mmf(mmf1 + eps, mmf2, beam, geo);
    let r2 = get_residuals_from_mmf(mmf1, mmf2 + eps, beam, geo);

    matrix![
        (r1[0] - r[0]) / eps, (r2[0] - r[0]) / eps;
        (r1[1] - r[1]) / eps, (r2[1] - r[1]) / eps;
    ]
}
