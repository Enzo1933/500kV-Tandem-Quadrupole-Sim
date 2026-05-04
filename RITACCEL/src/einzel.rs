use nalgebra::{Matrix2, SMatrix};

use crate::math_methods::sech;

pub struct EinzelGeometry {
    pub U_outer: f64, // Voltage of the outer electrodes
    pub U_mid: f64,   // Sitting voltage of the middle electrode
    pub L_mid: f64,   // Length of the middle electrode
    pub R: f64,       // Radius of the three cylinders
}

impl EinzelGeometry {
    /// Constructs a new einzel lens
    pub fn new(U_outer: f64, U_mid: f64, L_mid: f64, R: f64) -> Self {
        Self {
            U_outer,
            U_mid,
            L_mid,
            R,
        }
    }

    /// Returns the potential (V) at a certain z value
    pub fn voltage(&self, z: f64) -> f64 {
        let k = 1.318 / self.R;
        let u1 = k * (z + 0.5 * self.L_mid);
        let u2 = k * (z - 0.5 * self.L_mid);

        let dU = self.U_mid - self.U_outer;

        self.U_mid + (dU / 2.0) * (u1.tanh() - u2.tanh())
    }

    /// Returns the E field (-V') at a certain z value
    pub fn e_field(&self, z: f64) -> f64 {
        let k = 1.318 / self.R;
        let u1 = k * (z + 0.5 * self.L_mid);
        let u2 = k * (z - 0.5 * self.L_mid);

        let dU = self.U_mid - self.U_outer;

        (-dU * k / 2.0) * (sech(u1).powi(2) - sech(u2).powi(2))
    }

    /// Returns the radial focusing gradient (V'') at a certain z value
    pub fn _rad_focusing_grad(&self, z: f64) -> f64 {
        let k = 1.318 / self.R;
        let u1 = k * (z + 0.5 * self.L_mid);
        let u2 = k * (z - 0.5 * self.L_mid);

        let dU = self.U_mid - self.U_outer;

        (-dU * k.powi(2)) * (sech(u1) * u1.tanh() - sech(u2) * u2.tanh())
    }

    /// Calculates the Einzel Lense Transfer Matrix
    pub fn transfer_matrix(omega: f64, dz: f64) -> SMatrix<f64, 2, 2> {
        if omega < 1e-10 {
            return Matrix2::new(1.0, dz, 0.0, 1.0);
        }

        // Otherwise, apply the standard harmonic oscillator matrix
        Matrix2::new(
            (omega * dz).cos(),
            (1.0 / omega) * (omega * dz).sin(),
            -omega * (omega * dz).sin(),
            (omega * dz).cos(),
        )
    }
}
