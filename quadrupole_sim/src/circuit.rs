use std::f64::consts::{PI, SQRT_2};

use crate::MU0;

pub struct MagnetGeometry {
    pub r_gap: f64,  // Radius of the beam pipe (m)
    pub l_mag: f64,  // Physical length of the magnet (m)
    pub w_pole: f64, // Width of the pole tip face (m)
    pub h_pole: f64, // Height of the pole piece (m)
    pub l_iron: f64, // Average path length through the iron yoke (m)
    pub a_iron: f64, // Average cross-sectional area of the iron (m^2)
}

/// Calculates the three critical reluctances for the matrix
impl MagnetGeometry {
    pub fn calculate_reluctances(&self, mu_r_iron: f64) -> (f64, f64, f64) {
        // Gap Reluctance (Air)
        let a_pole = self.l_mag * self.w_pole;
        let r_gap = self.r_gap / (MU0 * 1.0 * a_pole);

        // Iron Reluctance (Dynamic Material), I will change this to also use other materials later
        let r_iron = self.l_iron / (MU0 * mu_r_iron * self.a_iron);

        // Leakage Reluctance (Air path between adjacent poles)
        // Distance between adjacent poles is roughly sqrt(2) * gap_radius
        let l_leak = SQRT_2 * self.r_gap;
        // Area of the side of the pole leaking to its neighbor
        let a_leak = self.l_mag * self.h_pole;
        let r_leak = l_leak / (MU0 * 1.0 * a_leak);

        (r_gap, r_iron, r_leak)
    }
}
