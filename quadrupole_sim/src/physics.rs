use crate::PROTON_MASS;

/// Calculates the beam rigidity (B_rho)
/// Dimensions: MV/c
/// Parameters: ke_mev [the kinetic energy in MeV]
pub fn beam_rigidity(ke_mev: f32) -> f32 {
    let p = f32::sqrt((ke_mev + PROTON_MASS).powi(2) - PROTON_MASS.powi(2)); // Momentum

    p
}
