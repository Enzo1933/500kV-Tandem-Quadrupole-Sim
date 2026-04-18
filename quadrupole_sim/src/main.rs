mod physics;

/// Constants in terms of Natural Units
const PROTON_MASS: f32 = 938.7; // The mass of a proton in MeV/c^2

/// SI Units and Conversions
const MU0: f32 = 1.256_637_061_4e-6; // Magnetic permeability T*m / A
const C_TM: f32 = 299.792_458; // Speed of light c = 1 in Natural Units; conversion factor might be needed: MeV/T*m

const IN_TO_M: f32 = 0.0254; // Conversion factor 
const MM_TO_M: f32 = 1e-3; // Conversion factor

fn main() -> std::io::Result<()> {
    Ok(())
}
