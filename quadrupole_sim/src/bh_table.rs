pub const B_VALUES: [f64; 12] = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.5, 1.6, 1.7, 1.8];

pub const H_LOW_CARBON_STEEL: [f64; 12] = [
    0.0, 20.0, 40.0, 80.0, 150.0, 300.0, 800.0, 2000.0, 4000.0, 8000.0, 16000.0, 30000.0,
];

pub const H_SILICON_STEEL: [f64; 12] = [
    0.0, 10.0, 25.0, 50.0, 90.0, 180.0, 400.0, 1000.0, 2500.0, 6000.0, 15000.0, 35000.0,
];

pub const H_SOFT_IRON: [f64; 12] = [
    0.0, 5.0, 15.0, 40.0, 100.0, 250.0, 700.0, 1800.0, 3500.0, 7000.0, 15000.0, 30000.0,
];

pub const H_GENERIC_STEEL: [f64; 12] = [
    0.0, 15.0, 35.0, 70.0, 130.0, 250.0, 600.0, 1500.0, 3000.0, 7000.0, 15000.0, 30000.0,
];

pub struct BHCurve {
    b: &'static [f64],
    h: &'static [f64],
    m: Vec<f64>, // precomputed slopes (PCHIP)
}

/// This might be short lived, I don't know
impl BHCurve {
    /// Create curve and precompute PCHIP slopes
    pub fn new(b: &'static [f64], h: &'static [f64]) -> Self {
        assert!(b.len() == h.len());
        assert!(b.len() >= 2);

        let n = b.len();
        let mut m = vec![0.0; n];

        // Secant slopes
        let mut delta = vec![0.0; n - 1];
        for i in 0..n - 1 {
            let db = b[i + 1] - b[i];
            assert!(db > 0.0, "B must be strictly increasing");
            delta[i] = (h[i + 1] - h[i]) / db;
        }

        // Endpoints (simple one-sided)
        m[0] = delta[0];
        m[n - 1] = delta[n - 2];

        // Interior points (PCHIP harmonic mean)
        for i in 1..n - 1 {
            if delta[i - 1] * delta[i] <= 0.0 {
                m[i] = 0.0;
            } else {
                m[i] = 2.0 * delta[i - 1] * delta[i] / (delta[i - 1] + delta[i]);
            }
        }

        Self { b, h, m }
    }

    /// Calculates H and slope using Hermite Interpolation
    pub fn h_and_slope(&self, b: f64, mut idx: usize) -> (f64, f64) {
        let n = self.b.len();

        // Clamp low
        if b <= self.b[0] {
            return (self.h[0], self.m[0]);
        }

        // Clamp high
        if b >= self.b[n - 1] {
            return (self.h[n - 1], self.m[n - 1]);
        }

        // Try cached index first
        if idx >= n - 1 {
            idx = n - 2;
        }

        if !(self.b[idx] <= b && b <= self.b[idx + 1]) {
            // Binary search fallback
            idx = match self
                .b
                .binary_search_by(|probe| probe.partial_cmp(&b).unwrap())
            {
                Ok(i) => i.min(n - 2),
                Err(i) => i - 1,
            };
        }

        let b0 = self.b[idx];
        let b1 = self.b[idx + 1];
        let h0 = self.h[idx];
        let h1 = self.h[idx + 1];
        let m0 = self.m[idx];
        let m1 = self.m[idx + 1];

        let dx = b1 - b0;
        let t = (b - b0) / dx;

        let t2 = t * t;
        let t3 = t2 * t;

        // Hermite basis
        let h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
        let h10 = t3 - 2.0 * t2 + t;
        let h01 = -2.0 * t3 + 3.0 * t2;
        let h11 = t3 - t2;

        // Interpolated H
        let H = h00 * h0 + h10 * dx * m0 + h01 * h1 + h11 * dx * m1;

        // Derivative dH/dt
        let dh_dt = (6.0 * t2 - 6.0 * t) * h0
            + (3.0 * t2 - 4.0 * t + 1.0) * dx * m0
            + (-6.0 * t2 + 6.0 * t) * h1
            + (3.0 * t2 - 2.0 * t) * dx * m1;

        // Chain rule
        let dH_dB = dh_dt / dx;

        (H, dH_dB)
    }
}
