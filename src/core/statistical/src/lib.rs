use pyo3::prelude::*;
use structural::{Vec3, sin_cos_to_angle, to_radians, angles_to_cartesian};

#[pyclass]
pub struct Output
{
    #[pyo3(get, set)]
    vector: Vec<f64>,
    #[pyo3(get, set)]
    bond_length: f64,
}

#[pymethods]
impl Output
{
    #[new]
    pub fn new(vector: Vec<f64>, bond_length: f64) -> Self
    {
        Output{vector, bond_length}
    }

    pub fn n(&self) -> usize
    {
        (self.vector.len() / 3).try_into().unwrap()
    }

    pub fn alpha(&self) -> Vec<f64>
    {
        let start: usize = 0;
        let end: usize = self.n();

        let mut alpha: Vec<f64> = vec![];
        for i in start..end
        {
            alpha.push(self.vector[i] * 180.0); // because original data are normalized
        }
        alpha
    }

    pub fn sin_theta(&self) -> Vec<f64>
    {
        let start: usize = self.n();
        let end: usize = 2 * self.n();

        let mut sin_theta: Vec<f64> = vec![];
        for i in start..end
        {
            sin_theta.push(self.vector[i]);
        }
        sin_theta
    }

    pub fn cos_theta(&self) -> Vec<f64>
    {
        let start: usize = 2 * self.n();
        let end: usize = self.vector.len();

        let mut cos_theta: Vec<f64> = vec![];
        for i in start..end
        {
            cos_theta.push(self.vector[i]);
        }
        cos_theta
    }

    pub fn theta(&self) -> Vec<f64>
    {
        let sin_theta = self.sin_theta();
        let cos_theta = self.cos_theta();

        let mut theta: Vec<f64> = vec![];
        for i in 0..self.n()
        {
            let sin: f64 = sin_theta[i];
            let cos: f64 = cos_theta[i];
            let angle: f64 = sin_cos_to_angle(sin, cos);
            theta.push(angle);
        }
        theta
    }

    pub fn to_cartesian(&self) -> Vec<Vec3>
    {
        let alpha = self.alpha();
        let theta = self.theta();

        let x_1: f64 = self.bond_length * (1.0 + to_radians(180.0 - alpha[0]).cos());
        let y_1: f64 = self.bond_length * to_radians(180.0 - alpha[0]).sin();
        let z_1: f64 = 0.0;

        // overhang at the beggining 
        let c_a = Vec3{x: 0.0, y: 0.0, z: 0.0};
        let c_b = Vec3{x: 3.8, y: 0.0, z: 0.0};

        let c_1 = Vec3{x: x_1, y: y_1, z: z_1}; // first proper atom 

        let start: usize = 3;
        let end: usize = self.n() + 2;

        let mut atoms = vec![c_a, c_b, c_1];
        for i in start..end
        {
            let mut c_i = atoms[i-3].clone();
            let mut c_j = atoms[i-2].clone();
            let mut c_k = atoms[i-1].clone();

            let mut c_new = angles_to_cartesian(&mut c_i, &mut c_j, &mut c_k, self.bond_length, alpha[i-2], theta[i-2]);
            c_new.add(&c_k);
            atoms.push(c_new);
        }
        atoms
    }

    pub fn compute_r1n(&self) -> f64
    {
        let n = self.n();
        let coordinates = self.to_cartesian();
        let mut last_atom = coordinates[n+1].clone();
        last_atom.subtract(&coordinates[2]);
        last_atom.length()
    }

    pub fn to_pdb(&self) -> Vec<String>
    {
        let coordinates = self.to_cartesian();

        let start: usize = 2;
        let end: usize = self.n() + 2;

        let mut lines = Vec::new();
        for i in start..end
        {
            let atom = &coordinates[i];
            let index = i - 1;
            let x = atom.to_list()[0];
            let y = atom.to_list()[1];
            let z = atom.to_list()[2];
            let line = format!("ATOM {:>6} CA XXX X {:>16.3} {:>8.3} {:>8.3}", index.to_string(), x, y, z);
            lines.push(line);
        }
        lines
    }
}

#[pymodule]
fn statistical(_: Python, m: &PyModule) -> PyResult<()>
{
    m.add_class::<Output>()?;
    Ok(())
}