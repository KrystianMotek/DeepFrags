use std::clone::Clone;
use std::f64::consts::PI;
use pyo3::prelude::*;

#[pyclass]
#[derive(Clone)]
pub struct Vec3
{
    #[pyo3(get, set)]
    pub x: f64,
    #[pyo3(get, set)]
    pub y: f64,
    #[pyo3(get, set)]
    pub z: f64,
}

#[pyclass]
#[derive(Clone)]
pub struct Output
{
    #[pyo3(get, set)]
    vector: Vec<f64>,
}

#[pymethods]
impl Vec3
{
    #[new]
    pub fn new(x: f64, y: f64, z: f64) -> Self 
    {
        Vec3{x, y, z}
    } 

    pub fn to_list(&self) -> Vec<f64> 
    {
        let list = vec![self.x, self.y, self.z];
        list
    }

    pub fn length(&self) -> f64
    {   
        if self.x == 0.0 && self.y == 0.0 && self.z == 0.0
        {
            return 0.0;
        }
        else
        {
            return (self.x.powf(2.0) + self.y.powf(2.0) + self.z.powf(2.0)).sqrt();
        }
    }

    pub fn normalize(&mut self)
    {
        let length = self.length();
        self.x = self.x / length;
        self.y = self.y / length;
        self.z = self.z / length;
    }
    
    pub fn add(&mut self, vector: &Vec3)
    {
        self.x += vector.x;
        self.y += vector.y;
        self.z += vector.z;
    }

    pub fn subtract(&mut self, vector: &Vec3)
    {
        self.x -= vector.x;
        self.y -= vector.y;
        self.z -= vector.z;
    }
}

#[pymethods]
impl Output
{
    #[new]
    pub fn new(vector: Vec<f64>) -> Self
    {
        Output{vector}
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
}

#[pyfunction]
pub fn to_degrees(radians: f64) -> f64
{
    radians * 180.0 / PI
}

#[pyfunction]
pub fn to_radians(degrees: f64) -> f64
{
    degrees * PI / 180.0
}

#[pyfunction]
pub fn dot_product(a: &Vec3, b: &Vec3) -> f64
{
    a.x * b.x + a.y * b.y + a.z * b.z
}

#[pyfunction]
pub fn cross_product(a: &Vec3, b: &Vec3) -> Vec3
{
    let result_x: f64 = a.y * b.z - a.z * b.y;
    let result_y: f64 = a.z * b.x - a.x * b.z;
    let result_z: f64 = a.x * b.y - a.y * b.x;
    Vec3{x: result_x, y: result_y, z: result_z,}
}

#[pyfunction]
pub fn two_atoms_vector(atom_1: &Vec3, atom_2: &Vec3) -> Vec3
{
    let v_x: f64 = atom_2.x - atom_1.x;
    let v_y: f64 = atom_2.y - atom_1.y;
    let v_z: f64 = atom_2.z - atom_1.z;
    Vec3{x: v_x, y: v_y, z: v_z,}
}

#[pyfunction]
pub fn compute_planar(atom_1: &Vec3, atom_2: &Vec3, atom_3: &Vec3) -> f64
{
    let v_12 = two_atoms_vector(&atom_1, &atom_2);
    let v_23 = two_atoms_vector(&atom_2, &atom_3);
    let planar: f64 = (dot_product(&v_12, &v_23) / (v_12.length() * v_23.length())).acos();
    to_degrees(planar)
}

#[pyfunction]
pub fn atan2(y: f64, x: f64) -> f64
{
    if x > 0.0
    {
        return (y / x).atan();
    }
    if x < 0.0 && y >= 0.0
    {
        return (y / x).atan() + PI;
    }
    if x < 0.0 && y < 0.0
    {
        return (y / x).atan() - PI;
    }
    if x == 0.0 && y > 0.0
    {
        return PI / 2.0;
    }
    if x == 0.0 && y < 0.0
    {
        return (-1.0) * PI / 2.0;
    }
    else
    {
        return 999.9;
    }
}

#[pyfunction]
pub fn compute_dihedral(atom_1: &Vec3, atom_2: &Vec3, atom_3: &Vec3, atom_4: &Vec3) -> f64
{
    let mut v_12 = two_atoms_vector(&atom_1, &atom_2);
    let mut v_23 = two_atoms_vector(&atom_2, &atom_3);
    let mut v_34 = two_atoms_vector(&atom_3, &atom_4);
    
    v_12.normalize();
    v_23.normalize();
    v_34.normalize();

    // vectors orthogonal to both planes
    let k = cross_product(&v_12, &v_23);
    let l = cross_product(&v_23, &v_34);
    let x: f64 = dot_product(&k, &l);
    let y: f64 = dot_product(&v_12, &l);
    let dihedral: f64 = atan2(y, x);
    to_degrees(dihedral)
}

#[pyfunction]
pub fn sin_cos_to_angle(sin: f64, cos: f64) -> f64
{
    let k: f64 = (sin.powf(2.0) + cos.powf(2.0)).sqrt();
    let angle: f64 = atan2(sin / k, cos / k);
    to_degrees(angle)
}

#[pyfunction]
pub fn angles_to_cartesian(atom_1: &Vec3, atom_2: &Vec3, atom_3: &Vec3, bond_length: f64, alpha: f64, theta: f64) -> Vec3
{
    let sin_alpha: f64 = to_radians(alpha).sin();
    let cos_alpha: f64 = to_radians(alpha).cos();
    let sin_theta: f64 = to_radians(theta).sin();
    let cos_theta: f64 = to_radians(theta).cos();
    
    let x: f64 = bond_length * cos_alpha;
    let y: f64 = bond_length * sin_alpha * cos_theta;
    let z: f64 = bond_length * sin_alpha * sin_theta;

    let v_12 = two_atoms_vector(&atom_1, &atom_2);
    let mut v_23 = two_atoms_vector(&atom_2, &atom_3);
    v_23.normalize();

    let mut k = cross_product(&v_12, &v_23);
    k.normalize();

    let l = cross_product(&k, &v_23);

    let new_x: f64 = atom_3.x - v_23.x * x + l.x * y + k.x * z;
    let new_y: f64 = atom_3.y - v_23.y * x + l.y * y + k.y * z;
    let new_z: f64 = atom_3.z - v_23.z * x + l.z * y + k.z * z;
    Vec3{x: new_x, y: new_y, z: new_z} // displacement from the last atom
}

#[pyfunction]
pub fn build_fragment(c_1: Vec3, c_2: Vec3, c_3: Vec3, output: Output, bond_length: f64) -> Vec<Vec3>
{
    let n = output.n();
    let alpha = output.alpha();
    let theta = output.theta();

    let start: usize = 3;
    let end: usize = n + 3;

    let mut atoms = vec![c_1, c_2, c_3];
    for i in start..end
    {
        let mut c_i = atoms[i-3].clone();
        let mut c_j = atoms[i-2].clone();
        let mut c_k = atoms[i-1].clone();

        let mut c_new = angles_to_cartesian(&mut c_i, &mut c_j, &mut c_k, bond_length, alpha[i-3], theta[i-3]);
        c_new.add(&c_k);
        atoms.push(c_new);
    }
    atoms
}

#[pyfunction]
pub fn compute_rmsd(a: Vec<Vec3>, b: Vec<Vec3>) -> f64
{
    let mut square_displacements: Vec<f64> = vec![];
    for iterator in a.iter().zip(b.iter())
    {
        let (vector_a, vector_b) = iterator;
        let displacement: f64 = two_atoms_vector(vector_a, vector_b).length();
        square_displacements.push(displacement.powf(2.0));
    }
    
    let mut rmsd: f64 = square_displacements.iter().sum();
    rmsd = (rmsd / square_displacements.len() as f64).sqrt();
    rmsd
}