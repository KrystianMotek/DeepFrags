use std::clone::Clone;
use std::f64::consts::PI;
use pyo3::prelude::*;

#[pyclass]
#[derive(Clone)]
pub struct Vec3
{
    // representation of vector in three dimensional coordinates system
    #[pyo3(get, set)]
    pub x: f64,
    #[pyo3(get, set)]
    pub y: f64,
    #[pyo3(get, set)]
    pub z: f64,
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
pub fn compute_alpha(atom_1: &Vec3, atom_2: &Vec3, atom_3: &Vec3) -> f64
{
    let v_12 = two_atoms_vector(&atom_1, &atom_2);
    let v_23 = two_atoms_vector(&atom_2, &atom_3);
    let alpha: f64 = (dot_product(&v_12, &v_23) / (v_12.length() * v_23.length())).acos();
    to_degrees(alpha)
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
pub fn compute_theta(atom_1: &Vec3, atom_2: &Vec3, atom_3: &Vec3, atom_4: &Vec3) -> f64
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
    let theta: f64 = atan2(y, x);
    to_degrees(theta)
}

#[pyfunction]
pub fn sin_cos_to_angle(sin: f64, cos: f64) -> f64
{
    let k: f64 = (sin.powf(2.0) + cos.powf(2.0)).sqrt();
    let angle: f64 = atan2(sin / k, cos / k);
    to_degrees(angle)
}

#[pyfunction]
pub fn angles_to_cartesian(atom_1: &mut Vec3, atom_2: &mut Vec3, atom_3: &mut Vec3, bond_length: f64, alpha: f64, theta: f64) -> Vec3
{
    let sin_alpha: f64 = to_radians(180.0 - alpha).sin();
    let cos_alpha: f64 = to_radians(180.0 - alpha).cos();
    let sin_theta: f64 = to_radians(theta).sin();
    let cos_theta: f64 = to_radians(theta).cos();
    
    let mut new_x: f64 = bond_length * cos_alpha;
    let mut new_y: f64 = bond_length * sin_alpha * cos_theta;
    let mut new_z: f64 = bond_length * sin_alpha * sin_theta;

    // transform to local coordinates
    atom_1.subtract(atom_3);
    atom_2.subtract(atom_3);

    let d_2_xz: f64 = (atom_2.x.powf(2.0) + atom_2.z.powf(2.0)).sqrt();
    let d_2: f64 = atom_2.length();
    
    let (d_2_inverse, d_2_xz_inverse): (f64, f64);
    let (xx_3, x_2_o, y_2_o, z_2_o, xz_2_o): (f64, f64, f64, f64, f64);

    if d_2 < 0.001
    {
        d_2_inverse = 1.0 / 0.001;
    }
    else
    {
        d_2_inverse = 1.0 / d_2;
    }

    if d_2_xz < 0.001
    {
        xx_3 = atom_1.x;
        x_2_o = 1.0;
        z_2_o = 0.0;
    }
    else
    {
        d_2_xz_inverse = 1.0 / d_2_xz;
        x_2_o = atom_2.x * d_2_xz_inverse;
        z_2_o = atom_2.z * d_2_xz_inverse;
        xx_3 = atom_1.x * x_2_o + atom_1.z * z_2_o;
        atom_1.z = atom_1.z * x_2_o - atom_1.x * z_2_o;
    }

    xz_2_o = d_2_xz * d_2_inverse;
    y_2_o = atom_2.y * d_2_inverse;
    atom_1.x = (-1.0) * (xx_3 * xz_2_o - atom_1.y * y_2_o);
    atom_1.y = xx_3 * y_2_o - atom_1.y * xz_2_o;

    let (d_1_yz, d_1_yz_inverse): (f64, f64);
    let (y_3_o, z_3_o, yy_4, xx_4, zz_4): (f64, f64, f64, f64, f64);

    d_1_yz = (atom_1.y.powf(2.0) + atom_1.z.powf(2.0)).sqrt();
    d_1_yz_inverse = 1.0 / d_1_yz;
    y_3_o = atom_1.y * d_1_yz_inverse;
    z_3_o = atom_1.z * d_1_yz_inverse;
    yy_4 = y_3_o * new_y - z_3_o * new_z;
    zz_4 = y_3_o * new_z + z_3_o * new_y;
    xx_4 = y_2_o * yy_4 - xz_2_o * new_x;

    new_y = (-1.0) * (xz_2_o * yy_4 + y_2_o * new_x);
    new_x = x_2_o * xx_4 - z_2_o * zz_4;
    new_z = z_2_o * xx_4 + x_2_o * zz_4;

    // displacement from the last atom
    let new_atom = Vec3{x: new_x, y: new_y, z: new_z};
    new_atom
}

#[pyfunction]
pub fn euler_angles(row_x: Vec<f64>, row_y: Vec<f64>, row_z: Vec<f64>) -> Vec<f64>
{
    let angles = vec![(-1.0) * row_z[0].asin(), atan2(row_z[1], row_z[2]), atan2(row_y[0], row_x[0])];
    angles
}

#[pymodule]
fn structural(_: Python, m: &PyModule) -> PyResult<()>
{
    m.add_class::<Vec3>()?;
    m.add_function(wrap_pyfunction!(to_degrees, m)?).unwrap();
    m.add_function(wrap_pyfunction!(to_radians, m)?).unwrap();  
    m.add_function(wrap_pyfunction!(dot_product, m)?).unwrap();
    m.add_function(wrap_pyfunction!(cross_product, m)?).unwrap();
    m.add_function(wrap_pyfunction!(two_atoms_vector, m)?).unwrap();
    m.add_function(wrap_pyfunction!(compute_alpha, m)?).unwrap();
    m.add_function(wrap_pyfunction!(atan2, m)?).unwrap();
    m.add_function(wrap_pyfunction!(compute_theta, m)?).unwrap();
    m.add_function(wrap_pyfunction!(sin_cos_to_angle, m)?).unwrap();
    m.add_function(wrap_pyfunction!(angles_to_cartesian, m)?).unwrap();
    m.add_function(wrap_pyfunction!(euler_angles, m)?).unwrap();
    Ok(())
}