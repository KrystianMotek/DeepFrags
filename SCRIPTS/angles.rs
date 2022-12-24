pub const PI: f64 = 3.14;

// Vec3 is a simple representation of vector in three dimensional coordinate system
// it is based on BioShell class Vec3I  

pub struct Vec3
{
    x: f64,
    y: f64,
    z: f64,
}

impl Vec3
{
    pub fn length(&self) -> f64
    {
        (self.x.powf(2.0) + self.y.powf(2.0) + self.z.powf(2.0)).sqrt()
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

pub fn z_matrix_to_cartesian(atom_1: &mut Vec3, atom_2: &mut Vec3, atom_3: &mut Vec3, bond_length: f64, planar_angle: f64, dihedral_angle: f64, mut new_atom: &mut Vec3)
{
    let sin_planar: f64 = (PI - planar_angle).sin();
    let cos_planar: f64 = (PI - planar_angle).cos();
    let sin_dih: f64 = dihedral_angle.sin();
    let cos_dih: f64 = dihedral_angle.cos();
    let bs: f64 = bond_length * sin_planar;
    new_atom.x = bond_length * cos_planar;
    new_atom.y = bs * cos_dih;
    new_atom.z = bs * sin_dih;

    // translate atoms
    atom_3.subtract(atom_1);
    atom_2.subtract(atom_1);

    let d2_xz: f64 = atom_2.x.powf(2.0) + atom_2.z.powf(2.0);
    let d2: f64 = (d2_xz + atom_2.y.powf(2.0)).sqrt();
    let dxz = d2_xz.sqrt();

    let (d2_inverse, dxz_inverse): (f64, f64);
    let (xx1, x2o, y2o, z2o, xz2o): (f64, f64, f64, f64, f64);

    if d2 < 0.001
    {
        d2_inverse = 1.0 / 0.001;
    }
    else
    {
        d2_inverse = 1.0 / d2;
    }

    if dxz < 0.001
    {
        xx1 = atom_3.x;
        x2o = 1.0;
        z2o = 0.0;
    }
    else
    {
        dxz_inverse = 1.0 / dxz;
        x2o = atom_2.x * dxz_inverse;
        z2o = atom_2.z * dxz_inverse;
        xx1 = atom_3.x * x2o + atom_3.z * z2o;
        atom_3.z = atom_3.z * x2o - atom_3.x * z2o;
    }

    xz2o = dxz * d2_inverse;
    y2o = atom_2.y * d2_inverse;
    atom_3.x = xx1 * (-1.0) * xz2o - atom_3.y * y2o;
    atom_3.y = xx1 * y2o - atom_3.y * xz2o;

    let (dyz, dyz_inverse): (f64, f64);
    let (y1o, z1o, yy4, xx4, zz4): (f64, f64, f64, f64, f64);

    dyz = (atom_3.y.powf(2.0) + atom_3.z.powf(2.0)).sqrt();
    dyz_inverse = 1.0 / dyz;
    y1o = atom_3.y * dyz_inverse;
    z1o = atom_3.z * dyz_inverse;
    yy4 = y1o * new_atom.y - z1o * new_atom.z;
    zz4 = y1o * new_atom.z + z1o * new_atom.y;
    xx4 = y2o * yy4 - xz2o * new_atom.x;

    new_atom.y = -xz2o * yy4 - y2o * new_atom.x;
    new_atom.x = x2o * xx4 - z2o * zz4;
    new_atom.z = z2o * xx4 + x2o * zz4;

    new_atom.add(atom_1); // finally update position of joined atom 
}

fn main()
{   
    // Test
    let r: f64 = 3.8;
    let alpha: f64 = 85.0;
    let theta: f64 = 160.0;
    let mut _c1 = Vec3{x: 0.0, y: 0.0, z: 0.0};
    let mut _c2 = Vec3{x: 3.8, y: 0.0, z: 0.0};
    let mut _c3 = Vec3{x: 3.8, y: 3.8, z: 0.0};
    let mut _c4 = Vec3{x: 0.0, y: 0.0, z: 0.0}; // coordinates of new atom can be random during initialization 

    z_matrix_to_cartesian(&mut _c1, &mut _c2, &mut _c3, r, alpha.to_radians(), theta.to_radians(), &mut _c4);
    println!("x {:.3} y {:.3} z {:.3}", _c4.x, _c4.y, _c4.z);
    println!("length {:.3}", _c4.length())
}
