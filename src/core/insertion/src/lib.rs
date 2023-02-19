use pyo3::prelude::*;
use structural::{Vec3};

#[pyclass]
#[derive(Clone)]
pub struct Fragment
{
    #[pyo3(get, set)]
    pub atoms: Vec<Vec3>,
}

#[pymethods]
impl Fragment
{
    #[new]
    pub fn new(atoms: Vec<Vec3>) -> Self
    {
        Fragment{atoms}
    }

    pub fn to_local(&mut self)
    {
        let n: usize = self.atoms.len();
        for i in 0..(n-1)
        {
            let start = self.atoms[0].clone();
            self.atoms[i].subtract(&start);
        }
    }

    pub fn spanning_vector(&self) -> Vec3
    {
        // from the first atom to the last one
        let n: usize = self.atoms.len();
        let start = self.atoms[0].clone();
        let mut end = self.atoms[n-1].clone();
        end.subtract(&start);
        let vector = end.clone();
        vector
    }
}

#[pymodule]
fn insertion(_: Python, m: &PyModule) -> PyResult<()>
{
    m.add_class::<Fragment>()?;
    Ok(())
}