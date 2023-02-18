use pyo3::prelude::*;
use structural::{Vec3};

#[pyclass]
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
}

#[pymodule]
fn insertion(_: Python, m: &PyModule) -> PyResult<()>
{
    m.add_class::<Fragment>()?;
    Ok(())
}