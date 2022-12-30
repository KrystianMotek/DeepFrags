use pyo3::prelude::*;

// some functionalities for data generating
// extract random fragments from proteins collection

#[pyfunction]
pub fn check_if_correct(line: str) -> bool
{
    let ok = true;
    let split = line.split_whitespace();
    let elements: Vec<&str> = split.collect();
}

#[pymodule]
fn scraper(_: Python, m: &PyModule) -> PyResult<()>
{
    Ok(())
}
