use std::io::BufReader;
use std::fs::read_dir;
use std::io::BufRead;
use std::path::Path;
use std::ffi::OsStr;
use std::fs::File;
use pyo3::prelude::*;

// some functionalities for data generating
// extract random fragments from proteins collection

static AA_CODES: &'static [char] = &['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'];
static SS_CODES: &'static [char] = &['H', 'E', 'C'];

#[pyfunction]
pub fn get_extension(file: &str) -> Option<&str>
{
    Path::new(file).extension().and_then(OsStr::to_str)
}

#[pyfunction]
pub fn check_if_correct(line: &str) -> bool
{
    // split line into vector of strings
    let elements: Vec<&str> = line.split_whitespace().collect();

    // get amino acids and secondary structure vectors
    let aa: Vec<char> = elements[4].chars().collect();
    let ss: Vec<char> = elements[5].chars().collect();

    let mut correct = true;
    for aa in aa
    {
        if !AA_CODES.contains(&aa)
        {
            correct = false;
        }
    }
    for ss in ss
    {
        if !SS_CODES.contains(&ss)
        {
            correct = false;
        }
    }
    correct
}

#[pyfunction]
pub fn collect_files(directory: &str) -> Vec<String>
{
    // get path for each file
    let all = read_dir(directory).unwrap();
    let mut files = Vec::new();
    for item in all
    {
        let file = item.unwrap().path().display().to_string();
        if Some(get_extension(&file)) == Some(Some("dat"))
        {
            files.push(file);
        }
    }
    files
}

#[pyfunction]
pub fn read_lines(file: &str) -> Vec<String>
{
    let stream = File::open(Path::new(file)).unwrap();
    let reader = BufReader::new(&stream);
    let all: Vec<String> = reader.lines().collect::<Result<_, _>>().unwrap();
    let mut lines = Vec::new();
    for line in all
    {
        // headlines are removed automatically
        if check_if_correct(&line) == true
        {
            lines.push(line);
        }
    }
    lines
}

#[pyfunction]
pub fn all_samples(directory: &str) -> Vec<String>
{
    // iterate over all files and get lines with data
    let files: Vec<String> = collect_files(directory);
    let mut samples = Vec::new();
    for file in files
    {
        let lines: Vec<String> = read_lines(&file);
        for line in lines
        {
            samples.push(line);
        }
    }
    samples
}

#[pymodule]
fn data(_: Python, m: &PyModule) -> PyResult<()>
{
    m.add_function(wrap_pyfunction!(check_if_correct, m)?).unwrap();
    m.add_function(wrap_pyfunction!(get_extension, m)?).unwrap();
    m.add_function(wrap_pyfunction!(collect_files, m)?).unwrap();
    m.add_function(wrap_pyfunction!(read_lines, m)?).unwrap();
    m.add_function(wrap_pyfunction!(all_samples, m)?).unwrap();
    Ok(())
}