use pyo3::prelude::*;

#[pyclass]
pub struct Histogram1D
{
    #[pyo3(get, set)]
    data: Vec<f64>,
    #[pyo3(get, set)]
    bins: i32,
    #[pyo3(get, set)]
    min: f64,
    #[pyo3(get, set)]
    max: f64,
}

#[pyclass]
pub struct Histogram2D
{
    #[pyo3(get, set)]
    x: Vec<f64>,
    #[pyo3(get, set)]
    y: Vec<f64>,
    #[pyo3(get, set)]
    bins: i32,
    #[pyo3(get, set)]
    x_min: f64,
    #[pyo3(get, set)]
    x_max: f64,
    #[pyo3(get, set)]
    y_min: f64,
    #[pyo3(get, set)]
    y_max: f64,
}

#[pymethods]
impl Histogram1D
{
    #[new]
    pub fn new(data: Vec<f64>, bins: i32, min: f64, max: f64) -> Self
    {
        Histogram1D{data, bins, min, max}
    }

    pub fn split(&self) -> Vec<(f64, f64, i32)>
    {
        let pointer: f64 = (self.max - self.min) / (self.bins as f64);
        let mut results = Vec::new();
        let mut current_value: f64 = self.min;
        while current_value < self.max
        {
            let mut n: i32 = 0;
            for element in self.data.iter()
            {
                if element > &current_value && element <= &(current_value + pointer)
                {
                    n += 1;
                }
                if current_value == self.min && element == &self.min
                {
                    n += 1;
                }
            }
            results.push((current_value, current_value + pointer, n));
            current_value += pointer;
        }
        results
    }
}

#[pymethods]
impl Histogram2D
{
    #[new]
    pub fn new(x: Vec<f64>, y: Vec<f64>, bins: i32, x_min: f64, x_max: f64, y_min: f64, y_max: f64) -> Self
    {
        Histogram2D{x, y, bins, x_min, x_max, y_min, y_max}
    }

    pub fn split(&self) -> Vec<(f64, f64, f64, f64, i32)>
    {
        let x_pointer: f64 = (self.x_max - self.x_min) / (self.bins as f64);
        let y_pointer: f64 = (self.y_max - self.y_min) / (self.bins as f64);
        let mut results = Vec::new();
        let mut current_x: f64 = self.x_min;
        while current_x < self.x_max 
        {
            let mut current_y: f64 = self.y_min;
            while current_y < self.y_max
            {
                let mut n: i32 = 0;
                for x in self.x.iter()
                {
                    let index = self.x.iter().position(|&k| k == *x).unwrap();
                    if x > &current_x && x <= &(current_x + x_pointer)
                    {
                        if self.y[index] > current_y && self.y[index] <= (current_y + y_pointer)
                        {
                            n += 1;
                        }
                        if current_y == self.y_min && self.y[index] == self.y_min
                        {
                            n += 1;
                        }
                    }
                    if current_x == self.x_min && x == &self.x_min
                    {
                        if self.y[index] > current_y && self.y[index] <= (current_y + y_pointer)
                        {
                            n += 1;
                        }
                        if current_y == self.y_min && self.y[index] == self.y_min
                        {
                            n += 1;
                        }
                    }
                }
                results.push((current_x, current_x + x_pointer, current_y, current_y + y_pointer, n));
                current_y += y_pointer;
            }
            current_x += x_pointer;
        }
        results
    }
}