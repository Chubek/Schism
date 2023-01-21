#![allow(unused)]

use std::borrow::BorrowMut;

#[derive(Clone)]
pub struct ConductanceMatrix {
    n: usize,
    m: usize,
    matrix: Vec<Vec<f64>>
}

impl ConductanceMatrix {
    pub fn new(n: usize, m: usize) -> Self {
        let matrix = vec![vec![0.0f64; m]; n];

        Self { n, m, matrix }
    }

    pub fn set_at(&mut self, n_prime: usize, m_prime: usize, conductance: f64) {
        let row = self.matrix.get_mut(n_prime).expect("Number of rows exceed n");
        let row_col = row.get_mut(m_prime).expect("Number or cols exceeds m");

        *row_col = conductance;
    }

    pub fn get_at(&self, n_prime: usize, m_prime: usize) -> &f64 {
        let row = self.matrix.get(n_prime).expect("Number exceeds n");
        let row_col = row.get(m_prime).expect("Number exceeds M");

        row_col
    }
}

impl std::ops::Index<(usize, usize)> for ConductanceMatrix {
    type Output = f64;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        let (n_prime, m_prime) = index;

        self.get_at(n_prime, m_prime)
    }
}

impl std::ops::IndexMut<(usize, usize)> for ConductanceMatrix {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        let (n_prime, m_prime) = index;

        let row = self.matrix.get_mut(n_prime).expect("Number exceeds N");
        let row_col = row.get_mut(m_prime).expect("Number exceeds M");

        row_col
    }
}