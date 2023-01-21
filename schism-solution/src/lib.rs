#![allow(unused)]

use num_traits::Num;
use rayon::prelude::*;
use std::{borrow::BorrowMut, marker::PhantomData};

pub enum ElementWiseOperation {
    Divide,
    Multiply,
    Modulo,
    Add,
    Subtract,
}

#[derive(Clone)]
pub struct Matrix<T: Clone + Default + Num> {
    n: usize,
    m: usize,
    matrix: Vec<Vec<T>>,
}

impl<T: Clone + Default + Num> Matrix<T> {
    pub fn new(n: usize, m: usize) -> Self {
        let matrix = vec![vec![T::default(); m]; n];

        Self { n, m, matrix }
    }

    pub fn set_at(&mut self, n_prime: usize, m_prime: usize, conductance: T) {
        let row = self
            .matrix
            .get_mut(n_prime)
            .expect("Number of rows exceed n");
        let row_col = row.get_mut(m_prime).expect("Number or cols exceeds m");

        *row_col = conductance;
    }

    pub fn get_at(&self, n_prime: usize, m_prime: usize) -> &T {
        let row = self.matrix.get(n_prime).expect("Number exceeds n");
        let row_col = row.get(m_prime).expect("Number exceeds M");

        row_col
    }

    pub fn from(matrix: Vec<Vec<T>>) -> Self {
        let n = matrix.len();
        let m = matrix.get(0).expect("Matrix empty").len();

        let mut cm = Self::new(n, m);
        cm.matrix = matrix;

        cm
    }

    pub fn from_one_element(n: usize, m: usize, element: T) -> Self {
        let matrix = vec![vec![element; m]; n];

        Self { n, m, matrix }
    }

    pub fn set_one_element(&mut self, element: T) {
        self.matrix.iter_mut().for_each(|v| {
            v.iter_mut().for_each(|e| *e = element.clone());
        });
    }

    fn multiply_row_i_with_this_col(&mut self, i: usize, col: Vec<T>) {
        (0..col.len())
            .into_iter()
            .for_each(|j| self[(i, j)] = col[j].clone());
    }

    pub fn get_col(&self, j: usize) -> Vec<T> {
        self.matrix
            .clone()
            .into_iter()
            .map(|v| v.get(j).expect("Column larger than M").clone())
            .collect()
    }

    pub fn product(&self, rhs: Self) -> Self {
        assert!(
            self.n == rhs.m && self.m == rhs.n,
            "Rows must be equal to the multipliee's columns and vice versa!"
        );

        let mut result = self.clone();

        (0..result.n).into_iter().for_each(|i| {
            (0..result.m).into_iter().for_each(|j| {
                let col = rhs.get_col(j);

                result.multiply_row_i_with_this_col(i, col);
            })
        });

        result
    }

    pub fn element_wise_op(&self, rhs: Self, op: ElementWiseOperation) -> Self {
        assert!(
            self.n == rhs.n && self.m == rhs.m,
            "Numbr of rows and cols must be equal for element-wise operation."
        );

        let mut result = self.clone();

        (0..self.n).into_iter().for_each(|i| {
            (0..self.m).into_iter().for_each(|j| match op {
                ElementWiseOperation::Divide => {
                    result[(i, j)] = result[(i, j)].clone() / rhs[(i, j)].clone()
                }
                ElementWiseOperation::Multiply => {
                    result[(i, j)] = result[(i, j)].clone() * rhs[(i, j)].clone()
                }
                ElementWiseOperation::Modulo => {
                    result[(i, j)] = result[(i, j)].clone() % rhs[(i, j)].clone()
                }
                ElementWiseOperation::Add => {
                    result[(i, j)] = result[(i, j)].clone() + rhs[(i, j)].clone()
                }
                ElementWiseOperation::Subtract => {
                    result[(i, j)] = result[(i, j)].clone() - rhs[(i, j)].clone()
                }
            })
        });

        result
    }

    pub fn element_wise_op_scalar(&self, rhs: T, op: ElementWiseOperation) -> Self {
        let mut result = self.clone();

        (0..self.n).into_iter().for_each(|i| {
            (0..self.m).into_iter().for_each(|j| match op {
                ElementWiseOperation::Divide => {
                    result[(i, j)] = result[(i, j)].clone() / rhs.clone()
                }
                ElementWiseOperation::Multiply => {
                    result[(i, j)] = result[(i, j)].clone() * rhs.clone()
                }
                ElementWiseOperation::Modulo => {
                    result[(i, j)] = result[(i, j)].clone() % rhs.clone()
                }
                ElementWiseOperation::Add => result[(i, j)] = result[(i, j)].clone() + rhs.clone(),
                ElementWiseOperation::Subtract => {
                    result[(i, j)] = result[(i, j)].clone() - rhs.clone()
                }
            })
        });

        result
    }

    pub fn get_slice(&self, rows_range: (usize, usize), cols_range: (usize, usize)) -> Self {
        let (row_start, row_end) = rows_range;
        let (col_start, col_end) = cols_range;

        let matrix = self
            .matrix
            .iter()
            .cloned()
            .enumerate()
            .filter(|(i, _)| *i > row_start && *i < row_end)
            .map(|(i, v)| v[col_start..col_end].to_vec())
            .collect::<Vec<_>>();

        Self::from(matrix)
    }

    pub fn get_row(&self, idx: usize) -> &Vec<T> {
        self.matrix.get(idx).expect("Size exceeds index")
    }

    pub fn get_row_mut(&mut self, idx: usize) -> &mut Vec<T> {
        self.matrix.get_mut(idx).expect("Size exceeds index")
    }

    pub fn swap_rows(&mut self, this_row: usize, with_this_row: usize) {
        self.matrix.swap(this_row, with_this_row);
    }

    pub fn elementwise_on_row(&mut self, row_idx: usize, rhs: Vec<T>, op: ElementWiseOperation) {
        let row = self
            .matrix
            .get_mut(row_idx)
            .expect("Row index for elementwise smaller than expected");

        assert!(
            row.len() == rhs.len(),
            "Length of RHS not the same as matrix vector for row element wise"
        );

        for (lhs, rhs_) in row.iter_mut().zip(rhs) {
            let cln = lhs.clone();

            match op {
                ElementWiseOperation::Divide => *lhs = cln / rhs_,
                ElementWiseOperation::Multiply => *lhs = cln * rhs_,
                ElementWiseOperation::Modulo => *lhs = cln % rhs_,
                ElementWiseOperation::Add => *lhs = cln + rhs_,
                ElementWiseOperation::Subtract => *lhs = cln - rhs_,
            }
        }
    }

    pub fn elementwise_on_row_scalar(&mut self, row_idx: usize, rhs: T, op: ElementWiseOperation) {
        let row = self
            .matrix
            .get_mut(row_idx)
            .expect("Row index for elementwise smaller than expected");

        for lhs in row.iter_mut() {
            let cln_lhs = lhs.clone();
            let cln_rhs = rhs.clone();

            match op {
                ElementWiseOperation::Divide => *lhs = cln_lhs / cln_rhs,
                ElementWiseOperation::Multiply => *lhs = cln_lhs * cln_rhs,
                ElementWiseOperation::Modulo => *lhs = cln_lhs % cln_rhs,
                ElementWiseOperation::Add => *lhs = cln_lhs + cln_rhs,
                ElementWiseOperation::Subtract => *lhs = cln_lhs - cln_rhs,
            }
        }
    }
}

impl<T: Clone + Default + Num> std::ops::Index<usize> for Matrix<T> {
    type Output = Vec<T>;

    fn index(&self, index: usize) -> &Self::Output {
        self.get_row(index)
    }
}

impl<T: Clone + Default + Num> std::ops::IndexMut<usize> for Matrix<T> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        self.get_row_mut(index)
    }
}

impl<T: Clone + Default + Num> std::ops::Index<(usize, usize)> for Matrix<T> {
    type Output = T;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        let (n_prime, m_prime) = index;

        self.get_at(n_prime, m_prime)
    }
}

impl<T: Clone + Default + Num> std::ops::IndexMut<(usize, usize)> for Matrix<T> {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        let (n_prime, m_prime) = index;

        let row = self.matrix.get_mut(n_prime).expect("Number exceeds N");
        let row_col = row.get_mut(m_prime).expect("Number exceeds M");

        row_col
    }
}

impl From<(usize, usize, f64)> for Matrix<f64> {
    fn from(value: (usize, usize, f64)) -> Self {
        let (m, n, element) = value;

        Self::from_one_element(n, m, element)
    }
}

impl std::ops::Mul<Matrix<f64>> for Matrix<f64> {
    type Output = Matrix<f64>;

    fn mul(self, rhs: Matrix<f64>) -> Self::Output {
        self.element_wise_op(rhs, ElementWiseOperation::Multiply)
    }
}

impl std::ops::Div<Matrix<f64>> for Matrix<f64> {
    type Output = Matrix<f64>;

    fn div(self, rhs: Matrix<f64>) -> Self::Output {
        self.element_wise_op(rhs, ElementWiseOperation::Divide)
    }
}

impl std::ops::Add<Matrix<f64>> for Matrix<f64> {
    type Output = Matrix<f64>;

    fn add(self, rhs: Matrix<f64>) -> Self::Output {
        self.element_wise_op(rhs, ElementWiseOperation::Add)
    }
}

impl std::ops::Sub<Matrix<f64>> for Matrix<f64> {
    type Output = Matrix<f64>;

    fn sub(self, rhs: Matrix<f64>) -> Self::Output {
        self.element_wise_op(rhs, ElementWiseOperation::Subtract)
    }
}

impl std::ops::Rem<Matrix<f64>> for Matrix<f64> {
    type Output = Matrix<f64>;

    fn rem(self, rhs: Matrix<f64>) -> Self::Output {
        self.element_wise_op(rhs, ElementWiseOperation::Modulo)
    }
}

impl std::ops::Mul<f64> for Matrix<f64> {
    type Output = Matrix<f64>;

    fn mul(self, rhs: f64) -> Self::Output {
        self.element_wise_op_scalar(rhs, ElementWiseOperation::Multiply)
    }
}

impl std::ops::Div<f64> for Matrix<f64> {
    type Output = Matrix<f64>;

    fn div(self, rhs: f64) -> Self::Output {
        self.element_wise_op_scalar(rhs, ElementWiseOperation::Divide)
    }
}

impl std::ops::Add<f64> for Matrix<f64> {
    type Output = Matrix<f64>;

    fn add(self, rhs: f64) -> Self::Output {
        self.element_wise_op_scalar(rhs, ElementWiseOperation::Add)
    }
}

impl std::ops::Sub<f64> for Matrix<f64> {
    type Output = Matrix<f64>;

    fn sub(self, rhs: f64) -> Self::Output {
        self.element_wise_op_scalar(rhs, ElementWiseOperation::Subtract)
    }
}

impl std::ops::Rem<f64> for Matrix<f64> {
    type Output = Matrix<f64>;

    fn rem(self, rhs: f64) -> Self::Output {
        self.element_wise_op_scalar(rhs, ElementWiseOperation::Modulo)
    }
}

pub fn lu_factorize_pivoted(matrix: Matrix<f64>) -> (Matrix<f64>, Matrix<f64>) {
    assert!(
        matrix.n == matrix.m,
        "Matrix must be square for LU Factorization."
    );

    let n = matrix.n;

    let mut a = matrix.clone();
    let mut u = matrix.clone();
    let mut l = Matrix::from_one_element(n, n, 1.0f64);
    let mut p = Matrix::from_one_element(n, n, 1.0f64);

    (0..n).into_iter().for_each(|i| {
        let mut k = i;

        loop {
            if u[(i, i)] != 0.0f64 {
                break;
            }

            u.swap_rows(i, k + 1);
            p.swap_rows(i, k + 1);

            k += 1;
        }

        (i + 1..n).into_iter().for_each(|j| {
            l[(j, i)] = u[(j, i)] / u[(i, i)];

            u.elementwise_on_row_scalar(i, l[(i, j)], ElementWiseOperation::Multiply);
            u.elementwise_on_row(j, u[i].clone(), ElementWiseOperation::Subtract);
        })
    });

    (l, u)
}
