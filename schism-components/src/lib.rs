#![allow(unused)]

use std::marker::PhantomData;

type NonLinearFunction<'a> = &'a dyn FnOnce(&[f64]) -> f64;

pub trait BasicComponent { }


#[derive(Clone)]
pub struct AcSource {
    amplitude: f64,
    freq: f64,
    init: f64,
}

impl BasicComponent for AcSource { }

#[derive(Clone)]
pub enum DcSource {
    Current(f64),
    Voltage(f64),
    CurrentByVoltage(f64, f64),
    VoltageByCurrent(f64, f64)
}

impl BasicComponent for DcSource { }

#[derive(Clone)]
pub enum Resistor<'a> {
    Linear(f64),
    NonLinear(NonLinearFunction<'a>),
}

impl<'a> BasicComponent for Resistor<'a> { }

#[derive(Clone)]
pub enum Capacitor<'a> {
    Linear(f64),
    NonLinear(NonLinearFunction<'a>),
}

impl<'a> BasicComponent for Capacitor<'a> { }

#[derive(Clone)]
pub enum Inductor<'a> {
    Linear(f64),
    NonLinear(NonLinearFunction<'a>)
}

impl<'a> BasicComponent for Inductor<'a> { }

pub struct NonLinearElement<'a>(NonLinearFunction<'a>);


pub struct Circuit<'a> {
    schematic: PhantomData<&'a u8>,
    mapper: PhantomData<&'a u8>,
    nodal_analytics: PhantomData<&'a u8>
}