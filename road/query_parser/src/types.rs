use pyo3::{create_exception, PyAny};
use serde::{Serialize, Deserialize};

create_exception!(query_parser, QueryParserError, pyo3::exceptions::PyException);

pub type Q = PyAny;

#[derive(Serialize, Deserialize)]
#[serde(rename_all = "lowercase", tag = "type")]
pub enum Modifier {
    And {
        queries: Vec<Query>,
    },
    Or {
        queries: Vec<Query>,
    },
    Not{
        query: Box<Query>,
    },
}

#[derive(Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum StructureOp {
    Equal,
    HasSubstruct,
    IsSubstruct,
}

#[derive(Serialize, Deserialize)]
#[serde(rename_all = "lowercase", tag = "type")]
pub enum AMW {
    Equal {
        value: f64,
        tolerance: f64,
    },
    GreaterThan {
        value: f64,
    },
    GreaterThanOrEqual {
        value: f64,
    },
    LessThan {
        value: f64,
    },
    LessThanOrEqual {
        value: f64,
    },
}

#[derive(Serialize, Deserialize)]
#[serde(rename_all = "lowercase", tag = "type")]
pub enum Quantity {
    // Operator enum, value
    Structure {
        op: StructureOp,
        value: Molecule
    },
    AMW {
        amw: AMW,
    },
}

#[derive(Serialize, Deserialize)]
#[serde(rename_all = "lowercase", tag = "type")]
pub enum Query {
    Modifier {
        query: Modifier,
    },
    Quantity {
        query: Quantity,
    },
}

#[derive(Serialize, Deserialize)]
#[serde(rename_all = "lowercase", tag = "type")]
pub enum Molecule {
    Smiles {
        value: String,
    },
}