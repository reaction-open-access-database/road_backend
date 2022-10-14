use pyo3::{create_exception, PyAny};
use serde::{Serialize, Deserialize, Serializer, de};
use std::collections::HashMap;
use periodic_table_on_an_enum::Element;

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
    Not {
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
pub enum MolecularWeight {
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
    MolecularWeight {
        molecular_weight: MolecularWeight,
    },
    MolecularFormula {
        atoms: HashMap<SerializableElement, u32>, // Element, count
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

#[derive(Hash, Eq, PartialEq)]
pub struct SerializableElement(pub Element);

impl<'de> Deserialize<'de> for SerializableElement {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let s = match String::deserialize(deserializer) {
            Ok(s) => s,
            Err(e) => return Err(de::Error::custom(e)),
        };
        let element = match Element::from_symbol(&s) {
            Some(element) => element,
            None => return Err(de::Error::custom("Invalid element symbol")),
        };
        Ok(SerializableElement(element))
    }
}

impl Serialize for SerializableElement {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_str(&self.0.get_symbol())
    }
}