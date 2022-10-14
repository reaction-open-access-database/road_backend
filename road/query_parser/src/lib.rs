mod types;
mod builder;

use std::collections::HashMap;
use pyo3::prelude::*;
use crate::types::{QueryParserError, SerializableElement};
use periodic_table_on_an_enum::Element;

#[pyfunction]
fn build_molecule_query<'a>(input: &'a str, q: &'a PyAny) -> PyResult<&'a PyAny> {
    let parsed_query = match serde_json::from_str(input) {
        Ok(query) => query,
        Err(e) => return Err(types::QueryParserError::new_err(e.to_string())),
    };
    builder::build_query(parsed_query, q)
}

#[pyfunction]
fn generate_example_json() -> PyResult<String> {
    let mut atom_hashmap = HashMap::new();
    atom_hashmap.insert(SerializableElement(Element::Carbon), 1);
    atom_hashmap.insert(SerializableElement(Element::Hydrogen), 4);

    let query = types::Query::Modifier {
        query: types::Modifier::And {
            queries: vec![
                types::Query::Quantity {
                    query: types::Quantity::Structure {
                        op: types::StructureOp::Equal,
                        value: types::Molecule::Smiles { value: "C".to_string() },
                    },
                },
                types::Query::Quantity {
                    query: types::Quantity::Structure {
                        op: types::StructureOp::HasSubstruct,
                        value: types::Molecule::Smiles { value: "C".to_string() },
                    },
                },
                types::Query::Quantity {
                    query: types::Quantity::MolecularWeight {
                        molecular_weight: types::MolecularWeight::Equal {
                            value: 0.0,
                            tolerance: 0.0
                        },
                    },
                },
                types::Query::Quantity {
                    query: types::Quantity::MolecularFormula {
                        atoms: atom_hashmap,
                    }
                }
            ],
        },
    };
    Ok(serde_json::to_string(&query).unwrap())
}

#[pymodule]
fn query_parser(py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(build_molecule_query, m)?)?;
    m.add_function(wrap_pyfunction!(generate_example_json, m)?)?;
    m.add("QueryParserError", py.get_type::<QueryParserError>())?;
    Ok(())
}