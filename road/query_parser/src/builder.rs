use pyo3::{PyResult, ToPyObject};
use pyo3::types::PyDict;
use crate::types::{Modifier, Molecule, Quantity, Query, QueryParserError, StructureOp, Q, AMW};

pub fn build_query(query: Query, q: &Q) -> PyResult<&Q> {
    match query {
        Query::Modifier { query } => build_modifier(query, q),
        Query::Quantity { query} => build_quantity(query, q),
    }
}

fn build_modifier(modifier: Modifier, q: &Q) -> PyResult<&Q> {
    match modifier {
        Modifier::And { queries } => {
            reduce_queries(queries, and, q)
        }
        Modifier::Or{ queries } => {
            reduce_queries(queries, or, q)
        }
        Modifier::Not { query } => {
            not(build_query(*query, q)?)
        }
    }
}

fn reduce_queries<'a, F>(queries: Vec<Query>, mut f: F, q: &'a Q) -> PyResult<&'a Q>
    where F: FnMut(&'a Q, &'a Q) -> PyResult<&'a Q> {
    let reduced_query = queries
        .into_iter()
        .map(|query| build_query(query, q)) // Build each child query
        .reduce(|a, b| f(a?, b?));

    match reduced_query {
        Some(result) => result,
        None => Err(QueryParserError::new_err("No arguments provided")),
    }
}

fn build_quantity(quantity: Quantity, q: &Q) -> PyResult<&Q> {
    match quantity {
        Quantity::Structure {op, value} => {
            let value = match value {
                Molecule::Smiles { value: smiles } => smiles,
            };
            create_q(match op {
                StructureOp::Equal => "molecule__equal",
                StructureOp::HasSubstruct => "molecule__has_substruct",
                StructureOp::IsSubstruct => "molecule__is_substruct",
            }, value, q)
        },
        Quantity::AMW { amw } => {
            match amw {
                AMW::Equal { value, tolerance } => {
                    let lower = value - tolerance;
                    let upper = value + tolerance;
                    and(
                        create_q("molecule__amw__gte", lower, q)?,
                        create_q("molecule__amw__lte", upper, q)?
                    )
                },
                AMW::GreaterThan { value } => {
                    create_q("molecule__amw__gt", value, q)
                },
                AMW::GreaterThanOrEqual { value } => {
                    create_q("molecule__amw__gte", value, q)
                },
                AMW::LessThan { value } => {
                    create_q("molecule__amw__lt", value, q)
                },
                AMW::LessThanOrEqual { value } => {
                    create_q("molecule__amw__lte", value, q)
                },
            }
        }
    }
}

fn create_q<'a, T>(key: &'a str, value: T, q: &'a Q) -> PyResult<&'a Q>
    where T: ToPyObject {
    let pydict = PyDict::new(q.py());
    pydict.set_item(key, value)?;
    q.call((), Some(pydict))
}

// Modifier functions
fn and<'a>(x: &'a Q, y: &'a Q) -> PyResult<&'a Q> {
    x.call_method1("__and__", (y,))
}

fn or<'a>(x: &'a Q, y: &'a Q) -> PyResult<&'a Q> {
    x.call_method1("__or__", (y,))
}

fn not(x: &Q) -> PyResult<&Q> {
    x.call_method0("__invert__")
}
