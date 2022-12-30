use crate::types::{
    Modifier, MolecularWeight, Molecule, Quantity, Query, QueryParserError, StructureOp, Q,
};
use periodic_table_on_an_enum::Element;
use pyo3::types::PyDict;
use pyo3::{PyResult, ToPyObject};
use std::cmp::Ordering;
use std::collections::HashMap;

pub fn build_query(query: Query, q: &Q) -> PyResult<&Q> {
    match query {
        Query::Modifier { query } => build_modifier(query, q),
        Query::Quantity { query } => build_quantity(query, q),
    }
}

fn build_modifier(modifier: Modifier, q: &Q) -> PyResult<&Q> {
    match modifier {
        Modifier::And { queries } => reduce_queries(queries, and, q),
        Modifier::Or { queries } => reduce_queries(queries, or, q),
        Modifier::Not { query } => not(build_query(*query, q)?),
    }
}

fn reduce_queries<'a, F>(queries: Vec<Query>, mut f: F, q: &'a Q) -> PyResult<&'a Q>
where
    F: FnMut(&'a Q, &'a Q) -> PyResult<&'a Q>,
{
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
        Quantity::Structure { op, value } => {
            let value = match value {
                Molecule::Smiles { value: smiles } => smiles,
            };
            create_q(
                match op {
                    StructureOp::Equal => "molecule",
                    StructureOp::HasSubstruct => "molecule__has_substruct",
                    StructureOp::IsSubstruct => "molecule__is_substruct",
                },
                value,
                q,
            )
        }
        Quantity::MolecularWeight {
            molecular_weight: amw,
        } => match amw {
            MolecularWeight::Equal { value, tolerance } => {
                let lower = value - tolerance;
                let upper = value + tolerance;
                and(
                    create_q("molecule__amw__gte", lower, q)?,
                    create_q("molecule__amw__lte", upper, q)?,
                )
            }
            MolecularWeight::GreaterThan { value } => create_q("molecule__amw__gt", value, q),
            MolecularWeight::GreaterThanOrEqual { value } => {
                create_q("molecule__amw__gte", value, q)
            }
            MolecularWeight::LessThan { value } => create_q("molecule__amw__lt", value, q),
            MolecularWeight::LessThanOrEqual { value } => create_q("molecule__amw__lte", value, q),
        },
        Quantity::MolecularFormula { atoms } => {
            let atoms = atoms.into_iter().map(|(e, c)| (e.0, c)).collect();
            let molecular_formula = create_molecular_formula(atoms);
            create_q("molecular_formula", molecular_formula, q)
        }
    }
}

fn create_molecular_formula(atoms: HashMap<Element, u32>) -> String {
    let mut atom_list: Vec<(Element, u32)> = atoms.into_iter().collect();
    atom_list.sort_by(|a, b| {
        if a.0 == Element::Carbon {
            Ordering::Less
        } else if b.0 == Element::Carbon {
            Ordering::Greater
        } else {
            a.0.cmp(&b.0)
        }
    });
    atom_list
        .iter()
        .map(|(element, count)| {
            if *count == 1 {
                element.get_symbol().to_string()
            } else {
                format!("{}{}", element.get_symbol(), count)
            }
        })
        .collect::<Vec<String>>()
        .join("")
}

fn create_q<'a, T>(key: &'a str, value: T, q: &'a Q) -> PyResult<&'a Q>
where
    T: ToPyObject,
{
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

#[cfg(test)]
mod tests {
    #[test]
    fn test_create_molecular_formula() {
        use super::create_molecular_formula;
        use periodic_table_on_an_enum::Element;
        use std::collections::HashMap;

        // Methane
        let mut atoms = HashMap::new();
        atoms.insert(Element::Carbon, 1);
        atoms.insert(Element::Hydrogen, 4);
        assert_eq!(create_molecular_formula(atoms), "CH4");

        // Ethane
        let mut atoms = HashMap::new();
        atoms.insert(Element::Carbon, 2);
        atoms.insert(Element::Hydrogen, 6);
        assert_eq!(create_molecular_formula(atoms), "C2H6");

        let mut atoms = HashMap::new();
        atoms.insert(Element::Carbon, 1);
        atoms.insert(Element::Hydrogen, 3);
        atoms.insert(Element::Nitrogen, 1);
        atoms.insert(Element::Oxygen, 1);
        assert_eq!(create_molecular_formula(atoms), "CH3NO");

        let mut atoms = HashMap::new();
        atoms.insert(Element::Carbon, 1);
        atoms.insert(Element::Oxygen, 2);
        assert_eq!(create_molecular_formula(atoms), "CO2");
    }
}
