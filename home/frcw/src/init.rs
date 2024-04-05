//! Utility functions for loading graph and partition data.
use crate::graph::{Edge, Graph};
use crate::partition::Partition;
use serde_json::Result as SerdeResult;
use serde_json::Value;
use std::collections::HashMap;
use std::fs;

/// Loads graph and partition data in the NetworkX `adjacency_data` format
/// used by [GerryChain](https://github.com/mggg/gerrychain). Returns a
/// [serde_json::Result] containing a [graph::Graph] and
/// a [partition::Partition] upon a successful load.
///
/// # Arguments
///
/// * `path` - the path of the graph JSON file.
/// * `pop_col` - The column in the graph JSON corresponding to total node
///    population. This column should be integer-valued.
/// * `assignment_col` - A column in the graph JSON corresponding to a
///    a seed partition. This column should be integer-valued and 1-indexed.
/// * `columns` - The metadata columns to sum over (per district).
pub fn from_networkx(
    path: &str,
    pop_col: &str,
    assignment_col: &str,
    columns: Vec<String>,
) -> SerdeResult<(Graph, Partition)> {
    let (graph, data) = match graph_from_networkx(path, pop_col, columns) {
        Ok(v) => v,
        Err(e) => return Err(e),
    };

    let raw_nodes = data["nodes"].as_array().unwrap();
    let assignments: Vec<u32> = raw_nodes
        .iter()
        .enumerate()
        .map(|(i, node)| match &node[assignment_col] {
            serde_json::Value::Number(num) => num.as_u64().unwrap() as u32,
            serde_json::Value::String(ref s) => s.parse::<u32>().unwrap(),
            _ => panic!(
                "{}{}{}",
                "Unexpected entry type in assignment column. ",
                format!(
                    "Found {:?} at index {:?} in column {:?}.",
                    node[assignment_col], i, assignment_col
                ),
                "Please make sure that all entries can be interpreted as positive numbers."
            ),
        })
        .collect();
    let partition = Partition::from_assignments(&graph, &assignments).unwrap();
    return Ok((graph, partition));
}

/// Loads graph data in the NetworkX `adjacency_data` format
/// used by [GerryChain](https://github.com/mggg/gerrychain). Returns a
/// [serde_json::Result] containing a [graph::Graph] and the raw
/// graph JSON tree upon a successful load.
///
/// # Arguments
///
/// * `path` - the path of the graph JSON file.
/// * `pop_col` - The column in the graph JSON corresponding to total node
///    population. This column should be integer-valued.
/// * `columns` - The metadata columns to sum over (per district).
pub fn graph_from_networkx(
    path: &str,
    pop_col: &str,
    columns: Vec<String>,
) -> SerdeResult<(Graph, Value)> {
    // TODO: should load from a generic buffer.
    let raw = fs::read_to_string(path).expect("Could not load graph");
    let data: Value = serde_json::from_str(&raw)?;

    let raw_nodes = data["nodes"].as_array().unwrap();
    let raw_adj = data["adjacency"].as_array().unwrap();
    let num_nodes = raw_nodes.len();
    let mut pops = Vec::<u32>::with_capacity(num_nodes);
    let mut neighbors = Vec::<Vec<usize>>::with_capacity(num_nodes);
    let mut edges = Vec::<Edge>::new();
    let mut edges_start = vec![0 as usize; num_nodes];
    let mut attr = HashMap::new();
    for col in columns.to_vec().into_iter() {
        attr.insert(col, Vec::<String>::with_capacity(num_nodes));
    }

    for (index, (node, adj)) in raw_nodes.iter().zip(raw_adj.iter()).enumerate() {
        edges_start[index] = edges.len();
        let node_neighbors: Vec<usize> = adj
            .as_array()
            .unwrap()
            .into_iter()
            .map(|n| n.as_object().unwrap()["id"].as_u64().unwrap() as usize)
            .collect();
        for col in columns.iter() {
            if let Some(data) = attr.get_mut(col) {
                match node.get(col) {
                    Some(value) => data.push(value.to_string()),
                    None => {
                        eprintln!(
                            "Failed to unwrap at column '{}', value {:?}",
                            col, node[col]
                        );
                        panic!("Unexpected None while unwrapping.");
                    }
                }
            }
        }
        let new_pop = match &node[pop_col] {
            serde_json::Value::Number(num) => num.as_u64().unwrap() as u32,
            serde_json::Value::String(ref s) => s.parse::<u32>().unwrap(),
            _ => panic!(
                "{}{}{}",
                "Unexpected entry type in population column. ",
                format!(
                    "Found {:?} at index {:?} in column {:?}.",
                    node[pop_col], index, pop_col
                ),
                "Please make sure that all entries can be interpreted as positive numbers."
            ),
        };
        pops.push(new_pop);
        neighbors.push(node_neighbors.clone());

        for neighbor in &node_neighbors {
            if neighbor > &index {
                let edge = Edge(index, *neighbor);
                edges.push(edge.clone());
            }
        }
    }

    let total_pop = pops.iter().sum();
    let graph = Graph {
        pops: pops,
        neighbors: neighbors,
        edges: edges.clone(),
        edges_start: edges_start.clone(),
        total_pop: total_pop,
        attr: attr,
    };
    return Ok((graph, data));
}
