use snafu::prelude::*;

/// Data structures for partitionings (districting plans).
use crate::buffers::SubgraphBuffer;
use crate::graph::{Edge, Graph};
use crate::recom::RecomProposal;

#[derive(Debug, PartialEq, Snafu)]
pub enum PartitionError {
    #[snafu(display("Empty assignment vector"))]
    ErrEmptyAssignmentVector,
    #[snafu(display("Assignments must be 1-indexed"))]
    ErrAssignmentVectorNotOneIndexed,
    #[snafu(display("District {district_number} has no nodes"))]
    ErrDistrictHasNoNodes { district_number: usize },
    #[snafu(display(
        "Mismatch: graph has {graph_nodes} nodes, assignment vector has {vector_nodes} nodes"
    ))]
    ErrGraphMismatchVector {
        graph_nodes: usize,
        vector_nodes: usize,
    },
    #[snafu(display("Could not parse assignments: {parse_error}"))]
    ErrParseAssignments { parse_error: String },
}

/// A partitioning (districting plan) on top of a [Graph].
/// The graph is referenced implicitly (we don't store a reference to it).
#[derive(Clone, Debug)]
pub struct Partition {
    /// The number of districts (parts) in the partitioning.
    pub num_dists: u32,
    /// An assignment vector mapping nodes in the graph to
    /// district labels.
    pub assignments: Vec<u32>,
    /// The nodes in each district (a list-of-lists representation).
    /// This should be consistent with `assignments`.
    pub dist_nodes: Vec<Vec<usize>>,
    /// The population in each district.
    pub dist_pops: Vec<u32>,
    /// The cut edges (that is, edges that connect nodes in different
    /// districts) in the partitioning.
    /// This should be consistent with `dist_nodes`.
    /// Computed lazily.
    cut_edges: Option<Vec<usize>>,
    /// A flattened district adjacency matrix. A district pair's entry
    /// is the number of cut edges between the pair; nonadjacency is
    /// represented by a cut edge count of 0.
    /// Computed lazily.
    dist_adj: Option<Vec<u32>>,
}

impl Partition {
    /// Updates a [Partition] with an underlying `graph` to reflect a `proposal`.
    pub fn update(&mut self, proposal: &RecomProposal) {
        // Move nodes.
        self.dist_nodes[proposal.a_label] = proposal.a_nodes.clone();
        self.dist_nodes[proposal.b_label] = proposal.b_nodes.clone();
        self.dist_pops[proposal.a_label] = proposal.a_pop;
        self.dist_pops[proposal.b_label] = proposal.b_pop;
        for &node in proposal.a_nodes.iter() {
            self.assignments[node] = proposal.a_label as u32;
        }
        for &node in proposal.b_nodes.iter() {
            self.assignments[node] = proposal.b_label as u32;
        }
        // Reset lazily computed derived properties.
        self.cut_edges = None;
        self.dist_adj = None;
    }

    /// Computes the partition's cut edges.
    pub fn cut_edges(&mut self, graph: &Graph) -> &Vec<usize> {
        if self.cut_edges.is_none() {
            let mut cut_edges = Vec::<usize>::new();
            for (index, edge) in graph.edges.iter().enumerate() {
                let dist_a = self.assignments[edge.0 as usize];
                let dist_b = self.assignments[edge.1 as usize];
                if dist_a != dist_b {
                    cut_edges.push(index);
                }
            }
            self.cut_edges = Some(cut_edges);
        }
        self.cut_edges.as_ref().unwrap()
    }

    /// Computes the partition's district adjacency matrix.
    pub fn dist_adj(&mut self, graph: &Graph) -> &Vec<u32> {
        if self.dist_adj.is_none() {
            let mut dist_adj = vec![0 as u32; (self.num_dists * self.num_dists) as usize];
            for edge in graph.edges.iter() {
                let dist_a = self.assignments[edge.0 as usize];
                let dist_b = self.assignments[edge.1 as usize];
                if dist_a != dist_b {
                    dist_adj[((dist_a * self.num_dists) + dist_b) as usize] += 1;
                    dist_adj[((dist_b * self.num_dists) + dist_a) as usize] += 1;
                }
            }
            self.dist_adj = Some(dist_adj);
        }
        self.dist_adj.as_ref().unwrap()
    }

    /// Copies the subgraph induced by the union of districts `a` and `b`
    /// into a buffer. (Node attributes are omitted.)
    ///
    /// The resulting subgraph has relabeled node IDs: nodes
    /// [0..# of nodes in district `a`] are from district `a`, and the
    /// remaining nodes are from district `b`. The `node_to_idx` member
    /// of the subgraph buffer contains a mapping between the node IDs
    /// of the parent graph and these new node IDs.
    ///
    /// # Arguments
    ///
    /// * `graph` - The underlying graph of the [Partition].
    /// * `buf` - The buffer to copy the nodes into.
    /// * `a` - The label of the `a`-district.
    /// * `b` - The label of the `b`-district.
    pub fn subgraph(&self, graph: &Graph, buf: &mut SubgraphBuffer, a: usize, b: usize) {
        buf.clear();
        buf.raw_nodes.clone_from(&self.dist_nodes[a]);
        buf.raw_nodes.extend_from_slice(&self.dist_nodes[b]);
        for (idx, &node) in buf.raw_nodes.iter().enumerate() {
            buf.node_to_idx[node] = idx as i64;
        }
        let mut edge_pos = 0;
        for (idx, &node) in buf.raw_nodes.iter().enumerate() {
            buf.graph.edges_start[idx] = edge_pos;
            for &neighbor in graph.neighbors[node].iter() {
                if buf.node_to_idx[neighbor] >= 0 {
                    let neighbor_idx = buf.node_to_idx[neighbor] as usize;
                    buf.graph.neighbors[idx].push(neighbor_idx as usize);
                    if neighbor_idx > idx as usize {
                        buf.graph.edges.push(Edge(idx, neighbor_idx as usize));
                        edge_pos += 1;
                    }
                }
            }
            buf.graph.pops.push(graph.pops[node]);
        }
        buf.graph.total_pop = self.dist_pops[a] + self.dist_pops[b];
    }

    /// Copies the subgraph induced by the union of districts `a` and `b`
    /// into a buffer. Similar to `subgraph`, but *selected* node attributes
    /// are also copied.
    ///
    /// # Arguments
    ///
    /// * `graph` - The underlying graph of the [Partition].
    /// * `buf` - The buffer to copy the nodes into.
    /// * `attrs` - The node attributes to copy.
    /// * `a` - The label of the `a`-district.
    /// * `b` - The label of the `b`-district.
    pub fn subgraph_with_attr_subset<'a>(
        &self,
        graph: &Graph,
        buf: &mut SubgraphBuffer,
        attrs: impl Iterator<Item = &'a String>,
        a: usize,
        b: usize,
    ) {
        self.subgraph(graph, buf, a, b);
        for key in attrs {
            let vals = graph.attr.get(key).unwrap();
            if buf.graph.attr.contains_key(key) {
                buf.graph.attr.get_mut(key).unwrap().clear();
            } else {
                buf.graph.attr.insert(
                    key.to_string(),
                    Vec::<String>::with_capacity(buf.raw_nodes.len()),
                );
            }
            let buf_vals = buf.graph.attr.get_mut(key).unwrap();
            for &node in buf.raw_nodes.iter() {
                buf_vals.push(vals[node].clone());
            }
        }
    }

    /// Copies the subgraph induced by the union of districts `a` and `b`
    /// into a buffer. Similar to `subgraph`, but *all* node attributes are
    /// also copied.
    ///
    /// # Arguments
    ///
    /// * `graph` - The underlying graph of the [Partition].
    /// * `buf` - The buffer to copy the nodes into.
    /// * `a` - The label of the `a`-district.
    /// * `b` - The label of the `b`-district.
    pub fn subgraph_with_attr(&self, graph: &Graph, buf: &mut SubgraphBuffer, a: usize, b: usize) {
        self.subgraph_with_attr_subset(graph, buf, graph.attr.keys(), a, b);
    }

    /// Builds a partition from a 1-indexed assignment vector.
    pub fn from_assignments(
        graph: &Graph,
        assignments: &Vec<u32>,
    ) -> Result<Partition, PartitionError> {
        match assignments.iter().min() {
            None => return Err(PartitionError::ErrEmptyAssignmentVector),
            Some(1) => (),
            Some(_) => return Err(PartitionError::ErrAssignmentVectorNotOneIndexed),
        };

        if assignments.len() != graph.neighbors.len() {
            return Err(PartitionError::ErrGraphMismatchVector {
                graph_nodes: graph.neighbors.len(),
                vector_nodes: assignments.len(),
            });
        }

        let num_dists = *assignments.iter().max().unwrap(); // guaranteed nonempty
        let mut dist_nodes = vec![Vec::<usize>::new(); num_dists as usize];
        let mut dist_pops = vec![0 as u32; num_dists as usize];
        let assignments_zeroed = assignments.iter().map(|a| a - 1).collect::<Vec<u32>>();
        for (node, &assignment) in assignments_zeroed.iter().enumerate() {
            assert!(assignment < num_dists);
            dist_nodes[assignment as usize].push(node);
            dist_pops[assignment as usize] += graph.pops[node];
        }
        for (dist, nodes) in dist_nodes.iter().enumerate() {
            if nodes.is_empty() {
                return Err(PartitionError::ErrDistrictHasNoNodes {
                    district_number: dist + 1,
                });
            }
        }
        let partition = Partition {
            num_dists: num_dists,
            assignments: assignments_zeroed,
            cut_edges: None,
            dist_adj: None,
            dist_pops: dist_pops,
            dist_nodes: dist_nodes,
        };
        Ok(partition)
    }

    /// Builds a partition from a space-delimited string representing a
    /// 1-indexed assignment vector.
    pub fn from_assignment_str(
        graph: &Graph,
        assignments: &str,
    ) -> Result<Partition, PartitionError> {
        match assignments
            .replace('\n', "")
            .split(' ')
            .map(|a| a.parse::<u32>())
            .collect()
        {
            Ok(vec) => Partition::from_assignments(graph, &vec),
            Err(err) => Err(PartitionError::ErrParseAssignments {
                parse_error: err.to_string(),
            }),
        }
    }
}

mod tests {
    #[allow(unused_imports)]
    use super::*;

    #[test]
    fn from_assignments_rect_grid_2x2() {
        let grid = Graph::rect_grid(2, 2);
        let assignments = vec![1, 1, 1, 2];
        let mut partition = Partition::from_assignments(&grid, &assignments).unwrap();
        assert_eq!(partition.num_dists, 2);
        assert_eq!(partition.assignments, vec![0, 0, 0, 1]);
        assert_eq!(partition.dist_pops, vec![3, 1]);
        assert_eq!(partition.dist_nodes, vec![vec![0, 1, 2], vec![3]]);
        assert_eq!(*partition.dist_adj(&grid), vec![0, 2, 2, 0]);
        assert_eq!(*partition.cut_edges(&grid), vec![2, 3]);
    }

    #[test]
    fn from_assignments_zero_indexed() {
        let grid = Graph::rect_grid(2, 2);
        let assignments = vec![0, 0, 0, 1];
        assert_eq!(
            Partition::from_assignments(&grid, &assignments).unwrap_err(),
            PartitionError::ErrAssignmentVectorNotOneIndexed
        );
    }

    #[test]
    fn from_assignments_two_indexed() {
        let grid = Graph::rect_grid(2, 2);
        let assignments = vec![2, 2, 2, 3];
        assert_eq!(
            Partition::from_assignments(&grid, &assignments).unwrap_err(),
            PartitionError::ErrAssignmentVectorNotOneIndexed
        );
    }

    #[test]
    fn from_assignments_missing_district() {
        let grid = Graph::rect_grid(2, 2);
        let assignments = vec![1, 1, 1, 3];
        assert_eq!(
            Partition::from_assignments(&grid, &assignments).unwrap_err(),
            PartitionError::ErrDistrictHasNoNodes { district_number: 2 }
        );
    }

    #[test]
    fn from_assignments_length_mismatch() {
        let grid = Graph::rect_grid(2, 2);
        let assignments = vec![1, 1, 3];
        assert_eq!(
            Partition::from_assignments(&grid, &assignments).unwrap_err(),
            PartitionError::ErrGraphMismatchVector {
                graph_nodes: 4,
                vector_nodes: 3
            }
        );
    }

    #[test]
    fn from_assignments_empty() {
        let grid = Graph::rect_grid(2, 2);
        let assignments = vec![];
        assert_eq!(
            Partition::from_assignments(&grid, &assignments).unwrap_err(),
            PartitionError::ErrEmptyAssignmentVector
        );
    }

    #[test]
    fn from_assignment_str_rect_grid_2x2() {
        let grid = Graph::rect_grid(2, 2);
        let assignments = "1 1 1 2";
        let mut partition = Partition::from_assignment_str(&grid, &assignments).unwrap();
        assert_eq!(partition.num_dists, 2);
        assert_eq!(partition.assignments, vec![0, 0, 0, 1]);
        assert_eq!(partition.dist_pops, vec![3, 1]);
        assert_eq!(partition.dist_nodes, vec![vec![0, 1, 2], vec![3]]);
        assert_eq!(*partition.dist_adj(&grid), vec![0, 2, 2, 0]);
        assert_eq!(*partition.cut_edges(&grid), vec![2, 3]);
    }

    #[test]
    fn from_assignment_str_bad_char() {
        let grid = Graph::rect_grid(2, 2);
        let assignments = "1 1 1 a";
        assert_eq!(
            Partition::from_assignment_str(&grid, &assignments).unwrap_err(),
            PartitionError::ErrParseAssignments {
                parse_error: "invalid digit found in string".to_string()
            }
        );
    }
}
