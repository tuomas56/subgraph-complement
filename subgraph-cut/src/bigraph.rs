use crate::rank::RankDecomposition;
use metis::VertexSeparator;
use petgraph as px;
use px::stable_graph::{NodeIndex, StableUnGraph};
use quizx::linalg::Mat2;
use std::collections::{HashMap, HashSet};

#[derive(Clone)]
pub struct BiGraph<N, E> {
    pub graph: px::stable_graph::StableUnGraph<N, E>,
    pub left: Vec<px::graph::NodeIndex>,
    pub right: Vec<px::graph::NodeIndex>,
}

impl<N: Clone, E: Clone> BiGraph<N, E> {
    pub fn from_sep(g: &StableUnGraph<N, E>, sep: &VertexSeparator<StableUnGraph<N, E>>) -> Self {
        let mut sep = sep.clone();
        if sep.left.len() < sep.right.len() {
            let temp = sep.right;
            sep.right = sep.left;
            sep.left = temp;
        }

        let mut subgraph = g.clone();
        subgraph.retain_edges(|g, e| {
            let (a, b) = g.edge_endpoints(e).unwrap();
            (sep.left.contains(&a) && sep.cut.contains(&b))
                || (sep.left.contains(&b) && sep.cut.contains(&a))
        });
        subgraph.retain_nodes(|g, n| g.neighbors(n).count() > 0);
        let new_left = subgraph
            .node_indices()
            .filter(|n| sep.left.contains(n))
            .collect::<Vec<_>>();
        let new_right = subgraph
            .node_indices()
            .filter(|n| sep.cut.contains(n))
            .collect::<Vec<_>>();

        BiGraph {
            graph: subgraph,
            left: new_left,
            right: new_right,
        }
    }

    pub fn min_vertex_cover(&self) -> Vec<px::graph::NodeIndex> {
        use rs_graph::{
            maxflow::MaxFlow, traits::GraphIterator, traits::GraphSize, traits::Undirected,
            vecgraph::VecGraphBuilder, Buildable, Builder,
        };
        let mut builder: VecGraphBuilder<usize> = rs_graph::VecGraph::new_builder();
        let source = builder.add_node();
        let sink = builder.add_node();
        let mut nodes_from = HashMap::new();
        let mut left_nodes = HashSet::new();
        let left_to = self
            .left
            .iter()
            .map(|&n| {
                let node = builder.add_node();
                builder.add_edge(source, node);
                nodes_from.insert(node, n);
                left_nodes.insert(node);
                (n, node)
            })
            .collect::<HashMap<_, _>>();
        let mut right_nodes = HashSet::new();
        let right_to = self
            .right
            .iter()
            .map(|&n| {
                let node = builder.add_node();
                builder.add_edge(node, sink);
                nodes_from.insert(node, n);
                right_nodes.insert(node);
                (n, node)
            })
            .collect::<HashMap<_, _>>();

        for edge in self.graph.edge_indices() {
            if let Some((s, e)) = self.graph.edge_endpoints(edge) {
                if let Some(&s) = left_to.get(&s) {
                    if let Some(&e) = right_to.get(&e) {
                        builder.add_edge(s, e);
                    }
                }

                if let Some(&s) = right_to.get(&s) {
                    if let Some(&e) = left_to.get(&e) {
                        builder.add_edge(e, s);
                    }
                }
            }
        }

        let graph = builder.into_graph();
        let mut flow = rs_graph::maxflow::PushRelabel::new(&graph);
        flow.solve(source, sink, |e| {
            let (a, b) = graph.enodes(e);
            if a == source || b == sink || a == sink || b == source {
                ordered_float::OrderedFloat(1.0)
            } else {
                ordered_float::OrderedFloat(f32::INFINITY)
            }
        });

        let mut cut_source = flow.mincut().into_iter().collect::<HashSet<_>>();
        let mut cut_sink = graph
            .nodes_iter()
            .iter(&graph)
            .filter(|n| !cut_source.contains(n))
            .collect::<HashSet<_>>();

        if cut_sink.contains(&source) {
            std::mem::swap(&mut cut_source, &mut cut_sink);
        }

        let mut out = HashSet::new();
        out.extend(left_nodes.intersection(&cut_sink));
        out.extend(right_nodes.intersection(&cut_source));

        out.into_iter().map(|n| nodes_from[n]).collect::<Vec<_>>()
    }

    pub fn biadjacency(&self) -> Mat2 {
        let mut mat = vec![vec![0; self.left.len()]; self.right.len()];

        let lvmapping = self
            .left
            .clone()
            .into_iter()
            .scan(0, |i, s| {
                *i += 1;
                Some((s, *i - 1))
            })
            .collect::<HashMap<_, _>>();
        
        let rvmapping = self
            .right
            .clone()
            .into_iter()
            .scan(0, |i, s| {
                *i += 1;
                Some((s, *i - 1))
            })
            .collect::<HashMap<_, _>>();

        for &v in &self.left {
            for n in self.graph.neighbors(v) {
                mat[rvmapping[&n]][lvmapping[&v]] = 1;
            }
        }
        Mat2::new(mat)
    }

    pub fn random(n: usize, p: f64, frac: f32) -> BiGraph<(), ()> {
        let mut rng = rand::thread_rng();

        let g: px::graph::UnGraph<(), ()> = petgraph_gen::random_gnp_graph(&mut rng, n, p);
        let mut subgraph = petgraph::stable_graph::StableUnGraph::from(g);

        let vs: Vec<_> = subgraph.node_indices().collect();
        let (left, right) = vs.split_at((n as f32 * frac) as usize);
        let left = left.to_vec();

        let right = right.to_vec();

        subgraph.retain_edges(|g, e| {
            let (a, b) = g.edge_endpoints(e).unwrap();
            (left.contains(&a) && right.contains(&b)) || (left.contains(&b) && right.contains(&a))
        });

        BiGraph {
            graph: subgraph,
            left,
            right,
        }
    }

    pub fn complement_cover(&self) -> Vec<Vec<NodeIndex>> {
        fn getcol(m: &Mat2, c: usize) -> Vec<u8> {
            let mut col = vec![0u8; m.num_rows()];
            for i in 0..m.num_rows() {
                col[i] = m[i][c]
            }
            col
        }

        fn getrow(m: &Mat2, r: usize) -> Vec<u8> {
            let mut row = vec![0u8; m.num_cols()];
            for i in 0..m.num_cols() {
                row[i] = m[r][i]
            }
            row
        }

        let mat = self.biadjacency();
        let (c, r) = mat.rank_decomposition();
        let rank = c.num_cols();
        let sc: Vec<_> = (0..rank).map(|x| (getcol(&c, x), getrow(&r, x))).collect();

        let mut sgcs: Vec<Vec<NodeIndex>> = vec![];
        for sgci in sc {
            let mut sgc: Vec<NodeIndex> = vec![];

            for (i, &v) in sgci.0.iter().enumerate() {
                if v == 1 {
                    sgc.push(self.right[i as usize]);
                }
            }

            for (i, &v) in sgci.1.iter().enumerate() {
                if v == 1 {
                    sgc.push(self.left[i as usize]);
                }
            }

            sgcs.push(sgc);
        }

        sgcs
    }

    pub fn crossing_edges(&self) -> usize {
        self.graph
            .edge_indices()
            .filter(|&e| {
                let (a, b) = self.graph.edge_endpoints(e).unwrap();
                (self.left.contains(&a) && self.right.contains(&b))
                    || (self.left.contains(&b) && self.right.contains(&a))
            })
            .count()
    }
}
