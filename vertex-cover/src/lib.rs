use petgraph as px;
use std::collections::{HashMap, HashSet};

pub struct BiGraph<N, E> {
    pub graph: px::graph::UnGraph<N, E>,
    pub left: Vec<px::graph::NodeIndex>,
    pub right: Vec<px::graph::NodeIndex>,
}

impl<N, E> BiGraph<N, E> {
    pub fn min_vertex_cover(&self) -> Vec<px::graph::NodeIndex> {
        use rs_graph::{
            Buildable, Builder, maxflow::MaxFlow, traits::Undirected, 
            traits::GraphSize, traits::GraphIterator, vecgraph::VecGraphBuilder
        };
        let mut builder: VecGraphBuilder<usize> = rs_graph::VecGraph::new_builder();
        let source = builder.add_node();
        let sink = builder.add_node();
        let mut nodes_from = HashMap::new();
        let mut left_nodes = HashSet::new();
        let left_to = self.left.iter()
            .map(|&n| {
                let node = builder.add_node();
                builder.add_edge(source, node);
                nodes_from.insert(node, n);
                left_nodes.insert(node);
                (n, node)
            })
            .collect::<HashMap<_,_>>();
        let mut right_nodes = HashSet::new();
        let right_to = self.right.iter()
            .map(|&n| {
                let node = builder.add_node();
                builder.add_edge(node, sink);
                nodes_from.insert(node, n);
                right_nodes.insert(node);
                (n, node)
            })
            .collect::<HashMap<_,_>>();

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

        let mut cut_source = flow.mincut()
            .into_iter()
            .collect::<HashSet<_>>();
        let mut cut_sink = graph.nodes_iter()
            .iter(&graph)
            .filter(|n| !cut_source.contains(n))
            .collect::<HashSet<_>>();

        if cut_sink.contains(&source) {
            std::mem::swap(&mut cut_source, &mut cut_sink);
        }

        let mut out = HashSet::new();
        out.extend(left_nodes.intersection(&cut_sink));
        out.extend(right_nodes.intersection(&cut_source));

        out.into_iter()
            .map(|n| nodes_from[n])
            .collect::<Vec<_>>()
    }
}
