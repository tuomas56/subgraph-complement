#![allow(dead_code)]

use crate::anneal::{ComplementFinder, GeometricSeries};
use petgraph as px;
use px::stable_graph::StableUnGraph;
use quizx::{
    circuit::Circuit,
    hash_graph::GraphLike,
    vec_graph::{BasisElem, Graph},
};
use std::{collections::HashMap, vec};

mod anneal;
mod bigraph;
mod rank;

trait GraphUtils {
    type Node;
    type Edge;

    fn to_petgraph(&self) -> StableUnGraph<Self::Node, Self::Edge>;
    fn to_qasm(&self) -> String;
}

impl<G: quizx::graph::GraphLike> GraphUtils for G {
    type Node = ();
    type Edge = ();

    fn to_petgraph(&self) -> StableUnGraph<Self::Node, Self::Edge> {
        let vmapping = self
            .vertices()
            .scan(0, |i, s| {
                *i += 1;
                Some((s, *i - 1))
            })
            .collect::<HashMap<_, _>>();
        let edges = self
            .edges()
            .map(|x| (vmapping[&x.0], vmapping[&x.1]))
            .collect::<Vec<_>>();
        let graph = px::stable_graph::StableUnGraph::<(), ()>::from_edges(&edges);
        graph
    }

    fn to_qasm(&self) -> String {
        let mut g = self.clone();
        g.x_to_z();
        quizx::simplify::spider_simp(&mut g);

        let indices = g.vertices()
            .enumerate()
            .map(|(i, j)| (j, i))
            .collect::<HashMap<_, _>>();
        
        let mut c = quizx::circuit::Circuit::new(indices.len());
        
        for v in g.vertices() {
            c.add_gate("h", vec![indices[&v]]);
            c.add_gate_with_phase("rz", vec![indices[&v]], g.phase(v));
        }

        for (a, b, _) in g.edges() {
            c.add_gate("cz", vec![indices[&a], indices[&b]])
        }

        for v in g.vertices() {
            c.add_gate("h", vec![indices[&v]]);
        }

        c.to_qasm()
    }
}

fn main() {
    let c = Circuit::random()
        .qubits(60)
        .depth(2000)
        .seed(35125)
        .clifford_t(0.1)
        .build();

    let mut g: Graph = c.clone().to_graph();
    g.plug_outputs(&vec![BasisElem::Z0; c.num_qubits()]);
    g.plug_inputs(&vec![BasisElem::Z0; c.num_qubits()]);

    quizx::simplify::full_simp(&mut g);

    println!("{:?}", g.num_vertices());

    let mut gp = g.to_petgraph();
    let mut rng = rand::thread_rng();

    for i in 1..=10 {
        let mut finder = ComplementFinder::new(
            &gp,
            &mut rng,
            GeometricSeries::new(0.02, 0.001, 10000),
            i,
            10,
        );
        finder.run(false);
        gp = finder.graph;
        println!("--------------------------")
    }
}
