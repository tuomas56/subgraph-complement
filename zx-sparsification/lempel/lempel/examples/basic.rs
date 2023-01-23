use rand::Rng;
use zx::graph::{BasisElem, GraphLike};
use zx::hash_graph::Graph;
use zx::circuit::Circuit;
use lempel::MatF2;
use std::collections::HashMap;
use kahypar::{Context, Hypergraph};

struct Bijection<U, V> {
    forward: HashMap<U, V>,
    reverse: HashMap<V, U>
}

impl<U, V> std::iter::FromIterator<(U, V)> for Bijection<U, V> where 
    U: Clone + std::hash::Hash + PartialEq + Eq, 
    V: Clone + std::hash::Hash + PartialEq + Eq {
    fn from_iter<T: IntoIterator<Item = (U, V)>>(iter: T) -> Self {
        let mut forward = HashMap::new();
        let mut reverse = HashMap::new();
        for (u, v) in iter {
            forward.insert(u.clone(), v.clone());
            reverse.insert(v, u);
        }
        Bijection { forward, reverse }
    }
}

fn adjacency_matrix(g: &impl GraphLike) -> (MatF2, Bijection<usize, usize>) {
    let n = g.num_vertices();
    let map = g.vertices()
        .enumerate()
        .collect::<Bijection<_, _>>();
    let mat = MatF2::build(n, n, |i, j| {
        g.connected(map.forward[&i], map.forward[&j])
    });
    (mat, map)
}

fn build_hypergraph(hedges: &MatF2) -> Hypergraph {
    let mut incidence = HashMap::new();
    for vertex in 0..hedges.rows {
        for edge in 0..hedges.cols {
            incidence.insert((vertex, edge), hedges.get(vertex, edge));
        }
    }

    struct Incidence(HashMap<(usize, usize), bool>);
    impl std::ops::Index<(usize, usize)> for Incidence {
        type Output = bool;
        fn index(&self, index: (usize, usize)) -> &Self::Output {
            &self.0[&index]
        }
    }

    Hypergraph::from_incidence(
        2, hedges.rows, hedges.cols, 
        Incidence(incidence)
    ).build()
}

fn hamming_weight(m: &MatF2) -> usize {
    let mut weight = 0;
    for i in 0..m.rows {
        for j in 0..m.cols {
            weight += m.get(i, j) as usize;
        }
    }
    weight
}

pub fn main() {
    let mut g = Circuit::random()
        .depth(2000)
        .qubits(50)
        .clifford_t(0.1)
        .build()
        .to_graph::<Graph>();

    g.plug_outputs(&[BasisElem::X0; 50]);
    g.plug_inputs(&[BasisElem::X0; 50]);

    zx::simplify::full_simp(&mut g);

    let (m, map) = adjacency_matrix(&g);

    println!("m =\n{}{}x{}", m, m.rows, m.cols);
    let b = m.rank_decomposition(false);
    println!("b =\n{}{}x{}", b, b.rows, b.cols);

    assert_eq!(b.clone() * b.transpose(), m);

    println!("plain: {} {}", hamming_weight(&m), hamming_weight(&b));

    for edge in 0..b.cols {
        let mut bp = MatF2::build(b.rows, b.cols, |i, _| b.get(i, edge));
        let mut bpp = b.clone();
        bpp.add_from(&bp);
        println!("{} {}", hamming_weight(&b), hamming_weight(&bpp));
    }

}