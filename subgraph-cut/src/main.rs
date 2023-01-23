
use std::collections::HashMap;
use petgraph::{prelude::UnGraph};
use quizx::{vec_graph::Graph, hash_graph::GraphLike, circuit::Circuit};



fn main() {
    let graph: petgraph::graph::UnGraph<(), ()> = petgraph::Graph::new_undirected();
    let options = metis::Options::default();
    let sep = metis::Graph::new(&graph).vertex_separator(&options).unwrap();
    println!("{:?} {:?} {:?}", sep.left, sep.cut, sep.right);

    let bigraph = vertex_cover::BiGraph { graph, left: sep.left, right: sep.cut };
    println!("{:?}", bigraph.min_vertex_cover());


    //testing vertex separator of zx-diagram
    let c = Circuit::random_pauli_gadget()
    .qubits(20)
    .depth(20)
    .seed(1234)
    .min_weight(2)
    .max_weight(5)
    .build();

    let mut g: Graph = c.clone().to_graph();
    quizx::simplify::full_simp(&mut g);

    let sep = Quizx_vertex_cut_finder(&g);
    println!("vertex separator {:?}", sep);




}

pub fn Quizx_vertex_cut_finder(g :&Graph) -> Vec<usize> {

    if g.tcount()==0 {
        return vec![];
    }

    //let nb_spider = g.num_vertices();

    let vmapping = g.vertices().scan(0, |i, s| { *i += 1; Some((s, *i - 1)) }).collect::<HashMap<_,_>>();
    let reverse_vmapping = g.vertices().scan(0, |i, s| { *i += 1; Some((*i - 1, s)) }).collect::<HashMap<_,_>>();
    let edges = g.edges().map(|x| (vmapping[&x.0] ,vmapping[&x.1] )).collect::<Vec<_>>();

    //For weights

    //let mut weights = vec![];
    //for v in 0..nb_spider {
    //    if *g.phase(reverse_vmapping[&v]).denom() == 4 { weights.push(1); }
    //    else { weights.push(0);}
    //}


    let graph = UnGraph::<i32, ()>::from_edges(&edges);
    let options = metis::Options::default();
    let bad_index_separator = metis::Graph::new(&graph).vertex_separator(&options).unwrap();

    
    let mut separator = vec![];
    for x in bad_index_separator.cut{

        separator.push(reverse_vmapping[&x.index()]);

    }


    separator
}
