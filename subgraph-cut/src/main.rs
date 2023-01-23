
use std::collections::HashMap;
use petgraph as px;
use quizx::{vec_graph::Graph, hash_graph::GraphLike, circuit::Circuit};

mod anneal;

fn main() {
    //testing vertex separator of zx-diagram
    let c = Circuit::random_pauli_gadget()
    .qubits(10)
    .depth(10)
    .seed(1234)
    .min_weight(2)
    .max_weight(5)
    .build();

    let mut g: Graph = c.clone().to_graph();
    quizx::simplify::full_simp(&mut g);

    println!("{:?}", g.num_vertices());

    let gp = quizx_to_petgraph(&g);
    let sep = metis::Graph::new(&gp)
        .vertex_separator(&metis::Options::default()
            .max_imbalance(10))
        .unwrap();

    println!("left = {:?}\n cut={:?}\n right={:?}", sep.left, sep.cut, sep.right);
    println!("all = {:?}", gp.node_indices().collect::<Vec<_>>());
    
    let mut subgraph = px::stable_graph::StableUnGraph::from(gp);
    subgraph.retain_edges(|g, e| {
        let (a, b) = g.edge_endpoints(e).unwrap();
        (sep.left.contains(&a) && sep.cut.contains(&b)) || (sep.left.contains(&b) && sep.cut.contains(&a))
    });
    subgraph.retain_nodes(|g, n| {
        g.neighbors(n).count() > 0
    });
    let new_left = subgraph.node_indices()
        .filter(|n| sep.left.contains(n))
        .collect::<Vec<_>>();
    let new_right = subgraph.node_indices()
        .filter(|n| sep.cut.contains(n))
        .collect::<Vec<_>>();
    
    let gp = &subgraph;
    println!("left = {} cut = {} nodes = {} edges = {} density = {}", new_left.len(), new_right.len(), gp.node_count(), gp.edge_count(), 2.0 * gp.edge_count() as f32 / (gp.node_count() as f32 * (gp.node_count() as f32 - 1.0)));
    
    let bg = vertex_cover::BiGraph { graph: subgraph, left: new_left, right: new_right };
    
    println!("vertex cover = {}", bg.min_vertex_cover().len());
    
    let mut rng = rand::thread_rng();
    let mut finder = anneal::ComplementFinder::new(&bg, &mut rng, anneal::GeometricSeries::new(1000.0, 1.0, 100000));
    finder.run(false);
    
    println!("vertex cover = {}", finder.graph.min_vertex_cover().len());
    let gp = finder.graph.graph;
    println!("{}", 2.0 * gp.edge_count() as f32 / (gp.node_count() as f32 * (gp.node_count() as f32 - 1.0)));
}

fn quizx_to_petgraph(g: &Graph) -> px::stable_graph::StableUnGraph<(), ()> {
    let vmapping = g.vertices().scan(0, |i, s| { *i += 1; Some((s, *i - 1)) }).collect::<HashMap<_,_>>();
    let edges = g.edges().map(|x| (vmapping[&x.0] ,vmapping[&x.1] )).collect::<Vec<_>>();
    let graph = px::stable_graph::StableUnGraph::<(), ()>::from_edges(&edges);
    graph
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


    let graph = px::graph::UnGraph::<i32, ()>::from_edges(&edges);
    let options = metis::Options::default();
    let bad_index_separator = metis::Graph::new(&graph).vertex_separator(&options).unwrap();

    
    let mut separator = vec![];
    for x in bad_index_separator.cut{

        separator.push(reverse_vmapping[&x.index()]);

    }


    separator
}
