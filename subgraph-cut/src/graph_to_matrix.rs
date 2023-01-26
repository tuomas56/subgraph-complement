use std::collections::HashMap;

use petgraph::prelude::UnGraph;
use vertex_cover::BiGraph;
use quizx::linalg::Mat2;


pub fn bipartite_to_biadjaency(bg :&BiGraph<(),()>) -> Mat2{

    let mut mat = vec![vec![0; bg.left.len()]; bg.right.len()];

    let lvmapping = bg.left.clone().into_iter().scan(0, |i, s| { *i += 1; Some((s, *i - 1)) }).collect::<HashMap<_,_>>();
    let rvmapping = bg.right.clone().into_iter().scan(0, |i, s| { *i += 1; Some((s, *i - 1)) }).collect::<HashMap<_,_>>();

    for &v in &bg.left {
        for n in bg.graph.neighbors(v){
            mat[rvmapping[&n]][lvmapping[&v]]=1;
        }
    }
    Mat2::new(mat)
}

#[test]
fn test_matrix() {

    let bg = random_bigraph(9,0.5,0.5);

    println!("{:?}", bipartite_to_biadjaency(&bg));
    println!("{:?}",petgraph::dot::Dot::new(&bg.graph));
}

pub fn random_bigraph(n:usize,p:f64,frac:f32) ->BiGraph<(),()>{


    let mut rng = rand::thread_rng();

    let g: UnGraph<(), ()> = petgraph_gen::random_gnp_graph(&mut rng, n, p);
    let mut subgraph = petgraph::stable_graph::StableUnGraph::from(g);


    let vs: Vec<_> = subgraph.node_indices().collect();
    let (left,right) = vs.split_at((n as f32*frac)as usize);
    let left=left.to_vec();

    let right=right.to_vec();


    subgraph.retain_edges(|g, e| {
        let (a, b) = g.edge_endpoints(e).unwrap();
        (left.contains(&a) && right.contains(&b)) || (left.contains(&b) && right.contains(&a))
    });
    
    vertex_cover::BiGraph { graph: subgraph, left, right }
}