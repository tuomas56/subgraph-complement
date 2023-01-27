#![allow(dead_code)]

use std::{collections::HashMap, vec, cmp};
use graph_to_matrix::bipartite_to_biadjaency;
use metis::VertexSeparator;
use petgraph as px;
use px::stable_graph::{NodeIndex, self,StableUnGraph};
use quizx::{linalg::Mat2, vec_graph::{Graph, BasisElem}, hash_graph::GraphLike, circuit::Circuit};
use roots::{find_root_brent, SimpleConvergency};
use vertex_cover::BiGraph;

use crate::anneal::{ComplementFinder, GeometricSeries};

mod anneal;
mod graph_to_matrix;

fn main() {
    //testing vertex separator of zx-diagram
    let c = Circuit::random()
    .qubits(60)
    .depth(2000)
    .seed(35125)
    .clifford_t(0.1)
    .build();

    // 12 to 8
    // .qubits(30)
    // .depth(20)
    // .seed(541651)
    // .min_weight(2)
    // .max_weight(5)
    // .build();
    

    let mut g: Graph = c.clone().to_graph();
    g.plug_outputs(&vec![BasisElem::Z0;c.num_qubits()]);
    g.plug_inputs(&vec![BasisElem::Z0;c.num_qubits()]);
    
    quizx::simplify::full_simp(&mut g);

    println!("{:?}", g.num_vertices());

    let mut gp = quizx_to_petgraph(&g);
    let mut rng = rand::thread_rng();


    for i in 1..=10{
        let mut finder = ComplementFinder::new(&gp,&mut rng,GeometricSeries::new(0.02,0.001,10000),i);
        finder.run(false);
        gp=finder.graph;
        println!("--------------------------")
    }
    



}

fn rank_decomposition(m: &Mat2) -> (Mat2, Mat2) {
    let mut r = m.clone();
    r.gauss_x(true, m.num_cols(), &mut ());

    let mut pivots = Vec::new();
    let mut row = 0;
    let mut col = 0;
    while row < r.num_rows() && col < r.num_cols() {
        if r[(row, col)] != 0 {
            pivots.push(col);
            row += 1;
        }
        col += 1;
    }

    let mut c = Mat2::new(vec![vec![0; pivots.len()]; m.num_rows()]);
    for (i, p) in pivots.into_iter().enumerate() {
        for j in 0..m.num_rows() {
            c[(j, i)] = m[(j, p)];
        }
    }

    let nonzero = (0..r.num_rows())
        .filter(|&i| (0..r.num_cols()).any(|j| r[(i, j)] != 0))
        .collect::<Vec<_>>();

    let mut f = Mat2::new(vec![vec![0; m.num_cols()]; nonzero.len()]);
    for (i, p) in nonzero.into_iter().enumerate() {
        for j in 0..m.num_cols() {
            f[(i, j)] = r[(p, j)];
        }
    }

    (c, f)
}

#[test]
fn rank_decomposition_test() {
    for _ in 0..1000 {
        let mut m = Mat2::new(vec![vec![0; 10]; 10]);
        for i in 0..10 {
            for j in 0..10 {
                m[(i, j)] = (rand::random::<f32>() < 0.2) as u8;
            }
        }

        let (a, b) = rank_decomposition(&m);
        let k = m.rank();
        assert_eq!(a.num_cols(), k);
        assert_eq!(b.num_rows(), k);
        assert_eq!(a * b, m);
    }
}


fn subgraphs_complements_from_rank_decomposition(c: &Mat2, r: &Mat2)-> Vec<(std::vec::Vec<u8>, std::vec::Vec<u8>)>{

    let rank = c.num_cols();
    let sc: Vec<_> =  (0..rank).map(|x| (getcol(c,x),getrow(r,x))).collect();


    fn getcol(m: &Mat2,c:usize)->Vec<u8>{
        
        let mut col = vec![0u8; m.num_rows()];
        for i in 0..m.num_rows(){
            col[i] = m[i][c]
        }
        col
    }
    fn getrow(m: &Mat2,r:usize)->Vec<u8>{
        
        let mut row = vec![0u8; m.num_cols()];
        for i in 0..m.num_cols(){
            row[i] = m[r][i]
        }
        row
    }

    sc
    
}

#[test]
fn subgraphs_complements_from_rank_decomposition_test(){


    let mut m = Mat2::new(vec![vec![0; 5]; 10]);
    for i in 0..10 {
        for j in 0..5 {
            m[(i, j)] = (rand::random::<f32>() < 0.2) as u8;
        }
    }

    let (a, b) = rank_decomposition(&m);
    println!("{}",a);
    println!("{}",b);
    let scs= subgraphs_complements_from_rank_decomposition(&a, &b);
    for sc in scs{
        let (x,y) = sc;
        println!("{:?} {:?}",x,y);
    }


}


fn get_subgraphs_complement_cover(bg : &vertex_cover::BiGraph<(),()>)-> Vec<Vec<NodeIndex>>{

    let mat = bipartite_to_biadjaency(bg);
    let (c,r) = rank_decomposition(&mat);
    let sgcs_wrong_index = subgraphs_complements_from_rank_decomposition(&c,&r);
    let mut sgcs: Vec<Vec<NodeIndex>> = vec![];

    for sgc_wrong_index in sgcs_wrong_index{

        let mut sgc:Vec<NodeIndex> = vec![];

        for (i,&v) in sgc_wrong_index.0.iter().enumerate() {
            if v == 1 {
                sgc.push(bg.right[i as usize]);
            }
        }

        for (i,&v) in sgc_wrong_index.1.iter().enumerate() {
            if v == 1 {
                sgc.push(bg.left[i as usize]);
            }
        }
        sgcs.push(sgc);
    }
    sgcs
}




fn quizx_to_petgraph(g: &Graph) -> px::stable_graph::StableUnGraph<(), ()> {
    let vmapping = g.vertices().scan(0, |i, s| { *i += 1; Some((s, *i - 1)) }).collect::<HashMap<_,_>>();
    let edges = g.edges().map(|x| (vmapping[&x.0] ,vmapping[&x.1] )).collect::<Vec<_>>();
    let graph = px::stable_graph::StableUnGraph::<(), ()>::from_edges(&edges);
    graph
}

pub fn quizx_vertex_cut_finder(g :&Graph) -> Vec<usize> {

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
    let options = metis::Options::default().max_imbalance(100);
    let bad_index_separator = metis::Graph::new(&graph).vertex_separator(&options).unwrap();

    
    let mut separator = vec![];
    for x in bad_index_separator.cut{

        separator.push(reverse_vmapping[&x.index()]);

    }


    separator
}


fn sep_to_bigraph(g:&StableUnGraph<(),()>,sep: &VertexSeparator<StableUnGraph<(),()>>)->BiGraph<(),()>{

    let mut sep = sep.clone();
    if sep.left.len() < sep.right.len(){
        let temp = sep.right;
        sep.right = sep.left;
        sep.left = temp;
    }
    
    let mut subgraph = g.clone();
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
    

   vertex_cover::BiGraph { graph: subgraph, left: new_left, right: new_right }

}


