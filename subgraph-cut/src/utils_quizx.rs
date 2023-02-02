use std::collections::HashMap;

use quizx::vec_graph::{GraphLike, EType};
use num::{self, Rational};


// This does not work out the correct scalar factors!!!
pub fn subgraph_complement<G: GraphLike>(g : & G, vertices : & Vec<usize>) -> (G,G){

    let mut g =g.clone();
    
    for &v1 in vertices {
        for &v2 in vertices  {
            if v1==v2 {continue;}
            g.add_edge_smart(v1, v2, EType::H);
        }
    }

    let mut g2 = g.clone();

    for &v in vertices{
        g.add_to_phase(v, Rational::new(1,2));
        g2.add_to_phase(v,Rational::new(-1,2));
    }

    (g,g2)
}

pub fn indices_petgraph_to_quizx<G: GraphLike>(g : & G,indices: & Vec<usize>) -> Vec<usize>{

    let vmapping = g.vertices().scan(0, |i, s| { *i += 1; Some((*i - 1, s)) }).collect::<HashMap<_,_>>();

    indices.iter().map(|x| vmapping[x]).collect()


}