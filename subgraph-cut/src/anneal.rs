use std::collections::HashSet;
use rand::seq::IteratorRandom;
use petgraph as px;

use std::cmp;
use metis::VertexSeparator;
use px::stable_graph::{NodeIndex,StableUnGraph};
use quizx::{vec_graph::Graph, hash_graph::GraphLike};
use roots::{find_root_brent, SimpleConvergency};

use crate::{sep_to_bigraph, get_subgraphs_complement_cover};




#[derive(Debug, Clone)]
pub struct GeometricSeries {
    current: f32,
    scale: f32,
    steps: usize
}

impl GeometricSeries {
    pub fn new(start: f32, stop: f32, steps: usize) -> Self {
        GeometricSeries {
            current: start.ln(), steps,
            scale: (stop.ln() - start.ln()) / steps as f32
        }
    }
}

impl Iterator for GeometricSeries {
    type Item = f32;

    fn next(&mut self) -> Option<f32> {
        if self.steps == 0 {
            None
        } else {
            let out = self.current.exp();
            self.current += self.scale;
            self.steps -= 1;
            Some(out)
        }
    }
}

pub struct ComplementFinder<R: rand::Rng, T: Iterator<Item = f32>> {
    pub graph: StableUnGraph<(), ()>,
    pub current: HashSet<px::graph::NodeIndex>,
    pub fitness: f32,
    rng: R,
    temperature: T,
    depth: usize
}

impl<R: rand::Rng, T: Iterator<Item = f32>> ComplementFinder<R, T> {
    pub fn new(graph: &StableUnGraph<(),()>, mut rng: R, temperature: T,depth:usize) -> Self {
        let complement = graph.node_indices()
                .choose_multiple(&mut rng, graph.node_count() / (2*depth));
        let mut finder = ComplementFinder {
            graph: graph.clone(),
            current: HashSet::new(),
            fitness: 0.0,
            rng, temperature,depth
        };
        for c in complement {
            finder.toggle_node(c);
        }
        finder.fitness = finder.fitness();
        finder
    }

    fn toggle_node(&mut self, node: px::graph::NodeIndex) {
        let mut present = false;
        for &other in &self.current {
            if other == node {
                present = true;
                continue
            }

            if let Some(e) = self.graph.find_edge(node, other) {
                self.graph.remove_edge(e);
            } else {
                self.graph.add_edge(node, other, ());
            }
        }

        if !present {
            self.current.insert(node);
        } else {
            self.current.remove(&node);
        }
    }

    fn step(&mut self, temp: f32) {
        let node = self.graph
            .node_indices()
            .choose(&mut self.rng)
            .unwrap();

        self.toggle_node(node);
        let new_fitness = self.fitness();

        if new_fitness < self.fitness {
            self.fitness = new_fitness;
            return
        }

        let prob = ((new_fitness - self.fitness) as f32 / temp).exp().recip();
        if self.rng.gen::<f32>() < prob {
            self.fitness = new_fitness;
        } else {
            self.toggle_node(node);
        }
    }

    fn fitness(&self) -> f32 {
        let sep = metis::Graph::new(&self.graph)
        .vertex_separator(&metis::Options::default()
            .max_imbalance(10))
        .unwrap();

        let bg = sep_to_bigraph(&self.graph,&sep);

        let sgcs = get_subgraphs_complement_cover(&bg);

        alpha_subgraph_complements(&self.graph,&sep,&sgcs,1)
        
    }

    pub fn run(&mut self, quiet: bool) {
        let mut step = 0;
        let original = self.fitness;
        while let Some(temp) = self.temperature.next() {
            if !quiet && step % 1000 == 0 {
                println!(
                    "step = {:?}, temp = {:.2?}, fitness = {:?}, ratio = {:.2?}", 
                    step, temp, self.fitness, self.fitness as f32 / original as f32
                );
            }
            self.step(temp);
            step += 1;
        }

        if !quiet {
            println!(
                "final: fitness = {:?}, ratio = {:.2?}", 
                self.fitness, self.fitness as f32 / original as f32
            );
        }
    }
}


pub fn crossing_edges(graph: &vertex_cover::BiGraph<(), ()>)-> usize{
    
    let edges = graph.graph.edge_indices();
    let mut crossing_count = 0;

    for e in edges{
        let (a,b) = graph.graph.edge_endpoints(e).unwrap();
        if  (graph.left.contains(&a) && graph.right.contains(&b)) || (graph.left.contains(&b) && graph.right.contains(&a)) {
            crossing_count +=1;
        }
    }

    crossing_count


}

fn alpha_subgraph_complements(g:&StableUnGraph<(),()>,sep :& VertexSeparator<StableUnGraph<(),()>>, subgraph_complements: &Vec<Vec<NodeIndex>>,nb_initials_subgraph_complements:usize) -> f32{
    // TODO : Take tcounts into account!
    
        let num_diag = subgraph_complements.len() + nb_initials_subgraph_complements;
        let n = g.node_count();
        let (d1,d2) = (n-cmp::min(sep.left.len(),sep.right.len())+ sep.cut.len(), n-cmp::max(sep.left.len(),sep.right.len()) );
        let (a,b) = (cmp::min(d1,d2),cmp::max(d1,d2));
    
    
        let f = |x:f32| {
            a as f32 * x.ln() - (num_diag as f32)/1.442695 + (1.0 +x.powi(a as i32 -b as i32)).ln_1p() 
            // the constant is lg(e)
        };

        let t = find_root_brent(1.0, 2.0, f, &mut SimpleConvergency{eps:0.001,max_iter:100}).unwrap();
        t.log2()
    }
    