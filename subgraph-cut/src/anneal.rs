use crate::bigraph::BiGraph;
use metis::VertexSeparator;
use petgraph as px;
use px::stable_graph::{NodeIndex, StableUnGraph};
use rand::seq::IteratorRandom;
use roots::{find_root_brent, SimpleConvergency};
use std::cmp;
use std::collections::HashSet;

#[derive(Debug, Clone)]
pub struct GeometricSeries {
    current: f32,
    scale: f32,
    steps: usize,
}

impl GeometricSeries {
    pub fn new(start: f32, stop: f32, steps: usize) -> Self {
        GeometricSeries {
            current: start.ln(),
            steps,
            scale: (stop.ln() - start.ln()) / steps as f32,
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
    imbalance: usize,
    rng: R,
    temperature: T,
    depth: usize,
    best_fitness: f32,
    best_current: HashSet<px::graph::NodeIndex>,
    best_graph: StableUnGraph<(), ()>,
}

impl<R: rand::Rng, T: Iterator<Item = f32>> ComplementFinder<R, T> {
    pub fn new(
        graph: &StableUnGraph<(), ()>,
        mut rng: R,
        temperature: T,
        depth: usize,
        imbalance: usize,
    ) -> Self {
        let complement = graph
            .node_indices()
            .choose_multiple(&mut rng, graph.node_count() / (2 * depth));
        let mut finder = ComplementFinder {
            graph: graph.clone(),
            current: HashSet::new(),
            fitness: 0.0,
            best_fitness: 0.0,
            best_graph: graph.clone(),
            best_current: HashSet::new(),
            rng,
            temperature,
            depth,
            imbalance,
        };
        for c in complement {
            finder.toggle_node(c);
        }
        finder.fitness = finder.fitness();
        finder.best_fitness = finder.fitness;
        finder
    }

    fn toggle_node(&mut self, node: px::graph::NodeIndex) {
        let mut present = false;
        for &other in &self.current {
            if other == node {
                present = true;
                continue;
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

    fn step(&mut self, temp: f32) -> f32 {
        let node = self.graph.node_indices().choose(&mut self.rng).unwrap();

        self.toggle_node(node);
        let new_fitness = self.fitness();

        if new_fitness < self.fitness {
            self.fitness = new_fitness;
            return 1.0;
        }

        let prob = ((new_fitness - self.fitness) as f32 / temp).exp().recip();
        if self.rng.gen::<f32>() < prob {
            self.fitness = new_fitness;
        } else {
            self.toggle_node(node);
        }

        if self.fitness < self.best_fitness {
            self.best_graph = self.graph.clone();
            self.best_current = self.current.clone();
            self.best_fitness = self.fitness;
        }

        prob.min(1.0)
    }

    fn complement_cover(&self) -> (VertexSeparator<StableUnGraph<(), ()>>, Vec<Vec<NodeIndex>>) {
        let sep = metis::Graph::new(&self.graph)
            .vertex_separator(&metis::Options::default().max_imbalance(self.imbalance))
            .unwrap();

        let bg = BiGraph::from_sep(&self.graph, &sep);

        let sgcs = bg.complement_cover();

        (sep, sgcs)
    }

    fn fitness(&self) -> f32 {
        let (sep, sgcs) = self.complement_cover();
        alpha_subgraph_complements(&self.graph, &sep, &sgcs, 1)
    }

    fn complements(&self) -> usize {
        let (_, sgcs) = self.complement_cover();
        sgcs.len() + self.depth
    }

    fn vertex_cut(&self) -> usize {
        let (sep, _) = self.complement_cover();
        self.depth + sep.cut.len()
    }

    pub fn run(&mut self, quiet: bool) {
        let mut step = 0;
        let mut prob = 1.0;
        let original = self.fitness;
        while let Some(temp) = self.temperature.next() {
            if !quiet && step % 1000 == 0 {
                println!(
                    "step = {:?}, temp = {:.2?}, fitness = {:?}, ratio = {:.2?}, prob = {:.2?}, complements = {}, cut size = {}", 
                    step, temp, self.fitness, self.fitness as f32 / original as f32, prob, self.complements(), self.vertex_cut()
                );
            }
            prob = self.step(temp);
            step += 1;
        }

        self.graph = self.best_graph.clone();
        self.current = self.best_current.clone();
        self.fitness = self.best_fitness;

        if !quiet {
            println!(
                "final: fitness = {:?}, ratio = {:.2?}",
                self.fitness,
                self.fitness as f32 / original as f32
            );
        }
    }
}

fn alpha_subgraph_complements(
    g: &StableUnGraph<(), ()>,
    sep: &VertexSeparator<StableUnGraph<(), ()>>,
    subgraph_complements: &Vec<Vec<NodeIndex>>,
    nb_initials_subgraph_complements: usize,
) -> f32 {
    // TODO : Take tcounts into account!

    let num_diag = subgraph_complements.len() + nb_initials_subgraph_complements;
    let n = g.node_count();
    let (d1, d2) = (
        n - cmp::min(sep.left.len(), sep.right.len()) + sep.cut.len(),
        n - cmp::max(sep.left.len(), sep.right.len()),
    );
    let (a, b) = (cmp::min(d1, d2), cmp::max(d1, d2));

    let f = |x: f32| {
        a as f32 * x.ln() - (num_diag as f32) / 1.442695
            + (1.0 + x.powi(a as i32 - b as i32)).ln_1p()
        // the constant is lg(e)
    };

    let t = find_root_brent(
        1.0,
        2.0,
        f,
        &mut SimpleConvergency {
            eps: 0.001,
            max_iter: 100,
        },
    )
    .unwrap();
    t.log2()
}
