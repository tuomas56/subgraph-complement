use std::collections::HashSet;
use rand::seq::IteratorRandom;
use petgraph as px;

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
    pub graph: vertex_cover::BiGraph<(), ()>,
    pub current: HashSet<px::graph::NodeIndex>,
    pub fitness: usize,
    rng: R,
    temperature: T
}

impl<R: rand::Rng, T: Iterator<Item = f32>> ComplementFinder<R, T> {
    pub fn new(graph: &vertex_cover::BiGraph<(), ()>, rng: R, temperature: T) -> Self {
        let mut finder = ComplementFinder {
            graph: graph.clone(),
            current: HashSet::new(),
            fitness: 0,
            rng, temperature
        };
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

            if let Some(e) = self.graph.graph.find_edge(node, other) {
                self.graph.graph.remove_edge(e);
            } else {
                self.graph.graph.add_edge(node, other, ());
            }
        }

        if !present {
            self.current.insert(node);
        } else {
            self.current.remove(&node);
        }
    }

    fn step(&mut self, temp: f32) {
        let node = self.graph.graph
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

    fn fitness(&self) -> usize {
        self.graph.min_vertex_cover().len() * self.graph.graph.edge_count()
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