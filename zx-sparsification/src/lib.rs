use std::collections::HashSet;
use rand::seq::IteratorRandom;
use zx::graph::{GraphLike, EType};
// use zx::hash_graph::Graph;

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

pub struct ComplementFinder<G: GraphLike, R: rand::Rng, T: Iterator<Item = f32>> {
    pub graph: G,
    pub current: HashSet<usize>,
    pub fitness: usize,
    rng: R,
    temperature: T
}

impl<G: GraphLike, R: rand::Rng, T: Iterator<Item = f32>> ComplementFinder<G, R, T> {
    pub fn new(graph: &G, rng: R, temperature: T) -> Self {
        ComplementFinder {
            graph: graph.clone(),
            current: HashSet::new(),
            fitness: graph.num_edges(),
            rng, temperature
        }
    }

    fn toggle_node(&mut self, node: usize) {
        let mut present = false;
        for &other in &self.current {
            if other == node {
                present = true;
                continue
            }

            self.graph.add_edge_smart(node, other, EType::H);
        }

        if !present {
            self.current.insert(node);
        } else {
            self.current.remove(&node);
        }
    }

    fn step(&mut self, temp: f32) {
        let node = self.graph
            .vertices()
            .choose(&mut self.rng)
            .unwrap();

        self.toggle_node(node);
        let new_fitness = self.graph.num_edges();

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

pub struct PivotFinder<G: GraphLike, R: rand::Rng, T: Iterator<Item = f32>> {
    pub graph: G,
    pub left: HashSet<usize>,
    pub right: HashSet<usize>,
    pub fitness: usize,
    rng: R,
    temperature: T
}

impl<G: GraphLike, R: rand::Rng, T: Iterator<Item = f32>> PivotFinder<G, R, T> {
    pub fn new(graph: &G, rng: R, temperature: T) -> Self {
        PivotFinder {
            graph: graph.clone(),
            left: HashSet::new(),
            right: HashSet::new(),
            fitness: graph.num_edges(),
            rng, temperature
        }
    }

    fn toggle_node(&mut self, node: usize, left: bool) {
        let (this_side, other_side) = if left {
            (&mut self.left, &self.right)
        } else {
            (&mut self.right, &self.left)
        };

        for &other in other_side {
            self.graph.add_edge_smart(node, other, EType::H);
        }

        if !this_side.contains(&node) {
            this_side.insert(node);
        } else {
            this_side.remove(&node);
        }
    }

    fn step(&mut self, temp: f32) {
        let left = self.rng.gen_bool(0.5);

        let node = loop {
            let prop = self.graph
                .vertices()
                .choose(&mut self.rng)
                .unwrap();

            if left && !self.right.contains(&prop) {
                break prop
            } else if !left && !self.left.contains(&prop) {
                break prop
            }
        };
        

        self.toggle_node(node, left);
        let new_fitness = self.graph.num_edges();

        if new_fitness < self.fitness {
            self.fitness = new_fitness;
            return
        }

        let prob = ((new_fitness - self.fitness) as f32 / temp).exp().recip();
        if self.rng.gen::<f32>() < prob {
            self.fitness = new_fitness;
        } else {
            self.toggle_node(node, left);
        }
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

#[derive(Debug)]
pub enum SparsifierMove {
    Complement(HashSet<usize>),
    Pivot(HashSet<usize>, HashSet<usize>)
}

pub struct Sparsifier<G: GraphLike, R: rand::Rng, T: Iterator<Item = f32> + Clone> {
    pub graph: G,
    pub moves: Vec<SparsifierMove>,
    rng: R,
    temperature: T,
    alpha: f32,
    beta: f32
}

impl<G: GraphLike, R: rand::Rng, T: Iterator<Item = f32> + Clone> Sparsifier<G, R, T > {
    pub fn new(graph: &G, temperature: T, rng: R, alpha: f32, beta: f32) -> Self {
        Sparsifier {
            graph: graph.clone(),
            moves: Vec::new(),
            temperature, rng, alpha, beta
        }
    }

    fn cost(&self, graph: &G) -> f32 {
        graph.num_vertices() as f32 * self.alpha + graph.num_edges() as f32 * self.beta
    }

    pub fn run(&mut self, quiet: bool) {
        let original = self.cost(&self.graph);
        let originale = self.graph.num_edges();
        loop {
            if !quiet {
                println!("=> cost: {:.2?}, edges: {:?}", self.cost(&self.graph), self.graph.num_edges());
            }

            let mut cfinder = ComplementFinder::new(
                &self.graph, &mut self.rng, self.temperature.clone()
            );
            cfinder.run(quiet);
            let cgraph = cfinder.graph;
            let cloc = cfinder.current;
            let complement_cost = 1.0 + self.cost(&cgraph);

            let mut pfinder = PivotFinder::new(
                &self.graph, &mut self.rng, self.temperature.clone()
            );
            pfinder.run(quiet);
            let pgraph = pfinder.graph;
            let ploc = (pfinder.left, pfinder.right);
            let pivot_cost = 2.0 + self.cost(&pgraph);

            if complement_cost.min(pivot_cost) >= self.cost(&self.graph) {
                if !quiet {
                    println!("=> no improvement found, stopping");
                }

                break
            }

            if complement_cost < pivot_cost {
                if !quiet {
                    println!("=> doing complement on {} vertices", cloc.len());
                }

                self.graph = cgraph;
                self.moves.push(SparsifierMove::Complement(cloc));
            } else {
                if !quiet {
                    println!("=> doing pivot between {} and {} vertices", ploc.0.len(), ploc.1.len());
                }

                self.graph = pgraph;
                self.moves.push(SparsifierMove::Pivot(ploc.0, ploc.1));
            }
        }

        if !quiet {
            println!(
                "=> cost ratio: {:?}, edge ratio: {:?}, move list:", 
                self.cost(&self.graph) / original, 
                self.graph.num_edges() as f32 / originale as f32
            );
            for m in &self.moves {
                match m {
                    SparsifierMove::Complement(c) => println!(
                        "  - complement on {:?} vertices", c.len()
                    ),
                    SparsifierMove::Pivot(a, b) => println!(
                        "  - pivot between {} and {} vertices", a.len(), b.len()
                    )
                }
            }
        }
    }
}

// fn gnp(n: usize, p: f32) -> Graph {
//     let mut graph = Graph::new();
    
//     let nodes = (0..n)
//         .map(|_| graph.add_vertex(VType::Z))
//         .collect::<Vec<_>>();

//     for &a in &nodes {
//         for &b in &nodes {
//             if a < b && rand::random::<f32>() < p {
//                 graph.add_edge_smart(a, b, EType::H);
//             }
//         }
//     }

//     graph
// }

// fn main() {
//     let graph = gnp(100, 0.8);

//     let mut finder = ComplementFinder::new(
//         &graph, 
//         rand::thread_rng(), 
//         GeometricSeries::new(3000.0, 0.1, 200000)
//     );

//     finder.run(false);

//     // let mut sparsifier = Sparsifier::new(
//     //     &graph, 
//     //     GeometricSeries::new(3000.0, 1.0, 100000),
//     //     rand::thread_rng(),
//     //     0.25, 0.03
//     // );

//     // let before = std::time::Instant::now();
//     // sparsifier.run(false);
//     // println!("{:?} elapsed", before.elapsed());
// }
