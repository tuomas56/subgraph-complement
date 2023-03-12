use std::{collections::{HashSet, HashMap}, path::PathBuf};
use rand::seq::IteratorRandom;
use zx::graph::{GraphLike, EType, VType};
use serde::Serialize;
use clap::Parser;

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
    pub fitness: f32,
    rng: R,
    temperature: T,
    cut: bool,
    alpha: f32,
    beta: f32
}

impl<G: GraphLike, R: rand::Rng, T: Iterator<Item = f32>> ComplementFinder<G, R, T> {
    pub fn new(graph: &G, rng: R, temperature: T, cut: bool, alpha: f32, beta: f32) -> Self {
        let mut finder = ComplementFinder {
            graph: graph.clone(),
            current: HashSet::new(),
            fitness: 0.0,
            rng, temperature, cut, alpha, beta
        };
        finder.fitness = finder.fitness();
        finder
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

    fn fitness(&self) -> f32 {
        let (edges, vertices) = if self.cut {
            (self.graph.num_edges(), self.graph.num_vertices())
        } else {
            (self.graph.num_edges() + self.current.len(), self.graph.num_vertices() + (self.current.len() > 0) as usize)
        };

        self.alpha * vertices as f32 + self.beta * edges as f32
    }

    fn step(&mut self, temp: f32) {
        let node = self.graph
            .vertices()
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

    pub fn run(&mut self, quiet: bool) {
        let mut step = 0;
        let original = self.fitness;
        while let Some(temp) = self.temperature.next() {
            if !quiet && step % 100000 == 0 {
                println!(
                    "    step = {:?}, temp = {:.2?}, fitness = {:?}, ratio = {:.2?}", 
                    step, temp, self.fitness, self.fitness as f32 / original as f32
                );
            }
            self.step(temp);
            step += 1;
        }

        if !quiet {
            println!(
                "  final: fitness = {:?}, ratio = {:.2?}", 
                self.fitness, self.fitness as f32 / original as f32
            );
        }
    }

    pub fn terms(&self) -> Vec<G> {
        if self.current.is_empty() {
            return vec![self.graph.clone()]
        }

        if self.cut {
            let mut a = self.graph.clone();
            let mut b = self.graph.clone();
            for &n in &self.current {
                a.add_to_phase(n, (1, 2).into());
                b.add_to_phase(n, (-1, 2).into());
            }
            b.scalar_mut().mul_phase((1, 2).into());
            vec![a, b]
        } else {
            let mut a = self.graph.clone();
            let z = a.add_vertex_with_phase(VType::Z, (1, 2).into());
            for &n in &self.current {
                a.add_to_phase(n, (1, 2).into());
                a.add_edge_smart(n, z, EType::H);
            }
            vec![a]
        }
    }
}

#[derive(Serialize)]
struct JsonGraph {
    phases: HashMap<usize, (isize, isize)>,
    edges: Vec<(usize, usize)>,
    scalar: (f64, f64),
    terms: usize
}

impl JsonGraph {
    fn new(g: &impl GraphLike, terms: usize) -> Self {
        let mut phases = HashMap::new();
        let mut edges = Vec::new();

        for n in g.vertices() {
            phases.insert(n, g.phase(n).into());
        }

        for (a, b, _) in g.edges() {
            edges.push((a, b));
        }

        let scalar = g.scalar().float_value();
        let scalar = (scalar.re, scalar.im);

        JsonGraph { phases, edges, scalar, terms }
    }
}

#[derive(Parser)]
struct Args {
    #[clap(required = true, help = "Input circuits to simplify and sparsify")]
    input: Vec<String>,
    #[clap(short, long, default_value_t = 3, help = "Number of rounds of sparsification")]
    rounds: usize,
    #[clap(short, long, default_value_t = 1000000, help = "Number of annealing steps")]
    steps: usize,
    #[clap(short = 'M', long, default_value_t = 3000.0, help = "Starting temperature for annealing")]
    max_temp: f32,
    #[clap(short = 'm', long, default_value_t = 0.1, help = "End temperature for annealing")]
    min_temp: f32,
    #[clap(short, long, help = "Don't cut the complement vertices")]
    no_cut: bool,
    #[clap(short, long, help = "Keep sparsifying forever until cost stops decreasing")]
    infinite: bool,
    #[clap(short, long, default_value_t = 0.0, help = "Weight of vertex count in fitness function")]
    alpha: f32,
    #[clap(short, long, default_value_t = 1.0, help = "Weight of edge count in fitness function")]
    beta: f32
}

fn cost_model(g: &impl GraphLike, alpha: f32, beta: f32) -> f32 {
    g.num_vertices() as f32 * alpha + g.num_edges() as f32 * beta
}

fn main() {
    let args = Args::parse();
    for path in args.input {
        println!("processing: {}", path);
        let c = zx::circuit::Circuit::from_file(&path).unwrap();

        let mut g: zx::hash_graph::Graph = c.to_graph();
        g.plug_outputs(&vec![zx::graph::BasisElem::Z0; c.num_qubits()]);
        g.plug_inputs(&vec![zx::graph::BasisElem::Z0; c.num_qubits()]);

        zx::simplify::full_simp(&mut g);
        let density = 2.0 * g.num_edges() as f32 / (g.num_vertices() as f32 * (g.num_vertices() as f32 - 1.0));
        println!("  vertices = {:?} edges = {:?} density = {:?} cost = {:?}", g.num_vertices(), g.num_edges(), density, cost_model(&g, args.alpha, args.beta));

        let mut rng = rand::thread_rng();
        let mut a = g.clone();
        let mut overall_ratio = 1.0;
        let mut do_round = || {
            let mut finder = ComplementFinder::new(
                &a, &mut rng, 
                GeometricSeries::new(args.max_temp, args.min_temp, args.steps), 
                !args.no_cut, args.alpha, args.beta
            );
            finder.run(false);
            let aa = finder.terms().pop().unwrap();
            let ratio = cost_model(&aa, args.alpha, args.beta) / cost_model(&a, args.alpha, args.beta);
            let done = cost_model(&aa, args.alpha, args.beta) >= cost_model(&a, args.alpha, args.beta);
            if !done { a = aa };
            let density = 2.0 * a.num_edges() as f32 / (a.num_vertices() as f32 * (a.num_vertices() as f32 - 1.0));
            println!("  vertices = {:?} edges = {:?} density = {:?} cost = {:?}", a.num_vertices(), a.num_edges(), density, cost_model(&a, args.alpha, args.beta));  
            overall_ratio *= ratio;
            done
        };

        let mut rounds = 0;
        if args.infinite {
            loop {
                if do_round() {
                    break
                } else {
                    rounds += 1;
                }
            }
        } else {
            for _ in 0..args.rounds {
                do_round();
                rounds += 1;
            }
        }

        println!("  overall cost ratio: {}", overall_ratio);
        
        let path = PathBuf::from(path);
        let aname = path.with_file_name(if args.infinite {
            format!("{}-sparsified-inf.json", path.file_stem().unwrap().to_string_lossy())
        } else {
            format!("{}-sparsified-{}.json", path.file_stem().unwrap().to_string_lossy(), args.rounds)
        }); 
        let gname = path.with_file_name(format!("{}-simplified.json", path.file_stem().unwrap().to_string_lossy()));
        let mut f = std::fs::File::create(&aname).unwrap();
        serde_json::to_writer(&mut f, &JsonGraph::new(&a, if args.no_cut { 1 } else {
            1 << rounds
        })).unwrap();
        let mut f = std::fs::File::create(&gname).unwrap();
        serde_json::to_writer(&mut f, &JsonGraph::new(&g, 0)).unwrap();
        println!("  wrote: `{}`", aname.display());
        println!("  wrote: `{}`", gname.display());
    }
}
