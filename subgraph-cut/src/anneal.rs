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
    max_nb_complement: usize,
    rng: R,
    temperature: T,
    depth: usize,
    best_fitness: f32,
    best_current: HashSet<px::graph::NodeIndex>,
    best_graph: StableUnGraph<(), ()>,
    current_cut: (VertexSeparator<StableUnGraph<(), ()>>, Vec<Vec<NodeIndex>>),
    best_cut: (VertexSeparator<StableUnGraph<(), ()>>, Vec<Vec<NodeIndex>>),
    // vertex_count:usize
}

impl<R: rand::Rng, T: Iterator<Item = f32>> ComplementFinder<R, T> {
    pub fn new(
        graph: &StableUnGraph<(), ()>,
        mut rng: R,
        temperature: T,
        depth: usize,
        imbalance: usize,
        max_nb_complement: usize,
    ) -> Self {
        let complement = graph
             .node_indices()
             .choose_multiple(&mut rng, 3);
            // .choose_multiple(&mut rng, graph.node_count() / 2);
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
            max_nb_complement,
            current_cut: (VertexSeparator{left : Vec::new(), cut :Vec::new(), right : Vec::new()},Vec::new()),
            // BiGraph { graph: StableUnGraph::default(), left: Vec::new(), right: Vec::new() },
            best_cut: (VertexSeparator{left : Vec::new(), cut :Vec::new(), right : Vec::new()},Vec::new()),
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
            self.best_cut = self.current_cut.clone();
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

    fn fitness(&mut self) -> f32 {
        let (sep, sgcs) = self.complement_cover();
        let barrier = log_barrier(sgcs.len(), self.max_nb_complement)/200.0;

        let temp = improved_alpha_subgraph_complements(&self.graph, &sep, &sgcs, self.depth) + barrier;

        self.current_cut = (sep,sgcs);

        temp
    }
    fn old_fitness(&self) -> f32 {
        let (sep, sgcs) = self.complement_cover();
        alpha_subgraph_complements(&self.graph, &sep, &sgcs, self.depth)
    }

    fn complements(&self) -> usize {
        self.current_cut.1.len() + self.depth
    }

    fn vertex_cut(&self) -> usize {
        self.depth + self.current_cut.0.cut.len()
    }

    pub fn run(&mut self, quiet: bool) {
        let mut step = 0;
        let mut prob = 1.0;
        let original = self.fitness;
        let initial_cut = self.complement_cover();

        let initial_cut_alpha = vertex_cut_alpha(&initial_cut,self.depth as isize -1);

        println!("Inital alpha with vertex cut : {} with {} vertices", initial_cut_alpha, initial_cut.1.len());



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
        self.current_cut=self.best_cut.clone();
        self.fitness = self.best_fitness;

        if !quiet {
            println!(
                "final: alpha = {:?}, nb_sugraph_complement : {} of which {} are cuts. \nTo do the same cut with vertex cut + rounds of sugraph complement : {}",
                improved_alpha_subgraph_complements(&self.graph, &self.best_cut.0, &self.best_cut.1, self.depth),
                self.complements(),
                vertex_cut_in_subgraph_complement_cut(&self.best_cut),
                self.vertex_cut()
            );
            println!("the best fitness is {}", self.best_fitness);

            let smaller = self.current_cut.0.left.len().min( self.current_cut.0.right.len());
            let bigger = self.current_cut.0.left.len().max( self.current_cut.0.right.len());
            println!(
                "The split is {} vs {}",
                smaller as isize + self.best_cut.1.len() as isize - vertex_cut_in_subgraph_complement_cut(&self.best_cut) as isize, 
                bigger
            );
        }
    }
    pub fn solution_found(&self) -> HashSet<px::graph::NodeIndex>{
        self.current.clone()
    }
}

fn alpha_subgraph_complements(
    g: &StableUnGraph<(), ()>,
    sep: &VertexSeparator<StableUnGraph<(), ()>>,
    subgraph_complements: &Vec<Vec<NodeIndex>>,
    nb_initials_subgraph_complements: usize
) -> f32 {
    // TODO : Take tcounts into account!

    // log2 of number of diagrams
    let num_diag = subgraph_complements.len() + nb_initials_subgraph_complements;
    let n = g.node_count();

    // d1 and d2 are the number of vertices in each sub instance
    let ( d1, d2) = (
        (cmp::min(sep.left.len(), sep.right.len()) + sep.cut.len()),
        cmp::max(sep.left.len(), sep.right.len()),
    );

    // a and b are the number of vertices removed
    let (a, b) = (n - cmp::max(d1, d2), n- cmp::min(d1, d2));

    // println!("a : {} b : {}, num_diag : {} ",a,b, num_diag);

    let f = |x: f32| {
        a as f32 * x.ln() - (num_diag as f32) / 1.442695
            + (1.0 + x.powi(a as i32 - b as i32)).ln_1p()
        // the constant is lg(e)
    };

    let t = find_root_brent(
        1.0,
        10.0,
        f,
        &mut SimpleConvergency {
            eps: 0.001,
            max_iter: 100,
        },
    )
    .unwrap();
    t.log2()
}



fn improved_alpha_subgraph_complements(
    g: &StableUnGraph<(), ()>,
    sep: &VertexSeparator<StableUnGraph<(), ()>>,
    subgraph_complements: &Vec<Vec<NodeIndex>>,
    nb_initials_subgraph_complements: usize,
) -> f32 {
    // TODO : Take tcounts into account!

    // log2 of number of diagrams
    let num_diag = subgraph_complements.len() + nb_initials_subgraph_complements;
    let n = g.node_count();


    // d1 and d2 are the number of vertices in each sub instance
    let ( mut d1, mut d2) = (
        (cmp::min(sep.left.len(), sep.right.len()) + sep.cut.len()),
        cmp::max(sep.left.len(), sep.right.len()),
    );

    for sgc in subgraph_complements{

        let subgraph = sgc.iter().collect::<HashSet<_>>();
        let intersect_cut = sep.cut.iter().collect::<HashSet<_>>().intersection(&subgraph).count();

        if intersect_cut == sgc.len() - 1 {
            d2 -= 1;
        } else if intersect_cut == 1  {
            d1 -= 1 ;
        }
    }

    // a and b are the number of vertices removed
    let (a, b) = (n - cmp::max(d1, d2), n- cmp::min(d1, d2));


    let f = |x: f32| {
        a as f32 * x.ln() - (num_diag as f32) / 1.442695
            + (1.0 + x.powi(a as i32 - b as i32)).ln_1p()
        // the constant is lg(e)
    };

    let t = find_root_brent(
        1.0,
        10.0,
        f,
        &mut SimpleConvergency {
            eps: 0.001,
            max_iter: 100,
        },
    )
    .unwrap();
    t.log2()
}

fn vertex_cut_in_subgraph_complement_cut((sep,sgcs):&(VertexSeparator<StableUnGraph<(), ()>>, Vec<Vec<NodeIndex>>)) -> usize{

    let mut nb = 0;
    for sgc in sgcs{

        let subgraph = sgc.iter().collect::<HashSet<_>>();
        let intersect_cut = sep.cut.iter().collect::<HashSet<_>>().intersection(&subgraph).count();

        if intersect_cut == sgc.len() - 1 {
            nb += 1;
        } else if intersect_cut == 1  {
            nb += 1 ;
        }
    }
    nb
}

fn vertex_cut_alpha((sep,_):&(VertexSeparator<StableUnGraph<(), ()>>, Vec<Vec<NodeIndex>>),depth:isize) -> f32{


    let n = (sep.left.len()+sep.right.len()+sep.cut.len()) as isize;
    let (nb_d1, nb_d2) = (cmp::max(sep.right.len(), sep.left.len()), cmp::min(sep.right.len(), sep.left.len()));

    // println!("nb_d1 : {} nb_d2 : {}, num_diag : {} ",nb_d1,nb_d2, sep.cut.len() as isize +depth);
    let (a,b) = (n  - nb_d1 as isize, n- nb_d2 as isize);
    // println!("a : {} b : {}, num_diag : {} ",a,b, sep.cut.len() as isize +depth);

    let f = |x: f32| {
        a as f32 * x.ln() - (sep.cut.len() as isize + depth) as f32 / 1.442695
            + (1.0 + x.powi(a as i32 - b as i32)).ln_1p()
        // the constant is lg(e)
    };

    let t = find_root_brent(
        1.0,
        10.0,
        f,
        &mut SimpleConvergency {
            eps: 0.001,
            max_iter: 100,
        },
    )
    .unwrap();
    t.log2()


}
fn log_barrier(x:usize,away_from:usize)->f32{
    if x >= away_from {
        2000.0 *(x as f32)/(away_from as f32)
    }else {
        -((away_from as f32 - x as f32)/away_from as f32).log2()
    }
}