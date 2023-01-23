fn main() {
    let mut rng = rand::thread_rng();
    let graph = petgraph_gen::random_gnp_graph::<_, petgraph::Undirected, _>(&mut rng, 100, 0.05);
    let mgraph = metis::Graph::new(&graph);
    let options = metis::Options::default();
    let sep = mgraph.vertex_separator(&options).unwrap();
    println!("{}", sep.cut.len());
}