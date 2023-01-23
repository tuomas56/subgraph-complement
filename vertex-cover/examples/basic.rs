use std::collections::HashSet;

fn main() {
    let mut total = std::time::Duration::ZERO;
    for _ in 0..10 {
        let mut rng = rand::thread_rng();
        let graph = petgraph_gen::random_gnp_graph(&mut rng, 1000, 0.01);
        let nodes = graph.node_indices().collect::<Vec<_>>();
        let (a, b) = nodes.split_at(nodes.len() / 2);
        let bigraph = vertex_cover::BiGraph { graph, left: a.to_vec(), right: b.to_vec() };
        let before = std::time::Instant::now();
        let cover = bigraph.min_vertex_cover();
        total += before.elapsed();

        let mut marked = HashSet::new();
        for &v in &cover {
            for e in bigraph.graph.neighbors(v) {
                if (bigraph.left.contains(&v) && bigraph.right.contains(&e)) || (bigraph.left.contains(&e) && bigraph.right.contains(&v)) {
                    marked.insert((v.min(e), v.max(e)));
                } 
            }
        }

        let mut all = 0;
        for &i in a {
            for &j in b {
                if bigraph.graph.contains_edge(i, j) {
                    all += 1;
                }    
            } 
        }

        println!("cover size = {:?}, covered edges = {:?}, total edges = {:?}", cover.len(), marked.len(), all);
        assert_eq!(marked.len(), all);
    }
    println!("avg time = {:?}", total/10);
}
