#[allow(nonstandard_style)] 
pub mod sys;

pub trait GraphLike {
    type Vertex: Clone;

    fn vertices(&self) -> Vec<Self::Vertex>;
    fn has_edge(&self, a: Self::Vertex, b: Self::Vertex) -> bool;
}

#[cfg(feature = "petgraph")]
impl<N, E, Ty: petgraph::EdgeType> GraphLike for petgraph::Graph<N, E, Ty> {
    type Vertex = petgraph::graph::NodeIndex;
    
    fn vertices(&self) -> Vec<Self::Vertex> {
        self.node_indices().collect()
    }

    fn has_edge(&self, a: Self::Vertex, b: Self::Vertex) -> bool {
        self.contains_edge(a, b)
    }
}

#[cfg(feature = "petgraph")]
impl<N, E, Ty: petgraph::EdgeType> GraphLike for petgraph::stable_graph::StableGraph<N, E, Ty> {
    type Vertex = petgraph::stable_graph::NodeIndex;
    
    fn vertices(&self) -> Vec<Self::Vertex> {
        self.node_indices().collect()
    }

    fn has_edge(&self, a: Self::Vertex, b: Self::Vertex) -> bool {
        self.contains_edge(a, b)
    }
}

pub struct Graph<G: GraphLike> {
    map: Vec<G::Vertex>,
    xadj: Vec<sys::idx_t>,
    adjncy: Vec<sys::idx_t>,
    nvtxs: sys::idx_t
}

impl<G: GraphLike> Graph<G> {
    pub fn new(graph: &G) -> Self {
        let map: Vec<_> = graph.vertices();
        let nvtxs = map.len() as sys::idx_t;
        let mut xadj = vec![0];
        let mut adjncy = Vec::new();
        
        for a in &map {
            for (j, b) in map.iter().enumerate() {
                if graph.has_edge(a.clone(), b.clone()) || graph.has_edge(b.clone(), a.clone()) {
                    adjncy.push(j as sys::idx_t);
                }
            }
            
            xadj.push(adjncy.len() as sys::idx_t);
        }

        Graph { map, xadj, adjncy, nvtxs }
    }

    pub fn vertex_separator(&self, options: &Options) -> Result<VertexSeparator<G>, Error> {
        let mut sepsize = 0;
        let mut part = vec![0; self.nvtxs as usize];
        let error = unsafe {
             sys::METIS_ComputeVertexSeparator(
                &self.nvtxs as *const _ as *mut _,
                self.xadj.as_ptr() as *mut _,
                self.adjncy.as_ptr() as *mut _,
                std::ptr::null_mut(),
                options.options.as_ptr() as *mut _, 
                &mut sepsize as *mut _,
                part.as_mut_ptr()
            )
        };

        if error == sys::rstatus_et_METIS_OK {
            let mut partition = VertexSeparator {
                left: Vec::new(), right: Vec::new(), cut: Vec::new()
            };

            for i in 0..self.nvtxs as usize {
                match part[i] {
                    0 => partition.left.push(self.map[i].clone()),
                    1 => partition.right.push(self.map[i].clone()),
                    2 => partition.cut.push(self.map[i].clone()),
                    _ => ()
                }
            }

            Ok(partition)
        } else {
            match error {
                sys::rstatus_et_METIS_ERROR_MEMORY => Err(Error::Memory),
                sys::rstatus_et_METIS_ERROR_INPUT => Err(Error::Input),
                _ => Err(Error::Other)
            }
        }
    }
}

#[derive(Debug)]
pub struct VertexSeparator<G: GraphLike> {
    pub left: Vec<G::Vertex>,
    pub right: Vec<G::Vertex>,
    pub cut: Vec<G::Vertex>
}

pub struct Options {
    options: [sys::idx_t; sys::METIS_NOPTIONS as usize]
}

impl Default for Options {
    fn default() -> Self {
        let mut options = [0; sys::METIS_NOPTIONS as usize];
        unsafe {
            sys::METIS_SetDefaultOptions(options.as_mut_ptr());
        }
        Options { options }
    }
}

impl Options {
    pub fn max_imbalance(mut self, factor: usize) -> Self {
        self.options[sys::moptions_et_METIS_OPTION_UFACTOR as usize] = factor as sys::idx_t;
        self
    }
}

#[derive(Debug)]
pub enum Error {
    Input,
    Memory,
    Other
}
