use std::{path::PathBuf, fs::File, io::Write};
use rand::Rng;
use clap::{Parser, ValueEnum};

#[derive(Debug, ValueEnum, Copy, Clone)]
enum CircuitType {
    PauliGadget,
    CliffordT
}

#[derive(Parser)]
struct Args {
    #[clap(help = "Output directory for the circuits")]
    output_dir: PathBuf,
    #[clap(short, long, default_value = "pauli-gadget", help = "Type of Circuit to generate")]
    r#type: CircuitType,
    #[clap(short, long, help = "Starting seed for generation, circuits are given sequential seeds")]
    seed: Option<u64>,
    #[clap(short, long, default_value_t = 1, help = "Number of circuits to generate")]
    count: usize,
    #[clap(short, long, default_value_t = 20, help = "Number of qubits for each circuit")]
    qubits: usize,
    #[clap(short, long, help = "Depth for each circuit (default is qubits^2 for Clifford+T, 2*qubits for Pauli gadgets)")]
    depth: Option<usize>,
    #[clap(short = 'm', long, default_value_t = 2, help = "Minimum weight of Pauli gadgets")]
    min_weight: usize,
    #[clap(short = 'M', long, default_value_t = 5, help = "Maximum weight of Pauli gadgets")]
    max_weight: usize,
    #[clap(short = 'f', long, default_value_t = 0.1, help = "Fraction of gates that are T-gates")]
    fraction_t: f32
}

fn main() {
    let args = Args::parse();

    let seed = args.seed.unwrap_or(rand::thread_rng().gen::<u64>());
    for i in 0..args.count {
        let (name, c) = match args.r#type {
            CircuitType::CliffordT => {
                let depth = args.depth.unwrap_or(args.qubits * args.qubits);
                let name = format!("circ-{}-{}-cliffordt-{}-{}.qasm", args.qubits, depth, args.fraction_t.to_string().replace('.', "_"), seed + i as u64);
                (name, zx::circuit::Circuit::random()
                    .qubits(args.qubits)
                    .depth(depth)
                    .clifford_t(args.fraction_t)
                    .seed(seed + i as u64)
                    .build()
                    .to_basic_gates())
            },
            CircuitType::PauliGadget => {
                let depth = args.depth.unwrap_or(args.qubits * 2);
                let name = format!("circ-{}-{}-pauligadget-{}-{}-{}.qasm", args.qubits, depth, args.min_weight, args.max_weight, seed + i as u64);
                (name, zx::circuit::Circuit::random_pauli_gadget()
                    .qubits(args.qubits)
                    .depth(depth)
                    .min_weight(args.min_weight)
                    .max_weight(args.max_weight)
                    .seed(seed + i as u64)
                    .build()
                    .to_basic_gates())
            }
        };

        let path = args.output_dir.join(name);
        let mut file = File::create(&path).unwrap();
        write!(&mut file, "{}", c.to_qasm()).unwrap();
        println!("{}/{}: wrote `{}`", i + 1, args.count, path.display());
    }
}