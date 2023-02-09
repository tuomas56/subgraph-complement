import pyzx as zx
import json
import quimb
import quimb.tensor
import numpy as np
import cotengra
from fractions import Fraction
import glob
import sys
import os

def from_json(j):
    g = zx.Graph()
    nodes = {}
    for node, phase in j["phases"].items():
        node = int(node)
        nodes[node] = g.add_vertex(zx.VertexType.Z, phase=Fraction(*phase))
    for a, b in j["edges"]:
        g.add_edge((nodes[a], nodes[b]), zx.EdgeType.HADAMARD)
    
    return g, j["terms"]


if __name__ == "__main__":
    dir = sys.argv[1]
    for file in glob.glob(f'{dir}/*.qasm'):
        name = file.removesuffix(".qasm")
        stats_file = name + '-stats.json'

        print(f"processing: {name}")
        if os.path.exists(stats_file):
            print("=> already done, skipping")
            continue

        circ = file
        circ_zxcirc = zx.Circuit.from_qasm_file(circ)
        circ_circ = quimb.tensor.Circuit(circ_zxcirc.qubits)
        for gate in circ_zxcirc.gates:
            qubits = []
            if hasattr(gate, 'phase'):
                qubits.append(float(gate.phase) * np.pi)
            for a in ["ctrl1","ctrl2", "control", "target"]:
                if hasattr(gate, a):
                    qubits.append(getattr(gate, a))
            circ_circ.apply_gate(gate.qasm_name, *qubits)
        circ_tn: quimb.tensor.TensorNetwork = circ_circ.amplitude_tn('0' * circ_circ.N)
        circ_tn = circ_tn.full_simplify(seq = 'ADCRS')
        print("=> finding contraction tree for circuit:")
        opt = cotengra.HyperOptimizer(progbar=True, max_time = 90.0, max_repeats=128, parallel=False)
        circ_path = circ_tn.contract(output_inds=(), get='path-info', optimize=opt)

        print("=> finding contraction tree for zx circuit:")
        circ_zx_tn = zx.quimb.to_quimb_tensor(circ_zxcirc.to_graph())
        circ_zx_tn.full_simplify(seq = 'ADCRS')
        opt = cotengra.HyperOptimizer(progbar=True, max_time = 90.0, max_repeats=128, parallel=False)
        circ_zx_path = circ_zx_tn.contract(output_inds=(), get='path-info', optimize=opt)

        simpl = name + '-simplified.json'
        simpl_graph, _ = from_json(json.load(open(simpl, "r")))
        simpl_tn: quimb.tensor.TensorNetwork = zx.quimb.to_quimb_tensor(simpl_graph)
        simpl_tn = simpl_tn.full_simplify(seq = 'ADCRS')
        
        print("=> finding contraction tree for simplified graph:")
        opt = cotengra.HyperOptimizer(progbar=True, max_time = 90.0, max_repeats=128, parallel=False)
        simpl_path = simpl_tn.contract(output_inds=(), get='path-info', optimize=opt)

        stats = {
            "circuit": {
                "flops": float(circ_path.opt_cost.log10()),
                "mem": float(circ_path.largest_intermediate.log10()) / np.log10(2)
            },
            "simplified": {
                "flops": float(simpl_path.opt_cost.log10()),
                "mem": float(simpl_path.largest_intermediate.log10()) / np.log10(2)
            },
            "sparsified": {}
        }

        for sparse in glob.glob(name + '-sparsified-*.json'):
            sparse_graph, sparse_terms = from_json(json.load(open(sparse, "r")))
            sparse_tn: quimb.tensor.TensorNetwork = zx.quimb.to_quimb_tensor(sparse_graph)
            sparse_tn = sparse_tn.full_simplify(seq = 'ADCRS')
            sparse_depth = int(np.log2(sparse_terms))
            print(f"=> finding contraction tree for sparsified graph (depth {sparse_depth}):")
            opt = cotengra.HyperOptimizer(progbar=True, max_time = 90.0, max_repeats=128, parallel=False)
            sparse_path = sparse_tn.contract(output_inds=(), get='path-info', optimize=opt)
        
            stats["sparsified"][sparse_depth] = {
                "flops": float((sparse_terms * sparse_path.opt_cost).log10()),
                "mem": float(sparse_path.largest_intermediate.log10()) / np.log10(2)
            }

            circ_flops_ratio = float(circ_path.opt_cost / (sparse_terms * sparse_path.opt_cost))
            circ_flops_percent = 100.0 * (1.0 - 1.0 / circ_flops_ratio)
            circ_mem_ratio = float(circ_path.largest_intermediate / sparse_path.largest_intermediate)
            circ_mem_percent = 100.0 * (1.0 - 1.0 / circ_mem_ratio)
            simpl_flops_ratio = float(simpl_path.opt_cost / (sparse_terms * sparse_path.opt_cost))
            simpl_flops_percent = 100.0 * (1.0 - 1.0 / simpl_flops_ratio)
            simpl_mem_ratio = float(simpl_path.largest_intermediate / sparse_path.largest_intermediate)
            simpl_mem_percent = 100.0 * (1.0 - 1.0 / simpl_mem_ratio)

            print(f"=> compared to circuit:\n    flops = {circ_flops_ratio:.2}x better, {circ_flops_percent:.2}% saved,\n    mem = {circ_mem_ratio}x better, {circ_mem_percent}% saved")
            print(f"=> compared to simplified:\n    flops = {simpl_flops_ratio:.2}x better, {simpl_flops_percent:.2}% saved,\n    mem = {simpl_mem_ratio}x better, {simpl_mem_percent}% saved")

        json.dump(stats, open(stats_file, 'w'))
