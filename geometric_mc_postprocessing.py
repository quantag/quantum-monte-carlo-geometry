import math
import numpy as np
import matplotlib.pyplot as plt

from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator


# -----------------------
# Scenario model
# -----------------------

def ry_for_prob(p):
    return 2.0 * math.asin(math.sqrt(p))


def build_A_gate(n, p_list):
    qc = QuantumCircuit(n, name="A")
    for i in range(n):
        qc.ry(ry_for_prob(p_list[i]), i)
    return qc.to_gate()


def build_Sf_mark_all_ones(n):
    qc = QuantumCircuit(n + 1, name="Sf")
    data = [qc.qubits[i] for i in range(n)]
    anc = qc.qubits[n]
    qc.mcx(data, anc)
    return qc.to_gate()


def build_S0_flip_all_zeros(n):
    qc = QuantumCircuit(n + 1, name="S0")
    data = [qc.qubits[i] for i in range(n)]
    anc = qc.qubits[n]
    for q in data:
        qc.x(q)
    qc.mcx(data, anc)
    for q in data:
        qc.x(q)
    return qc.to_gate()


def build_Q_gate(n, A_gate):
    Sf = build_Sf_mark_all_ones(n)
    S0 = build_S0_flip_all_zeros(n)

    qc = QuantumCircuit(n + 1, name="Q")
    data = list(range(n))
    qc.append(Sf, qc.qubits)
    qc.append(A_gate.inverse(), data)
    qc.append(S0, qc.qubits)
    qc.append(A_gate, data)
    return qc.to_gate()


def build_circuit(n, p_list, k):
    qc = QuantumCircuit(n + 1, n)
    data = [qc.qubits[i] for i in range(n)]
    anc = qc.qubits[n]

    qc.x(anc)
    qc.h(anc)

    A = build_A_gate(n, p_list)
    qc.append(A, data)

    Q = build_Q_gate(n, A)
    for _ in range(k):
        qc.append(Q, qc.qubits)

    qc.measure(data, list(range(n)))
    return qc


# -----------------------
# Estimation utilities
# -----------------------

def prob_all_ones_from_counts(counts, n):
    total = 0
    hit = 0
    target = "1" * n
    for key, v in counts.items():
        total += v
        if key[::-1] == target:
            hit += v
    return hit / float(total)


def run_prob(sim, circ_t, shots, n):
    res = sim.run(circ_t, shots=shots).result()
    counts = res.get_counts(0)
    return prob_all_ones_from_counts(counts, n)


def allocate_shots(total_shots, k_list):
    w = np.array([(2 * k + 1) ** 2 for k in k_list], dtype=float)
    w /= np.sum(w)
    shots = np.maximum(1, np.floor(w * total_shots).astype(int))
    diff = int(total_shots - np.sum(shots))
    i = 0
    while diff != 0:
        j = i % len(shots)
        if diff > 0:
            shots[j] += 1
            diff -= 1
        else:
            if shots[j] > 1:
                shots[j] -= 1
                diff += 1
        i += 1
    return shots.tolist()


def mle_theta_sin2(ks, pks, Nk):
    thetas = np.linspace(1e-10, (math.pi / 2) - 1e-10, 20001)
    ll = np.zeros_like(thetas)

    for k, pk, N in zip(ks, pks, Nk):
        a = 2 * k + 1
        p = np.sin(a * thetas) ** 2
        p = np.clip(p, 1e-12, 1.0 - 1e-12)
        pk = float(np.clip(pk, 1e-6, 1.0 - 1e-6))
        ll += (pk * N) * np.log(p) + ((1.0 - pk) * N) * np.log(1.0 - p)

    return float(thetas[np.argmax(ll)])


# -----------------------
# Benchmark + plots
# -----------------------

def main():
    n = 10
    p_list = [0.8, 0.8, 0.8, 0.8, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4]
    p_true = np.prod(p_list)

    k_list = [0, 1, 2, 4, 8]
    eps = 1e-4
    conf = 0.95
    trials = 200

    shot_grid = [
        200, 400, 800, 1600, 3200, 5000,
        8000, 12000, 20000, 50000,
        120000, 300000, 800000
    ]

    sim = AerSimulator()

    circ0_t = transpile(build_circuit(n, p_list, 0), sim, optimization_level=1)
    circ_t_map = {k: transpile(build_circuit(n, p_list, k), sim, optimization_level=1)
                  for k in set([0] + k_list)}

    naive_err = []
    geo_err = []
    naive_succ = []
    geo_succ = []

    for S in shot_grid:
        e_naive = []
        e_geo = []
        ok_n = 0
        ok_g = 0

        for _ in range(trials):
            p_n = run_prob(sim, circ0_t, S, n)
            e_naive.append(abs(p_n - p_true))
            if abs(p_n - p_true) <= eps:
                ok_n += 1

            Sk = allocate_shots(S, k_list)
            pks = [run_prob(sim, circ_t_map[k], sk, n) for k, sk in zip(k_list, Sk)]
            theta_hat = mle_theta_sin2(k_list, pks, Sk)
            p_g = math.sin(theta_hat) ** 2
            e_geo.append(abs(p_g - p_true))
            if abs(p_g - p_true) <= eps:
                ok_g += 1

        naive_err.append(np.mean(e_naive))
        geo_err.append(np.mean(e_geo))
        naive_succ.append(ok_n / trials)
        geo_succ.append(ok_g / trials)

    # -----------------------
    # Plot 1: error vs shots
    # -----------------------
    plt.figure(figsize=(6, 4))
    plt.loglog(shot_grid, naive_err, "o-", label="Naive Monte Carlo")
    plt.loglog(shot_grid, geo_err, "o-", label="Geometric AE")
    plt.axhline(eps, linestyle="--", color="gray")
    plt.xlabel("Total shots")
    plt.ylabel("Mean absolute error")
    plt.legend()
    plt.tight_layout()
    plt.savefig("mc_error_vs_shots.png", dpi=300)

    # -----------------------
    # Plot 2: success probability
    # -----------------------
    plt.figure(figsize=(6, 4))
    plt.semilogx(shot_grid, naive_succ, "o-", label="Naive Monte Carlo")
    plt.semilogx(shot_grid, geo_succ, "o-", label="Geometric AE")
    plt.axhline(conf, linestyle="--", color="gray")
    plt.xlabel("Total shots")
    plt.ylabel("Success probability")
    plt.legend()
    plt.tight_layout()
    plt.savefig("mc_success_vs_shots.png", dpi=300)

    print("Plots written:")
    print("  mc_error_vs_shots.png")
    print("  mc_success_vs_shots.png")
    print("p_true =", p_true)
    print("eps =", eps, "conf =", conf)


if __name__ == "__main__":
    main()
