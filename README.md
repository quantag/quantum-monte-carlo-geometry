Geometric Monte Carlo Post-Processing
====================================

This repository accompanies the article:

  "Geometric Quantum Amplitude Estimation for Monte Carlo Post-Processing"

It contains a single, self-contained Python script that reproduces all numerical
results and figures reported in the paper.

The purpose of this repository is reproducibility: running the provided script
generates the data and plots used to demonstrate a reduction in measurement
cost for a concrete Monte Carlo post-processing task.

----------------------------------------------------------------------
Problem description
----------------------------------------------------------------------

We consider a Monte Carlo scenario model implemented as a quantum circuit.

- The circuit prepares a probability distribution over binary variables
  representing independent Bernoulli risk factors.
- A rare event is defined as a specific predicate on the sampled bitstrings
  (in this benchmark: all variables simultaneously equal to one).
- The goal is to estimate the probability p of this rare event.

This task is representative of Monte Carlo post-processing problems such as:
- tail-risk estimation,
- rare-event probability estimation,
- constraint-violation probability estimation.

The benchmark compares two estimation strategies applied to the *same*
scenario generator:

1) Naive Monte Carlo sampling:
   - Sample the circuit output repeatedly.
   - Estimate p as the empirical frequency of the event.

2) Geometric amplitude-estimation-based post-processing:
   - Use amplitude amplification with several circuit depths.
   - Infer p from the geometric structure of the amplified probabilities
     using maximum-likelihood estimation.

The comparison is performed at fixed accuracy and confidence requirements.

----------------------------------------------------------------------
Benchmark criterion
----------------------------------------------------------------------

For both methods we measure:

- Absolute estimation error |p_hat - p_true|
- Success probability: fraction of trials satisfying |p_hat - p_true| <= eps

The benchmark reports the *minimal total number of shots* required to achieve:

- Absolute error <= eps
- With success probability >= conf

In the results reported in the paper:
- eps = 1e-4
- conf = 0.95
- Each data point is averaged over 200 independent trials

----------------------------------------------------------------------
Repository contents
----------------------------------------------------------------------

geometric_mc_postprocessing.py
  The main script. It:
  - constructs the quantum circuits,
  - runs both estimation methods,
  - evaluates accuracy and success probability,
  - generates publication-quality PNG figures.

LICENSE
  MIT license.

----------------------------------------------------------------------
Requirements
----------------------------------------------------------------------

Python 3.9+ is recommended.

Install dependencies with:

  pip install qiskit qiskit-aer numpy matplotlib

No additional packages are required.

----------------------------------------------------------------------
How to run
----------------------------------------------------------------------

Simply execute:

  python geometric_mc_postprocessing.py

The script will:

1) Run the Monte Carlo post-processing benchmark.
2) Generate two PNG files in the current directory:
     - mc_error_vs_shots.png
     - mc_success_vs_shots.png
3) Print the true probability p_true and the benchmark parameters.

The runtime depends on the shot budget and number of trials; on a typical
workstation it completes within a few minutes.

----------------------------------------------------------------------
Output figures
----------------------------------------------------------------------

mc_error_vs_shots.png
  Mean absolute error versus total number of shots for naive Monte Carlo
  sampling and the geometric estimator.

mc_success_vs_shots.png
  Probability of achieving the target absolute error as a function of
  the total number of shots. The 95% confidence threshold is indicated.

These figures are the ones used in the accompanying article.

----------------------------------------------------------------------
Scope and limitations
----------------------------------------------------------------------

This repository demonstrates a *post-processing advantage* for probability
estimation under the following conditions:

- The quantity of interest is a probability (or rare-event probability).
- The event can be implemented as a clean phase oracle.
- The amplitude amplification geometry is preserved.

The code does not claim a universal speedup for all quantum algorithms,
nor does it replace general-purpose solvers. It demonstrates a concrete,
reproducible advantage for a specific and practically relevant estimation task.

----------------------------------------------------------------------
Reproducibility
----------------------------------------------------------------------

All numerical results and figures in the article were generated using
this script without manual intervention.

