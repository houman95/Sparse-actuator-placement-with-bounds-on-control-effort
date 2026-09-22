# Sparse Actuator Placement Under Control-Effort Constraints

This repository contains MATLAB code for selecting a small set of actuator locations in linear dynamical systems. The objective is not only to make the system controllable. It is also to retain strong reachability while limiting the required control energy.

The project implements and numerically adapts the approximation algorithms in Tzoumas *et al.*, ["Minimal Actuator Placement with Bounds on Control Effort"](https://doi.org/10.1109/TCNS.2015.2444031). It was completed as a Sparse Control course project at the Department of Electrical Engineering, Sharif University of Technology.

**Authors:** Houman Asgari and Milad Farjadnasab

## Problem setting

Consider the continuous-time linear time-invariant system

$$
\dot{x}(t)=Ax(t)+B_\Delta u(t),
$$

where each state is a possible actuator location. An actuator set $\Delta\subseteq\{1,\ldots,n\}$ is represented by

$$
B_\Delta=\mathrm{diag}(\delta), \qquad
\delta_i=
\begin{cases}
1, & i\in\Delta,\\
0, & i\notin\Delta.
\end{cases}
$$

For a given set $\Delta$, the controllability Gramian can be written as

$$
W_\Delta=\sum_{i\in\Delta} W_i,
$$

where $W_i$ is the Gramian contribution from actuating state $i$ alone. If $W_\Delta$ is nonsingular, the unit-energy reachable set is

$$
\mathcal{R}_\Delta =
\lbrace x \mid x^\top W_\Delta^{-1}x \leq 1 \rbrace.
$$

Its volume is proportional to $\sqrt{\det(W_\Delta)}$. A larger log-determinant therefore means that more state-space directions can be reached with limited input energy. This gives the equivalent control-effort metric

$$
J(\Delta)=\log\det\left(W_\Delta^{-1}\right)
=-\log\det(W_\Delta).
$$

The project addresses two NP-hard placement problems.

### 1. Minimum-cardinality placement

Find the smallest actuator set that satisfies a prescribed control-effort bound:

$$
\min_{\Delta}|\Delta|
\quad\text{subject to}\quad
\log\det\left(W_\Delta^{-1}\right)\leq E.
$$

This formulation asks how few states must be actuated to obtain a desired reachable-set volume.

### 2. Placement with a fixed actuator budget

For a maximum of $r$ actuators, find a set with low control effort:

$$
\min_{\Delta}\log\det\left(W_\Delta^{-1}\right)
\quad\text{subject to}\quad
|\Delta|\leq r.
$$

Equivalently, this maximizes the volume of the unit-energy reachable set under a fixed actuator budget.

## Implemented approach

The code follows a three-level procedure.

1. **Greedy actuator selection (`alg_2.m`)**  
   Starting from an empty set, the algorithm adds the actuator that gives the largest decrease in the regularized objective

   $$
   \log\det\left(\widetilde W_\Delta+\epsilon I\right)^{-1}.
   $$

   The regularization makes the objective well-defined before the selected actuator set gives a full-rank Gramian.

2. **Regularization search (`alg_3.m`)**  
   A bisection search selects a sufficiently small $\epsilon$. This controls the gap between the regularized problem and the original log-determinant objective.

3. **Fixed-budget placement (`alg_4.m`)**  
   A second bisection search varies the admissible effort bound until the returned set contains at most $r$ actuators. The implementation refines the bisection tolerance when needed to recover the requested cardinality.

## Numerical implementation

Several changes were made to improve numerical behavior in MATLAB:

- **QR-based log-determinants:** Direct determinant evaluation can underflow for large or poorly conditioned Gramians. The code uses the diagonal of the $R$ factor from a QR decomposition to evaluate the log-determinant. In one large-scale run from the original project, this reduced the reported computation time from about 20 seconds with an eigenvalue calculation to about 9 seconds.
- **Regularization of singular Gramians:** Adding $\epsilon I$ permits the greedy search to compare actuator sets before full controllability has been reached.
- **Adaptive bisection accuracy:** The tolerance is reduced only when the current search does not return the requested number of actuators.
- **Support for unstable models:** The large-scale experiment uses `CtrGram.m`, a Gramian routine for unstable LTI systems. The routine requires that $A$ have no poles on the imaginary axis.
- **Floating-point safeguards:** When the required $\epsilon$ falls below the representable range, the implementation uses MATLAB's `realmin` as a numerical floor.

## Experiments

### Five-state chain

A small chain-structured system is used to check the algorithms across several control-effort bounds and actuator budgets. The selected sets show the expected diminishing returns. Adding the first few actuators produces a much larger increase in reachable-set volume than adding later actuators.

### 105-state epidemic-spreading network

The algorithms are also tested on a 105-state model built from the Pajek `GD99c` social-network graph. Node-dependent spreading and recovery parameters define the state matrix. This choice can produce an unstable system, so the unstable-system Gramian routine is used.

The experiment evaluates the log-determinant metric over different actuator budgets. It also compares the selected five-actuator set with 1,000 uniformly sampled five-actuator sets. The algorithm's set achieved a larger log-determinant than every sampled baseline. This is empirical evidence of strong placement quality, but it is not a proof of global optimality.

The scripts generate:

- the log-determinant of the selected Gramian versus the actuator budget;
- selected actuator sets for different effort bounds; and
- a histogram comparing the five-actuator solution with random placements.

## Core files

| File | Purpose |
| --- | --- |
| `alg_2.m` | Greedy solution of the regularized bound-constrained problem |
| `alg_3.m` | Bisection search for the regularization parameter $\epsilon$ |
| `alg_4.m` | Fixed-budget actuator placement through bisection on the effort bound |
| [`CtrGram.m`](https://www.mathworks.com/matlabcentral/fileexchange/38328-controllability-gramian-for-unstable-systems) | Third-party controllability-Gramian routine for unstable LTI systems |
| `GD99_c.mat` | Adjacency data for the 105-node network experiment |
| `alpha.mat`, `beta.mat` | Node-dependent parameters for the spreading model |

## Requirements and basic use

- MATLAB
- Control System Toolbox
- [`CtrGram.m`](https://www.mathworks.com/matlabcentral/fileexchange/38328-controllability-gramian-for-unstable-systems) for the unstable-network experiment
- The supplied `.mat` files for the 105-state case study

The core routines expect an array `W(:,:,i)` containing the individual Gramian contribution of every candidate actuator and a matrix `Wn` containing their sum. Typical calls are:

```matlab
% Minimum-cardinality placement for a prescribed effort bound
[delta, W_delta, epsilon] = alg_3(W, Wn, E_tilde, c, a0);

% Placement with a budget of r actuators
[delta, W_delta] = alg_4(W, Wn, r, a0_prime, a0_prime_min);

selected_actuators = find(delta);
```

Run the experiment script after adding the repository folder and the required data files to the MATLAB path.

## Scope and attribution

This repository is an implementation and experimental study of published actuator-placement methods. The project contribution consists of the MATLAB implementation, numerical adaptations, large-scale case study, and empirical evaluation. The approximation framework and its theoretical guarantees are due to the cited work.

## Reference

V. Tzoumas, M. A. Rahimian, G. J. Pappas, and A. Jadbabaie, "Minimal Actuator Placement With Bounds on Control Effort," *IEEE Transactions on Control of Network Systems*, vol. 3, no. 1, pp. 67–78, 2016. DOI: [10.1109/TCNS.2015.2444031](https://doi.org/10.1109/TCNS.2015.2444031).

```bibtex
@article{tzoumas2016minimal,
  author  = {Vasileios Tzoumas and Mohammad Amin Rahimian and
             George J. Pappas and Ali Jadbabaie},
  title   = {Minimal Actuator Placement With Bounds on Control Effort},
  journal = {IEEE Transactions on Control of Network Systems},
  volume  = {3},
  number  = {1},
  pages   = {67--78},
  year    = {2016},
  doi     = {10.1109/TCNS.2015.2444031}
}
```
