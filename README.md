# Travelling Salesman Problem: Simulated Annealing vs Genetic Algorithm

> Solving the Travelling Salesman Problem (TSP) for major Ukrainian cities using two stochastic metaheuristics — **Simulated Annealing (SA)** and a **Genetic Algorithm (GA)** — implemented in *Wolfram Mathematica*.

---

## Table of Contents

1. [Problem Statement](#problem-statement)
2. [Why Stochastic Algorithms?](#why-stochastic-algorithms)
3. [Algorithms](#algorithms)
   - [Simulated Annealing](#simulated-annealing)
   - [Genetic Algorithm](#genetic-algorithm)
4. [Literature](#literature)
5. [Implementation](#implementation)
6. [Input Parameters](#input-parameters)
7. [Results](#results)
8. [Repository Structure](#repository-structure)
9. [How to Run](#how-to-run)

---

## Problem Statement

The **Travelling Salesman Problem (TSP)** is one of the most studied combinatorial optimisation problems. Given a set of $n$ cities, the goal is to find the shortest open path visiting each city exactly once:

$$\min_{x \in \mathcal{P}_n} f(x) = \sum_{i=1}^{n-1} d(x_i,\, x_{i+1})$$

where $\mathcal{P}_n$ is the set of all permutations of $n$ cities and $d(\cdot,\cdot)$ is the road distance between two cities.

This implementation solves the TSP for the **regional capitals of Ukraine** (approx. 25 cities), using road travel distances. The initial random tour and the final optimised routes are shown below.

| Random initial tour | SA final tour | GA final tour |
|:-------------------:|:-------------:|:-------------:|
| ![rand](results/randPermutation.png) | ![sa](results/lastSA.png) | ![ga](results/lastGA.png) |

---

## Why Stochastic Algorithms?

TSP belongs to the class of **NP-hard** problems. For $n$ cities, the number of distinct open tours is $(n-1)!/2$, which grows super-exponentially. For $n = 25$ this is already $\approx 3 \times 10^{23}$ — brute-force enumeration is completely infeasible.

Exact methods (branch-and-bound, dynamic programming via Held–Karp in $O(n^2 2^n)$) remain tractable only for small $n$. For practical problem sizes, **metaheuristic** approaches are preferred:

- They **escape local optima** through randomisation (SA via temperature-controlled acceptance of worse solutions; GA via crossover and mutation introducing genetic diversity).
- They are **solution-space agnostic** — the same framework applies regardless of the specific distance metric or city layout.
- They provide **high-quality approximate solutions** in polynomial-time budget.
- They are straightforward to implement and parallelise.

Both SA and GA have long, well-documented track records on TSP instances of all scales.

---

## Algorithms

### Simulated Annealing

Simulated Annealing is inspired by the annealing process in metallurgy, where a material is slowly cooled to minimise its internal energy. The analogy: temperature $T$ controls the willingness to accept a *worse* solution, allowing escape from local optima.

**Algorithm:**

```
Input:  x⁰ — initial tour,  T₀ — initial temperature,
        α  — cooling factor (0 < α < 1),
        k_max — Metropolis steps per temperature level,
        T_min — stopping temperature

x = x⁰;  T = T₀
while T > T_min:
    for k = 1 … k_max:
        x_new = neighbour(x)          // perturbation operator
        ΔE    = f(x_new) − f(x)
        if Exp(−ΔE / T) ≥ Uniform[0,1]:
            x = x_new                 // Metropolis acceptance
    T = α · T                         // exponential cooling
return x
```

**Metropolis acceptance criterion:**

$$P(\Delta f, T) = \min\!\left\{1,\; e^{-\Delta f / T}\right\}$$

A deteriorating move ($\Delta f > 0$) is accepted with probability $e^{-\Delta f/T}$. At high $T$ almost all moves are accepted (exploration); as $T \to 0$ only improvements are accepted (exploitation).

**Perturbation operators** (neighbours of a tour $x$):

| Operator | Description |
|----------|-------------|
| (a) Adjacent swap | Swap two neighbouring cities in the permutation |
| (b) Random swap   | Swap any two cities in the permutation |
| (c) Segment reversal (2-opt) | Reverse the sub-path between two randomly chosen positions |

---

### Genetic Algorithm

The Genetic Algorithm mimics biological evolution: a **population** of candidate tours undergoes **selection**, **crossover** (recombination), and **mutation** across successive **generations**.

**Algorithm:**

```
Input:  N      — population size,
        N_elit — elite size,
        n_gen  — number of generations,
        P_m    — mutation probability

population = N random tours
evaluate fitness f(x) for each individual

for k = 2 … n_gen:
    new_gen = top-N_elit individuals  // elitism
    while |new_gen| < N:
        p1 = tournament_select(population)
        p2 = tournament_select(population)
        child = OrderCrossover(p1, p2)
        if Uniform[0,1] ≤ P_m:
            child = mutate(child)
        new_gen.append(child)
    population = new_gen
    evaluate fitness for each individual
    if improvement stagnates: break

return individual with highest fitness
```

**Fitness function:**

$$f(x) = \frac{1}{\text{tour length}(x)}$$

**Order Crossover (OX):**

1. Copy a random sub-segment from Parent 1 into the child.
2. Fill remaining positions with cities from Parent 2 in their original order, skipping already-placed cities.

**Example:**
```
Parent 1:  [1, 2, | 3, 4, 5, | 6]
Parent 2:  [3, 1,   5, 6, 2,   4]
                ↓ OX
Child:     [1, 6, | 3, 4, 5, | 2]
```

**Mutation:** same perturbation operators (a), (b), (c) as in SA.

**Selection:** Tournament selection — $k$ random individuals are drawn; the one with the highest fitness becomes a parent.

---

## Literature

1. **Knor M., Tomek L.** — *Optimalizácia 2*, STU Bratislava, 2019.
   [[PDF]](https://www.svf.stuba.sk/buxus/docs/dokumenty/skripta/Knor_M.__Tomek_L._Optimalizacia_2_web.pdf)
   *(Primary course reference; Sections 10.3 and 11 cover SA and GA with pseudocode and TSP application.)*

2. **Kirkpatrick S., Gelatt C. D., Vecchi M. P.** — *Optimization by Simulated Annealing*, Science, 220(4598):671–680, 1983.
   [[DOI]](https://doi.org/10.1126/science.220.4598.671)
   *(Seminal paper introducing SA; includes TSP experiments.)*

3. **Černý V.** — *Thermodynamical Approach to the Traveling Salesman Problem: An Efficient Simulation Algorithm*, Journal of Optimization Theory and Applications, 45(1):41–51, 1985.
   [[DOI]](https://doi.org/10.1007/BF00940812)
   *(Independent discovery of SA, with explicit TSP formulation.)*

4. **Holland J. H.** — *Adaptation in Natural and Artificial Systems*, University of Michigan Press, 1975; MIT Press (2nd ed.), 1992.
   *(Foundational work on genetic algorithms.)*

5. **Goldberg D. E.** — *Genetic Algorithms in Search, Optimization, and Machine Learning*, Addison-Wesley, 1989.
   *(Standard GA textbook; covers crossover operators, selection strategies, and TSP encoding.)*

6. **Davis L.** — *Applying Adaptive Algorithms to Epistatic Domains*, Proceedings of IJCAI, 1985.
   *(Introduced Order Crossover (OX) for permutation-based problems such as TSP.)*

7. **Johnson D. S., McGeoch L. A.** — *The Traveling Salesman Problem: A Case Study in Local Optimization*, in *Local Search in Combinatorial Optimization*, Wiley, 1997.
   *(Comprehensive empirical comparison of local-search and metaheuristic methods on TSP.)*

8. **Boyd S., Vandenberghe L.** — *Convex Optimization*, Cambridge University Press, 2004.
   [[Online]](https://web.stanford.edu/~boyd/cvxbook/)
   *(Background reference for Lagrangian duality and barrier methods discussed in the course.)*

---

## Implementation

The algorithms are implemented in **Wolfram Mathematica** (`.nb` notebooks).

Key implementation details:

- **Distance matrix** `L[i,j]` is precomputed once at startup (road distances via `TravelDistance` / geodetic formula) and reused in every fitness evaluation — avoids recomputing distances inside the inner loop.
- **Tour length** is evaluated in $O(n)$ by summing the appropriate entries of `L`.
- **SA neighbour generation**: all three operators (a)(b)(c) are available; best results obtained with segment reversal (2-opt).
- **GA crossover**: `OrderCrossover[p1, p2]` preserves relative city order from Parent 2 while keeping a contiguous block from Parent 1.
- **Elitism** in GA ensures the best solution is never lost between generations.
- **Initial temperature** $T_0$ for SA is chosen so that approximately 50 % of worsening moves are accepted — estimated by sampling random neighbours from the initial tour.

---

## Input Parameters

### Simulated Annealing

| Parameter | Symbol | Typical value | Description |
|-----------|--------|---------------|-------------|
| Initial temperature | $T_0$ | auto (≈ mean $\Delta f$ of random moves) | Controls initial acceptance of worse solutions |
| Cooling factor | $\alpha$ | `0.995` | Geometric cooling rate: $T_{i+1} = \alpha T_i$ |
| Metropolis steps / temperature | $k_{\max}$ | `500` | Inner-loop iterations before lowering $T$ |
| Minimum temperature | $T_{\min}$ | `0.1` km | Stopping criterion |
| Perturbation operator | — | segment reversal (2-opt) | How neighbours are generated |

### Genetic Algorithm

| Parameter | Symbol | Typical value | Description |
|-----------|--------|---------------|-------------|
| Population size | $N$ | `200` | Number of tours per generation |
| Elite size | $N_{\text{elit}}$ | `10` | Tours copied unchanged to next generation |
| Generations | $n_{\text{gen}}$ | `600` | Maximum number of generations |
| Mutation probability | $P_m$ | `0.05` | Probability of applying a mutation operator |
| Tournament size | $k_t$ | `5` | Candidates drawn per tournament selection |
| Crossover operator | — | Order Crossover (OX) | Recombination strategy |
| Mutation operator | — | segment reversal (2-opt) | Mutation strategy |

---

## Results

All experiments use **road travel distances** between Ukrainian regional capitals.

### Convergence

| SA convergence | GA convergence |
|:--------------:|:--------------:|
| ![SA convergence](results/sa_convergence_550steps_rescale.png) | ![GA convergence](results/ga_convergence.png) |

| Metric | SA | GA |
|--------|----|----|
| Initial tour length | ≈ 42 000 km (random) | ≈ 24 000 km (random pop. best) |
| Final tour length | ≈ 7 000 km | ≈ 6 500 km |
| Iterations to converge | 550 temperature steps | ~250 generations |
| Convergence character | Noisy, gradual descent | Rapid early drop, smooth plateau |

**Observations:**
- GA converges faster (fewer evaluations to reach plateau) and to a slightly shorter final tour.
- SA convergence is noisier due to stochastic acceptance of worse solutions; the noise decreases as temperature drops.
- GA benefits from population-level information sharing (crossover recombines good partial tours); SA explores from a single solution.
- Both significantly outperform the random initial permutation (≈ 6× improvement).

### Optimised routes

| SA final route | GA final route |
|:--------------:|:--------------:|
| ![SA final](results/lastSA.png) | ![GA final](results/lastGA.png) |

GA produces a visually cleaner route with fewer crossing edges, consistent with its lower final tour length.

### Video demonstrations

| Video | Link | Description |
|-------|------|-------------|
| **SA — Ukraine TSP** | [▶ YouTube](https://youtu.be/kA8am6jNm4M) | Simulated Annealing solving TSP on Ukrainian regional capitals — animated tour evolution over 550 cooling steps, converging from ≈ 42 000 km to ≈ 7 000 km |
| **GA — Ukraine TSP** | [▶ YouTube](https://youtu.be/tn9F32EiJSs) | Genetic Algorithm solving TSP on Ukrainian regional capitals — best-individual tour evolution across 600 generations, converging from ≈ 24 000 km to ≈ 6 500 km |

---

## Repository Structure

```
tsp-sa-vs-ga-ukraine/
│
├── README.md
│
├── notebooks/
│   └── cv11_SA_Cities_TravelDist.nb      # Mathematica notebook (SA + GA)
│
└── results/
    ├── randPermutation.png               # Initial random tour
    ├── lastSA.png                        # SA optimised tour
    ├── lastGA.png                        # GA optimised tour
    ├── sa_convergence_550steps_rescale.png  # SA tour-length vs step
    └── ga_convergence.png                # GA best-fitness vs generation
```

---

## How to Run

1. **Requirements:** Wolfram Mathematica 12+ (or Wolfram Engine with a valid licence).

2. **Open the notebook:**
   ```
   File → Open → notebooks/cv11_SA_Cities_TravelDist.nb
   ```

3. **Evaluate all cells** (`Evaluation → Evaluate Notebook`).

4. **Interactive controls** — use the `Manipulate` sliders to:
   - Step through SA / GA iterations frame by frame.
   - Adjust temperature, cooling rate, population size, mutation probability on the fly.

5. **Re-run with different parameters** — modify the parameter block at the top of each algorithm section and re-evaluate.

> The distance matrix is built from Mathematica's built-in `TravelDistance` function and may take 1–2 minutes to compute on first run. Subsequent runs use the cached matrix.

---

*Course: I1-OPT2 Optimalizácia 2 — STU Bratislava, Stavebná fakulta, B-MPM programme.*
