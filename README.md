# FastTreeID
This is an implementation of the algorithm from the paper 'Faster Generic Identification in Tree-Shaped Structural Causal Models' (2025) by Briefs and Bläser, based on the approach of Gupta and Bläser in 'Identification for Tree-Shaped Structural Causal Models in Polynomial Time' (2024).

It determines, for a structural causal model (SCM) whose directed edges form a tree, whether each parameter is unidentifiable, 1-identifiable or 2-identifiable (other cases cannot occur), using a randomized algorithm with provable running time $O(n^3 log^2 n)$.

This repository contains both an R package (fasttreeid) and a C++ version (in cpp-standalone/). The C++ source files shared by the R package and the C++ implementation are in src/.

## The R package fasttreeid
**Usage**: `fasttreeid_identify(bidirected, directed, seed = NULL, prime = NULL)`

**Arguments**:
- `bidirected`: A 2-column matrix in which each row specifies a bidirected edge (1-indexed)
- `directed`: The directed parents $p_i$ of nodes i from 2 to n ($1 \leq p_i < i$)
- `seed`: The seed to use as a decimal string. Picks a random seed by default. Must fit into an unsigned 64-bit integer
- `prime`: The prime to use as a decimal string. Picks a random prime by default. If no prime is given, the package `gmp` must be installed. The prime must be between $2^{58}$ and $2^{59}$

**Return value**: Returns a nested named list describing the identification result:
- `identification`: A list with $n$ entries for all nodes $i$ from 1 to $n$. For each node, it contains a list describing the result for the node. For all nodes, this list contains an entry `identifiability` that can take the values 0 (unidentifiable), 1 (1-identifiable) and 2 (2-identifiable). For all nodes that are not unidentifiable, this list contains an entry `type` that can take the following values:
   - `path`: This means that this node can be identified by missing bidirected edges of rank 2 to a node identified by a fraction or a cycle (Lemma 4 in the paper by Gupta and Bläser). In this case, the list additionally contains an entry `nodes` that describes the path to this node. The first node in the list is the current node and the type of the last node in the list is `fraction` or `cycle`
   - `fraction`: This means that this node can be identified by a missing bidirected edge to the root (Section 2.1 in the paper by Gupta and Bläser). In this case, the list additionally contains two entries `numerator` and `denominator`, each of the format `{"what": "sigma", "i": i, "j": j}`
   - `cycle`: This means that this node can be identified by a cycle of missing bidirected edges of rank 2 (Definition 10 in the paper by Gupta and Bläser). In this case, the list additionally contains an entry `nodes` that describes the cycle. The first and last node in the list are the current node. If the identifiability of the current node is 1, the list further contains an entry `reason` explaining why the identifiability is 1. It can take the following values:
      - `a_is_zero`: This means that $a = 0$ in the equation for some $\lambda_{i,j}$ (case 2 in Lemma 8 in the paper by Gupta and Bläser). In this case, the list additionally contains an entry `reason_edge` of the format `{"what": "lambda", "i": i, "j": j}`
      - `discriminant_is_zero`: This means that the discriminant is 0 (case 1 in Lemma 8 in the paper by Gupta and Bläser)
      - `only_one_option`: This means that the equation for a missing bidirected edge { $i, j$ } is only satisfied by one of the options (last step in Algorithm 2 in the paper by Gupta and Bläser). In this case, the list additionally contains an entry `reason_edge` of the format `{"what": "missing_bidirected", "i": i, "j": j}`
- `seed`: The seed used for identification for reproducibility
- `prime`: The prime used for identification for reproducibility

**Example**
```
library(fasttreeid)

### This defines the example in Figure 1 in the paper by Gupta and Bläser
bidirected <- matrix(c(2, 3), byrow = TRUE, ncol = 2)
directed <- c(1, 2)

fasttreeid_identify(bidirected, directed, seed="1909956540882291621", prime = "499562279898941057")
```

## The C++ implementation
The C++ program can be compiled by running `make` in cpp-standalone/. This requires OpenSSL to be installed (install by `apt install libssl-dev` in Ubuntu for example). Alternatively, using `make with-gmp`, the program can be compiled with GMP instead of OpenSSL. By `make no-prime-check`, the program can be compiled without OpenSSL and GMP, but then, a 59-bit prime has to be given manually every call. Otherwise, the program can generate a random prime using GMP or OpenSSL.

The executable is called `FastTreeID`. It reads the input from stdin and writes to stdout.

**Input format**

$n$  
$p_1$ $\dots$ $p_{n-1}$  
$m$  
$u_0~v_0$  
$\dots$  
$u_{m-1}$ $v_{m-1}$  

where $n$ is the number of nodes; $m$ is the number of edges;
the nodes are numbered $0, …, n-1$;
$p_1, …, p_{n-1}$ are the directed parents of $1, …, n-1$;
and ${u_0, v_0}, …, {u_{m-1},v_{m-1}}$ are the bidirected edges

**Options**
```
   --seed SEED        Seed the random with SEED. Generates a random seed by default.
   --prime PRIME      Use PRIME as the prime. PRIME should be between 2^58 and 2^59.
                      Selects a random prime in this range by default.
   --minimal          Only output n-1 integers i_1 … i_{n-1}, i_j∈{0,1,2}
                      0 means unidentifiable, 1 means 1-identifiable, 2 means 2-identifiable
   --help             Display this help and exit
```

**Example**
```
> ./FastTreeID < tests/small/small_figure1_in.txt
Identification completed with seed = 16301772305941192412 and prime = 441644412843059291
Result:
1 is 1-identifiable as (σ_0,1)/(σ_0,0)
2 is 1-identifiable as (σ_0,2)/(σ_0,1)

```

## General notes
The algorithm is randomized. The error probability is below $4.1e-6$ for $n = 200$ and below $0.0026$ for $n = 1000$ (the estimated running time for $n = 1000$ is more than two hours). These estimates are very conservative worst-case estimates, and the algorithm can be run multiple times with different seeds (and primes) to decrease the error probability exponentially.

If no seed is given, a seed is generated randomly. If the prime is not given, using the same seed (on the same system) twice will generate the same prime twice. Whether or not a prime is given doesn't influence the other random computations.

Note that the input is 1-indexed in R and 0-indexed in C++.

## References
Yasmine Briefs and Markus Bläser. Faster generic identification in tree-shaped structural causal models. Advances in Neural Information Processing Systems 29: Annual Conference on Neural Information Processing Systems 2025, NeurIPS 2025, December 2-7, 2025, San Diego, USA

Aaryan Gupta and Markus Bläser. Identification for tree-shaped structural causal models in polynomial time. Thirty-Eight AAAI Conference on Artificial Intelligence, AAAI 2024, Thirty-Sixth Conference on Innovative Applications of Artificial Intelligence, IAAI 2024, Fourteenth Symposium on Educational Advances in Artificial Intelligence, EAAI 2014, February 20-27, 2024, Vancouver, Canada

Benito van der Zander, Marcel Wienöbst, Markus Bläser, and Maciej Liskiewicz. Identification in tree-shaped linear structural causal models. International Conference on Artificial Intelligence and Statistics, AISTATS 2022, March 28-30, 2022, Virtual Event

