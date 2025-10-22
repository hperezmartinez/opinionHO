# Opinion formation model with higher-order interactions

Associated to the paper "Social polarization promoted by sparse higher-order interactions", doi: https://doi.org/10.48550/arXiv.2507.12325.

This repository provides the necessary codes to reproduce the results shown in the aforementioned paper:
- <code>opinions_HO.c</code> receives a given combination of parameters $\lambda$ and $\beta$, and generates a file that includes results for a fixed number of runs starting from random initial conditions. Saved results include parameter combination, number of run, number of steps performed, mean opinion and standard deviation, and exposure defined as in Eq. 4.
- After generation of such files for multiple parameter combinations, code <code>joiner.c</code> provides a unified txt file for a given network structure.


# Higher-Order Network Generators

This repository also provides Python implementations to generate **higher-order networks** with both pairwise and triangular interactions.  
Two models are included:

1. **Hypergraph with maximum inter-order hyperedge overlap**
2. **RSC (Random Simplicial Complex)**

Both models output two text files:
- one containing **pairwise links** (edges)
- one containing **triangular interactions** (3-node hyperedges)

For generating hypergraphs with variable inter-order hyperedge overlap as the ones utilized in Supplementary Fig. 2 we refer to the repository https://github.com/santiagolaot/hyperedge-overlap
