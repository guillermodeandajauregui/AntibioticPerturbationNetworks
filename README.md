# Antibiotic Resistance Network Model

This repository contains code to simulate the emergence and spread of
antibiotic resistance in microbial communities represented as networks.

The model represents microbes as nodes in a graph. Resistance can spread
through contacts between nodes, while antibiotic exposure perturbs the
network by removing susceptible nodes. Simulations are executed across
multiple network topologies and parameter regimes to explore how
antibiotic pressure reshapes microbial interaction networks.

## Repository structure

R/ graph_prep.R sim_resistance_helpers.R sim_resistance_steps.R
simulate_resistance_network.R sweep_helpers.R

scripts/ run_sweep_ER_raw.R run_sweep_WS_raw.R run_sweep_BA_raw.R
run_sweep_real_raw.R

analysis/ extract_sweep_topology_tables.R plot\_\*

data/ enterotype_graph.graphml

outputs/ sweep_er/ sweep_ws/ sweep_ba/ sweep_real/

## Model overview

Each simulation proceeds as follows:

1.  A network is generated or loaded.
2.  A fraction of nodes is initialized as resistant.
3.  Resistance spreads through the network with probability p.
4.  At time AB an antibiotic perturbation removes susceptible nodes.
5.  The simulation continues until time T_used.

## Output

Each simulation produces:

-   g_pre: network snapshot before antibiotic exposure
-   g_post: network snapshot after the simulation
-   counts_ts: time series of resistant and susceptible nodes
-   metadata: parameters used in the simulation

Outputs are stored inside the outputs directory and organized by
parameter sweep.

## Requirements

R packages:

-   tidyverse
-   igraph
-   vroom

## Running simulations

Example:

Rscript scripts/run_sweep_ER_raw.R

Equivalent scripts exist for Watts--Strogatz, Barabási--Albert, and
empirical networks.

## License

MIT
