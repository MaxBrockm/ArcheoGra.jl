# ArcheoGra.jl
Welcome to ArcheoGra.jl, a Julia package designed to reconstruct possible prehistoric human migration routes and social structures from archaeological data and calculate Least-Cost paths using adaptive Wavelet methods.

## Overview
For creating archaeological networks, we use archaeological sites as the baseline of confirmed human presence. From there, migration routes are modeled as edges in a spatial network, where the archaeological sites act as the vertices of the network. Vertices are connected probabilistically.

The second part, concerns the calculation of Least-Cost paths on topographic data. In order to improve the runtime of the Least-Cost path analysis, topographic data is represented adaptively. LCPs are calculated via the adpative Dijkstra's Algorithm.

## Modules and Key Features
The package is divided into two main submodules:
* Connect Network (ConNetw): Provides the foundational tools to construct networks from archaeological datasets, handling initial setup, global parameters, and fundamental data preparation required to model prehistoric mobility and social connectivity.
* Least Cost Path Analysis (LCPA): Calculates realistic travel distances between archaeological sites using topographical data. Instead of relying on straight-line distances, it models movement based on energy expenditure, taking elevation changes into account. This submodule employs adaptive Wavelet methods to find LCPs using an adaptive Dijkstra's Algorithm.

## Installation
You can install the package using Julia's package manager. From the Julia REPL, type ] to enter the Pkg prompt and run:

## Citation
We encourage you to cite our work if you have used it in your research. The latest citation information can always be found in the `CITATION.bib` file within the repository.
Key references include:
* Brockmann, Perlberg, Kunoth, "Efficient Path Reconstruction in Prehistoric Human Migration: An Adaptive Dijkstra's Algorithm Based on Wavelet Compression for Topographic Data", 2026, (preprint).
* Hilpert, Brockmann, Maier, et.al. An interdisciplinary approach to mobility network building in Archaeology (in preparation).
