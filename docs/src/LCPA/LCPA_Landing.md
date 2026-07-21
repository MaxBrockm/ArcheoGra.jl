# Least Cost Path Analysis (LCPA)

This submodule provides tools to calculate realistic travel distances between archaeological sites using topographical data. Instead of relying on straight-line distances (Haversine), the network models human movement based on energy expenditure, taking elevation changes into account. For references on LCPA, Wavelet methods, and adaptive Dijsktra's Algorithm see [Brockmann, Perlberg, Kunoth (2026)](https://arxiv.org/abs/2607.12408) and the references therein.

## Core Data Structures

To efficiently handle the spatial and topographical data, the following types are defined:

```@docs
TopographicalData
ConnectionFunction
```