# Shared funtions across all modules
These core functions are utilized globally across both the `ConNetw` and `LCPA` submodules to ensure consistent spatial calculations and reproducible stochastic processes.

```@docs
earth_radius
set_seed
distance(u::AbstractArray{Float64},v::AbstractArray{Float64})
```