[![linux-macos-windows](https://github.com/j-fu/VoronoiFVMDiffEq.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/j-fu/VoronoiFVMDiffEq.jl/actions/workflows/ci.yml)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://j-fu.github.io/VoronoiFVMDiffEq.jl/dev)


VoronoiFVMDiffEq.jl
===================

Glue package between [VoronoiFVM.jl](https://github.com/j-fu/VoronoiFVM.jl) and [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)

The package extends the constructors for `ODEFunction` and `ODEProblem` by methods taking in a `VoronoiFVM.System`:

```
problem=ODEProblem(sys::VoronoiFVM.System,inival,tspan,callback=DifferentialEquations.CallbackSet())
odesolution=DifferentialEquations.solve(problem)
voronoifvmsolution=reshape(odesolution, sys::VoronoiFVM.System)
```

It re-exports VoronoiFVM, so that it is sufficient to `use` this package instead of VoronoiFVM.

The package requires a recent Julia version (currently 1.8) due to significant advances in package precompilation.


This package is in pre-release state. Versions 0.0.x are registered in https://github.com/j-fu/PackageNursery.
If you trust this registry, you can issue once

```
pkg> registry add https://github.com/j-fu/PackageNursery
```
and add the package via  
```
pkg> add VoronoiFVMDiffEq
```
