[![linux-macos-windows](https://github.com/j-fu/VoronoiFVMDiffEq.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/j-fu/VoronoiFVMDiffEq.jl/actions/workflows/ci.yml)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://j-fu.github.io/VoronoiFVMDiffEq.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://j-fu.github.io/VoronoiFVMDiffEq.jl/dev)


VoronoiFVMDiffEq.jl
===================

Glue package between [VoronoiFVM.jl](https://github.com/j-fu/VoronoiFVM.jl) and [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)

The package extends the constructors for `ODEFunction` and `ODEProblem` by methods taking in a `VoronoiFVM.System`:
```
using VoronoiFVMDiffEq, DifferentialEquations
system = VoronoiFVM.System(...)
inival = unknowns(system)
problem = ODEProblem(system,inival,tspan)
odesolution = DifferentialEquations.solve(problem, QNDF2())
voronoifvmsolution = reshape(odesolution, system)
```
Instead of [`QNDF2`](https://sciml.ai/news/2021/05/24/QNDF/) you can try all mass matrix form capable stiff ode solvers form the  [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) suite.

The package re-exports all of VoronoiFVM, so that it is sufficient to `use` this package instead of [VoronoiFVM.jl](https://github.com/j-fu/VoronoiFVM.jl).

The package requires a recent Julia version (currently 1.8) due to significant advances in package precompilation.
