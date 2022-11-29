VoronoiFVMDiffEq.jl
==================

Glue package between [VoronoiFVM.jl](https://github.com/j-fu/VoronoiFVM.jl) and [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)

The package extends the constructors for `ODEFunction` and `ODEProblem` by methods taking in a `VoronoiFVM.System`:

```
problem=ODEProblem(sys::VoronoiFVM.System,inival,tspan,callback=DifferentialEquations.CallbackSet())
odesolution=DifferentialEquations.solve(problem)
voronoifvmsolution=reshape(odesolution, sys::VoronoiFVM.System)
```

It re-exports VoronoiFVM, so that it is sufficient to `use` this package instead of VoronoiFVM.


