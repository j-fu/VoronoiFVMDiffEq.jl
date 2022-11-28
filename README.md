VoronoiFVMDiffEq.jl
==================

Glue package between [VoronoiFVM.jl](https://github.com/j-fu/VoronoiFVM.jl) and [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)

The package extends the constructors for `ODEFunction` and `ODEProblem` by methods taking in a `VoronoiFVM.System`:

```
ODEFunction(sys::VoronoiFVM.System; jacval=unknowns(sys,0), tjac=0)
ODEProblem(sys::VoronoiFVM.System,inival,tspan,callback=DifferentialEquations.CallbackSet())
reshape(odesolution, sys::VoronoiFVM.System)
```


