````@eval
using Markdown
Markdown.parse("""
$(read("../../README.md",String))
""")
````

### Choice of ODE/DAE solver

As this package interfaces to  the PDE solver package [VoronoiFVM.jl](https://github.com/j-fu/VoronoiFVM.jl),
a general advise is to choose methods able to handle [stiff problems](https://diffeq.sciml.ai/stable/solvers/ode_solve/#Stiff-Problems).
Moreover, often, discretized PDE systems (e.g. containing elliptic equations) are differential agebraic equation (DAE) systems 
which should be solved by [DAE solvers](https://diffeq.sciml.ai/stable/solvers/dae_solve/).


### Methods
```@docs
SciMLBase.ODEProblem
SciMLBase.ODEFunction
Base.reshape
```

