module VoronoiFVMDiffEq


import DifferentialEquations
using RecursiveArrayTools
using Reexport

@reexport using VoronoiFVM


"""
     ODEFunction(system,inival=unknowns(system,inival=0),t0=0)
    
Create an [ODEPFunction](https://diffeq.sciml.ai/stable/basics/overview/#Defining-Problems) which can 
be handeled by ODE solvers from [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).

Parameters:
- `system`: A [`VoronoiFVM.System`](@ref)
- `jacval` (optional): Initial value. Consider to  pass a stationary solution at time `tjac`.
- `tjac` (optional): tjac

The `jacval` and `tjac` are passed  for a first evaluation of the Jacobian, allowing to detect
the sparsity pattern which is passed to the DifferentialEquations.jl solver.
"""
function DifferentialEquations.ODEFunction(sys::VoronoiFVM.AbstractSystem; jacval=unknowns(sys,0), tjac=0)
    DifferentialEquations.ODEFunction(eval_rhs!;
                                      jac=eval_jacobian!,
                                      jac_prototype=prepare_diffeq!(sys,jacval, tjac),
                                      mass_matrix=mass_matrix(sys))
end


"""
    ODEProblem(system,inival,tspan,callback=DifferentialEquations.CallbackSet())
    
Create an [ODEProblem](https://diffeq.sciml.ai/stable/basics/overview/#Defining-Problems) which can 
be handeled by ODE solvers from [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).

Parameters:
- `system`: A [`VoronoiFVM.System`](@ref)
- `inival`: Initial value. Consider to  pass a stationary solution at `tspan[1]<`.
- `tspan`: Time interval 
- `callback` : (optional) [callback](https://diffeq.sciml.ai/stable/features/callback_functions/#Using-Callbacks) for ODE solver 

The method returns an [ODEProblem](https://diffeq.sciml.ai/stable/basics/overview/#Defining-Problems) which can be solved
by [DifferentialEquations.solve()](https://diffeq.sciml.ai/stable/basics/common_solver_opts/).
"""
function DifferentialEquations.ODEProblem(sys::VoronoiFVM.AbstractSystem, inival, tspan, callback=DifferentialEquations.CallbackSet())
    odefunction=DifferentialEquations.ODEFunction(sys; jacval=inival,tjac=tspan[1])
    DifferentialEquations.ODEProblem(odefunction,vec(inival),tspan,sys,callback)
end

"""
    reshape(ode_solution, system)
Create a [`TransientSolution`](@ref) from the output of the ode solver.
"""
function Base.reshape(sol::AbstractDiffEqArray, sys::VoronoiFVM.AbstractSystem)
     TransientSolution([reshape(sol.u[i],sys) for i=1:length(sol.u)] ,sol.t)
end


end
