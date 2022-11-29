### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 6d467640-b19c-4f77-845d-f9b4aca62104
begin
	using Pkg
	push!(LOAD_PATH,joinpath("..")) # Stack with package environment
	Pkg.activate(joinpath("..","test"))
end

# ╔═╡ 7b4b635d-c600-4d6a-9cef-64bd97eb6a3e
begin
	using Test
	using Revise
	using Printf
	using VoronoiFVMDiffEq
	using DifferentialEquations
	using LinearAlgebra
	using PlutoUI
	using ExtendableGrids
	using DataStructures
	using GridVisualize,PlutoVista
	default_plotter!(PlutoVista);
end

# ╔═╡ adea3b1c-854f-4f9d-9013-ab6a2a6f7fd7
md"""
##  A Brusselator problem
"""

# ╔═╡ 9cefef6c-b487-4daa-a814-efad85f6a634
md"""
Two species diffusing and interacting via a reaction
```math
\begin{aligned}
  \partial_t u_1 - \nabla \cdot (D_1 \nabla u_1) &+ (B+1)u_1-A-u_1^2u_2  =0\\
  \partial_t u_2 - \nabla \cdot (D_2 \nabla u_2) &+ u_1^2u_2 -B u_1 =0\\
\end{aligned}
```
"""


# ╔═╡ 32107aac-050d-4f50-b95c-383e1bb38652
begin 
	const bruss_A=2.25
	const bruss_B=7.0
	const bruss_D_1=0.025
	const bruss_D_2=0.25
	const pert=0.1
	const bruss_tend=150
end;


# ╔═╡ 2fb1c53a-7ff8-4ac9-ae78-83bcbc57c926
function bruss_storage(f,u,node)
	f[1]=u[1]
	f[2]=u[2]
end;

# ╔═╡ 71b7e770-5cd4-4671-a76a-8e29eda04eec
function bruss_diffusion(f,u,edge)
	f[1]=bruss_D_1*(u[1,1]-u[1,2])
	f[2]=bruss_D_2*(u[2,1]-u[2,2])	
end;

# ╔═╡ f1e7a242-d631-4624-b0c0-ac44f139d77c
function bruss_reaction(f,u,node)
    f[1]= (bruss_B+1.0)*u[1]-bruss_A-u[1]^2*u[2]
    f[2]= u[1]^2*u[2]-bruss_B*u[1]
end;

# ╔═╡ 7e214c83-9c5c-40a9-8b00-db79dfec9a88
diffeqmethods=OrderedDict(
"Rosenbrock23 (Rosenbrock)" => Rosenbrock23,
"QNDF2 (Like matlab's ode15s)" =>  QNDF2,
"FBDF" => FBDF,
"Implicit Euler" => ImplicitEuler
)

# ╔═╡ 95be1da7-5f98-4a15-bd8e-7db1ee324768
begin
	
function ODESolver(system,inival,solver)
	problem = ODEProblem(system,inival,(0,bruss_tend))
	odesol = DifferentialEquations.solve(problem,
		                                 solver,
			                             dt=1.0e-5,reltol=1.0e-4)
    reshape(odesol,system)
end;

	sys0=VoronoiFVM.System(simplexgrid(0:0.1:1),species=[1,2],flux=bruss_diffusion, storage=bruss_storage, reaction=bruss_reaction);
	problem0 = ODEProblem(sys0,unknowns(sys0),(0,0.1))

	for method in diffeqmethods
		DifferentialEquations.solve(problem0,method.second())
   end
end

# ╔═╡ b35c0d6a-1b76-4396-992a-0e35d7d99cb1
html"""<hr>"""

# ╔═╡ dd91ee3d-0739-4ec0-a1cb-642be792594b
md"""
This notebook uses the test environment of the package and cannot be started outside of the `examples` subdirectory. To make it relocateable, make the next cell a markdown cell. Please be aware that this way, version information on dependencies gets lost.
"""

# ╔═╡ 2d439adf-c6b8-4e28-bb6a-58e00a9823ac
const Layout=PlutoUI.ExperimentalLayout;

# ╔═╡ 1462d783-93d3-4ad4-8701-90bde88c7553
Layout.hbox([
	md"""dim=$(@bind bruss_dim Scrubbable(1:2,default=1))""",
   md""" ``\quad``""",
	md""" method: $(@bind bruss_method Select([keys(diffeqmethods)...]))""",
    md""" ``\quad``""",
	md"""t=$(@bind t_bruss Slider(0:bruss_tend/1000:bruss_tend,show_value=true,default=bruss_tend))"""
])

# ╔═╡ d48ad585-9d0a-4b7e-a54b-3c76d8a5ca21
if bruss_dim==1
		bruss_X=-1:0.01:1
		bruss_grid=simplexgrid(bruss_X)
	else
		bruss_X=-1:0.1:1
		bruss_grid=simplexgrid(bruss_X,bruss_X)
end;

# ╔═╡ 719a15e1-7a69-4e70-b20e-d75fa448458e
bruss_system=VoronoiFVM.System(bruss_grid,species=[1,2],
			flux=bruss_diffusion, storage=bruss_storage, reaction=bruss_reaction);

# ╔═╡ 62a1dad1-b095-4df9-b1f8-e6a97084d8f8
begin
    inival=unknowns(bruss_system,inival=0)
	coord=bruss_grid[Coordinates]
	fpeak(x)=exp(-norm(10*x)^2)
	for i=1:size(inival,2)
   	 		inival[1,i]=fpeak(coord[:,i])
   	 		inival[2,i]=0
			#
    end
end

# ╔═╡ e71a2ed0-5f39-473f-87a0-6f61748f2793
t_run=@elapsed bruss_tsol=ODESolver(bruss_system,inival,diffeqmethods[bruss_method]());

# ╔═╡ c1da7d8e-2921-4366-91f0-dc8e1834595b
(t_run=t_run,details(bruss_system.history)...)

# ╔═╡ 4a4a2f78-4d4c-4b9b-883e-1e57de41c9a9
scalarplot(bruss_tsol.t[1:end-1],bruss_tsol.t[2:end]-bruss_tsol.t[1:end-1],yscale=:log,resolution=(500,200),xlabel="t",ylabel="Δt",title="timesteps")

# ╔═╡ e7a8aae1-8e7a-4b7d-8ce6-701ea586b89a
bruvis=(GridVisualizer(;size=(300,200),dim=bruss_dim),GridVisualizer(;size=(300,200),dim=bruss_dim));bruvis

# ╔═╡ 59049ff3-8f18-4adc-aded-2ffb752ce574
bruss_sol=bruss_tsol(t_bruss);

# ╔═╡ 52b2c6dc-650d-4a1f-bb94-757c9cff9ee6
scalarplot!(bruvis[1],bruss_grid,bruss_sol[1,:],limits=(0,10),show=true,colormap=:summer)

# ╔═╡ f6571e3b-ef2a-415a-a4e0-ee88d84e505e
scalarplot!(bruvis[2],bruss_grid,bruss_sol[2,:],limits=(0.5,3),show=true,colormap=:summer)

# ╔═╡ Cell order:
# ╠═7b4b635d-c600-4d6a-9cef-64bd97eb6a3e
# ╟─adea3b1c-854f-4f9d-9013-ab6a2a6f7fd7
# ╟─9cefef6c-b487-4daa-a814-efad85f6a634
# ╠═32107aac-050d-4f50-b95c-383e1bb38652
# ╠═2fb1c53a-7ff8-4ac9-ae78-83bcbc57c926
# ╠═71b7e770-5cd4-4671-a76a-8e29eda04eec
# ╠═f1e7a242-d631-4624-b0c0-ac44f139d77c
# ╠═95be1da7-5f98-4a15-bd8e-7db1ee324768
# ╠═7e214c83-9c5c-40a9-8b00-db79dfec9a88
# ╠═d48ad585-9d0a-4b7e-a54b-3c76d8a5ca21
# ╠═719a15e1-7a69-4e70-b20e-d75fa448458e
# ╠═62a1dad1-b095-4df9-b1f8-e6a97084d8f8
# ╠═e71a2ed0-5f39-473f-87a0-6f61748f2793
# ╟─c1da7d8e-2921-4366-91f0-dc8e1834595b
# ╟─1462d783-93d3-4ad4-8701-90bde88c7553
# ╟─e7a8aae1-8e7a-4b7d-8ce6-701ea586b89a
# ╠═4a4a2f78-4d4c-4b9b-883e-1e57de41c9a9
# ╠═59049ff3-8f18-4adc-aded-2ffb752ce574
# ╠═52b2c6dc-650d-4a1f-bb94-757c9cff9ee6
# ╠═f6571e3b-ef2a-415a-a4e0-ee88d84e505e
# ╟─b35c0d6a-1b76-4396-992a-0e35d7d99cb1
# ╟─dd91ee3d-0739-4ec0-a1cb-642be792594b
# ╠═6d467640-b19c-4f77-845d-f9b4aca62104
# ╟─2d439adf-c6b8-4e28-bb6a-58e00a9823ac
