### A Pluto.jl notebook ###
# v0.19.22

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
	Pkg.activate(joinpath(@__DIR__,"..","docs"))# activate test environment
    Pkg.instantiate()
end

# ╔═╡ 7b4b635d-c600-4d6a-9cef-64bd97eb6a3e
begin
	using Test	
	using Revise
	using Printf
	using VoronoiFVMDiffEq
	using DifferentialEquations
	using DASSL
	using LinearAlgebra
	using PlutoUI, HypertextLiteral,UUIDs
	using DataStructures
	using GridVisualize,PlutoVista
end

# ╔═╡ 57deb502-bac9-4093-86e3-e8d25eb02df2
md"""
## 1D Nonlinear Storage
"""

# ╔═╡ a2dd09f8-ffe8-4b34-b2bf-21266ad97a76
md"""
This equation comes from the transformation of the nonlinear diffusion equation
```math
\partial_t v - \Delta v^m = 0
```
to 
```math
\partial_t u^\frac{1}{m} -\Delta u = 0
```
in $\Omega=(-1,1)$ with homogeneous Neumann boundary conditions.
We can derive an exact solution from the Barenblatt solution of the equation for u.
"""

# ╔═╡ bf88f84f-8633-4ed3-aae6-5b937f8ea470
function barenblatt(x,t,m)
    tx=t^(-1.0/(m+1.0))
    xx=x*tx
    xx=xx*xx
    xx=1- xx*(m-1)/(2.0*m*(m+1));
    if xx<0.0
        xx=0.0
    end
    return tx*xx^(1.0/(m-1.0))
end

# ╔═╡ 46a0da67-0418-4ac2-abc5-9dba3abf306c
begin
	const m=2
	const ε=1.0e-10
	const n=50
	const t0=1.0e-3
	const tend=1.0e-2
end

# ╔═╡ f4eb15d1-4e0a-4738-bbc5-d5a7ced80289
X=collect(-1:2.0/n:1)

# ╔═╡ c4f6f616-996b-4984-82c3-fd3f1301fb7c
u0=map(x->barenblatt(x,t0,m)^m,X)

# ╔═╡ 27bd76a1-1739-4faf-9fd0-e095ee3f43f3
begin
    grid=VoronoiFVM.Grid(X)
end

# ╔═╡ 108cf7c1-17c7-4f56-b1f7-f38dfe8c857e
md"""
### Direct implementation with VoronoiFVM
"""

# ╔═╡ 54aef797-210b-4eea-95b4-e0bdcb478c1d
    function flux!(f,u,edge)
        f[1]=u[1,1]-u[1,2]
    end

# ╔═╡ c33cfa04-7418-4b1b-868d-5e69b16ea2d9
md"""
Storage term needs to be regularized as its derivativ at 0 is infinity:
"""

# ╔═╡ 3790d063-8a3d-429b-8a59-05cfc35b6878
    function storage!(f,u,node)
        f[1]=(ε+u[1])^(1.0/m)
    end

# ╔═╡ 62f5e4f1-3635-4de8-9052-1771fa3a31cf
begin
    physics=VoronoiFVM.Physics(
        flux=flux!,
        storage=storage!)
    sys=VoronoiFVM.System(grid,physics,species=1)
	    inival=unknowns(sys)
	inival[1,:]=u0
	control=VoronoiFVM.SolverControl()
    tsol=VoronoiFVM.solve(sys;inival,times=(t0,tend),Δt_min=1.0e-4,Δt=1.0e-4,Δu_opt=0.1,force_first_step=true,log=true)
	summary(sys.history)
end

# ╔═╡ 101754f3-8e05-4cf0-8770-e17efbbf4f82
md"""
### Implementation as DAE
"""

# ╔═╡ 1d3e90c6-ba1a-4d3c-bc09-dc78e8aeafe6
md"""
If we want to solve the problem with DifferentialEquations.jl, we see that the problem structure does not fit into the setting of that package due to the nonlinearity under the time derivative. Here we propose a reformulation to a DAE as a way to achieve this possibility:

```math
\begin{cases}
	\partial_t w -\Delta u &= 0\\
               w^m - u &=0
\end{cases}
```
"""

# ╔═╡ d495a088-69dd-4f4e-964c-88638958799a
function dae_storage!(y,u,node)
	y[1]=u[2]
end

# ╔═╡ 4c09ec49-c4d6-4c5a-8f73-8a736acaff61
function dae_reaction!(y,u,node)
	y[2]= u[2]^m-u[1]
end

# ╔═╡ 39e4fd5f-70a6-41e9-a02d-b11709762c19
md"""
First, we test this with the implicit Euler method of VoronoiFVM
"""

# ╔═╡ 89212312-7846-4a2f-a4b8-bd7b61235cf9
begin
    dae_physics=VoronoiFVM.Physics(
        flux=flux!,
       storage=dae_storage!,
		reaction=dae_reaction!
	)
    dae_sys=VoronoiFVM.System(grid,dae_physics,species=[1,2])
	dae_inival=unknowns(dae_sys)
	dae_inival[1,:].=u0
	dae_inival[2,:].=u0.^(1/m)
	dae_control=VoronoiFVM.SolverControl()
    dae_tsol=VoronoiFVM.solve(dae_sys;inival=dae_inival,times=(t0,tend),Δt_min=1.0e-4,Δt=1.0e-4,Δu_opt=0.1,force_first_step=true,log=true)
	summary(dae_sys.history)
end

# ╔═╡ b4f4baa2-64a2-464c-9e19-a39b490d210a
md"""
### Implementation via DifferentialEquations.jl
"""

# ╔═╡ 8079ba59-7595-4109-805a-32d135d383f9
diffeqmethods=OrderedDict(
"QNDF2 (Like matlab's ode15s)" =>  QNDF2,
"Rodas5" => Rodas5,
"Rosenbrock23 (Rosenbrock)" => Rosenbrock23,
"FBDF" => FBDF,
"Implicit Euler" => ImplicitEuler
)


# ╔═╡ 5e09155b-1b9e-43a2-983f-4be15557b513
md"""
method: $(@bind method Select([keys(diffeqmethods)...]))
"""


# ╔═╡ da7645c2-d254-4886-b2b6-28289368fc22
begin
	de_sys=VoronoiFVM.System(grid,dae_physics,species=[1,2])
	problem = ODEProblem(de_sys,dae_inival,(t0,tend));
    de_odesol=DifferentialEquations.solve(problem,
		diffeqmethods[method](linsolve=UMFPACKFactorization(;reuse_symbolic=false)),
		adaptive=true,
        reltol=1.0e-3,
		abstol=1.0e-3,
		initializealg=NoInit()
		)          
		de_tsol=reshape(de_odesol,de_sys)
end;

# ╔═╡ b35c0d6a-1b76-4396-992a-0e35d7d99cb1
html"""<hr>"""

# ╔═╡ 0488cf35-1e59-4761-99ed-e91c75259403
function myaside(x;top=1)
	uuid=uuid1()
	@htl("""
		<style>
		
		
		@media (min-width: calc(700px + 30px + 300px)) {
			aside.plutoui-aside-wrapper-$(uuid) {

	color: var(--pluto-output-color);
	position:fixed;
	right: 1rem;
	top: $(top)px;
	width: 400px;
	padding: 10px;
	border: 3px solid rgba(0, 0, 0, 0.15);
	border-radius: 10px;
	box-shadow: 0 0 11px 0px #00000010;
	/* That is, viewport minus top minus Live Docs */
	max-height: calc(100vh - 5rem - 56px);
	overflow: auto;
	z-index: 40;
	background-color: var(--main-bg-color);
	transition: transform 300ms cubic-bezier(0.18, 0.89, 0.45, 1.12);
	
			}
			aside.plutoui-aside-wrapper > div {
#				width: 300px;
			}
		}
		</style>
		
		<aside class="plutoui-aside-wrapper-$(uuid)">
		<div>
		$(x)
		</div>
		</aside>
		
		""")
end


# ╔═╡ 9340c3a2-12f9-4f0e-9e5b-3c960388f9cc
myaside(md"""$(vis=GridVisualizer(resolution=(380,200),dim=1,Plotter=PlutoVista,legend=:lt);)""",top=300)

# ╔═╡ 1a5c5b50-8aea-47b9-9961-5d74e93e9d69
myaside(md"""
t=$(@bind t  PlutoUI.Slider(range(t0,tend,length=10001),default=2*t0;show_value=true)) 
""",top=550)

# ╔═╡ 7fdef034-feac-49e4-9b98-bf579ac5fd94
let
	u=tsol(t)
	u_dae=dae_tsol(t)
	u_de=de_tsol(t)
	scalarplot!(vis,X,map(x->barenblatt(x,t,m).^m,X),clear=true,color=:red,linestyle=:solid,flimits=(0,100),label="exact")
	scalarplot!(vis,grid,u_dae[1,:],clear=false,color=:green, linestyle=:solid,label="vfvm_dae")
	scalarplot!(vis,grid,u_de[1,:],clear=false,color=:blue, markershape=:cross,linestyle=:dot,label="vfvm_diffeq")
	scalarplot!(vis,grid,u[1,:],clear=false,color=:black,markershape=:none, linestyle=:dash,title="t=$(t)",label="vfvm_default")
	reveal(vis)
end

# ╔═╡ 785999c3-62d9-49ee-a890-70ec745211c1
TableOfContents()

# ╔═╡ dd91ee3d-0739-4ec0-a1cb-642be792594b
md"""
This notebook uses the documentation environment of the package and cannot be started outside of the `examples` subdirectory. To make it relocateable, make the next cell a markdown cell and restart the notebook. Please be aware that this way, version information is rebuilt from scratch.
"""

# ╔═╡ Cell order:
# ╠═7b4b635d-c600-4d6a-9cef-64bd97eb6a3e
# ╟─57deb502-bac9-4093-86e3-e8d25eb02df2
# ╟─a2dd09f8-ffe8-4b34-b2bf-21266ad97a76
# ╠═bf88f84f-8633-4ed3-aae6-5b937f8ea470
# ╠═46a0da67-0418-4ac2-abc5-9dba3abf306c
# ╠═f4eb15d1-4e0a-4738-bbc5-d5a7ced80289
# ╠═c4f6f616-996b-4984-82c3-fd3f1301fb7c
# ╠═27bd76a1-1739-4faf-9fd0-e095ee3f43f3
# ╟─108cf7c1-17c7-4f56-b1f7-f38dfe8c857e
# ╠═54aef797-210b-4eea-95b4-e0bdcb478c1d
# ╟─c33cfa04-7418-4b1b-868d-5e69b16ea2d9
# ╠═3790d063-8a3d-429b-8a59-05cfc35b6878
# ╠═62f5e4f1-3635-4de8-9052-1771fa3a31cf
# ╠═9340c3a2-12f9-4f0e-9e5b-3c960388f9cc
# ╠═1a5c5b50-8aea-47b9-9961-5d74e93e9d69
# ╠═7fdef034-feac-49e4-9b98-bf579ac5fd94
# ╟─101754f3-8e05-4cf0-8770-e17efbbf4f82
# ╟─1d3e90c6-ba1a-4d3c-bc09-dc78e8aeafe6
# ╠═d495a088-69dd-4f4e-964c-88638958799a
# ╠═4c09ec49-c4d6-4c5a-8f73-8a736acaff61
# ╟─39e4fd5f-70a6-41e9-a02d-b11709762c19
# ╠═89212312-7846-4a2f-a4b8-bd7b61235cf9
# ╟─b4f4baa2-64a2-464c-9e19-a39b490d210a
# ╟─8079ba59-7595-4109-805a-32d135d383f9
# ╟─5e09155b-1b9e-43a2-983f-4be15557b513
# ╠═da7645c2-d254-4886-b2b6-28289368fc22
# ╟─b35c0d6a-1b76-4396-992a-0e35d7d99cb1
# ╟─0488cf35-1e59-4761-99ed-e91c75259403
# ╟─785999c3-62d9-49ee-a890-70ec745211c1
# ╟─dd91ee3d-0739-4ec0-a1cb-642be792594b
# ╠═6d467640-b19c-4f77-845d-f9b4aca62104
