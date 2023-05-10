### A Pluto.jl notebook ###
# v0.19.25

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

# ╔═╡ e9f8097f-8dcc-4924-9eee-d8b0f49b0e4d
begin
    using Pkg
    Pkg.activate(joinpath(@__DIR__,"..","docs"))# activate test environment
    Pkg.instantiate()
end


# ╔═╡ 5e728bfa-7571-11ed-09a3-e3fe7605ab0d
begin
	using Test
	using Revise
	using Printf
	using VoronoiFVMDiffEq
	using DifferentialEquations
	using ExtendableGrids
	using LinearAlgebra
	using PlutoUI, HypertextLiteral,UUIDs
	using GridVisualize,PlutoVista
	using DataStructures
end

# ╔═╡ fcdf6695-c96c-4579-ae57-53ec332a60ae
md"""

Consider the standard van Roosbroeck system in a n-doped region

```math
\begin{aligned}
- \nabla \cdot (\varepsilon_s \nabla \psi) &= (z_n n_n + z_p n_p - z_n C_n)\\
     z_n \partial_t n_n + \nabla \cdot \mathbf{j}_n &= 0,\\
z_p \partial_t n_p + \nabla \cdot \mathbf{j}_p &= 0
\end{aligned}
```
with the boundary conditions

```math
\begin{aligned}
\psi &= \psi_s + \Delta u\\
    \varphi_n &= \Delta u,\\
	\varphi_p &= \Delta u.
\end{aligned}
```
As ``\mathbf{j}_n, \mathbf{j}_p `` we consider nonlinear drift-diffusion fluxes.

The system can be rewritten in terms of quasi Fermi potentials via

``\begin{aligned}
    n_\alpha = \mathcal{F}_\alpha(z_\alpha (\varphi_\alpha - \psi)), \alpha = n,p,
\end{aligned}
``
where ``\mathcal{F}_\alpha `` is the statistics function.
"""

# ╔═╡ 2d35f805-34f1-4f99-b226-dad8924f5bd3
begin
	const h            = 5.0

	const zn           = - 1;  const zp           =   1
	const εr           = 1.0e-2
	const barrierLeft  = 1.0e-2;  const barrierRight = 1.0e-2

	const Nd           = 1.0;
end

# ╔═╡ 88d0c7bd-6d8a-47e1-8064-3d1de6f46a3f
begin
	const iphin   = 1
	const iphip   = 2
	const ipsi    = 3
	const inn     = 4
	const inp     = 5
	const species = [iphin, iphip, ipsi, inn, inp]
end

# ╔═╡ 73339ba6-71fc-4c25-9b9f-50a780d84ac2
begin
	## time mesh
	const endTime       = 9.0
	const amplitude     = 15.0
	const scanrate      = 3*amplitude/endTime
	const endTimeSim    = 0.5
end

# ╔═╡ c2e90dcc-35a3-441f-8cee-2a1a9fb1e391
begin
	function statistics(x::Real)
		if x < 1.6107
			z=log1p(exp(x))
			return ( 1 + 0.16 * z ) * z
		elseif 1.6107 <= x <= 344.7
			z = log1p(exp( x^(3/4)) )
			return 0.3258 - (0.0321  - 0.7523 * z ) * z
		else
			z = x^(3/4)
			return 0.3258 - (0.0321 -  0.7523 * z ) * z
		end
	end

	function biasValue(t)

		biasVal = 0.0
		if      t <= 3.0
			biasVal = 0.0 + scanrate * t
		elseif  t >= 3.0  && t <= 6.0
			biasVal = 3 * amplitude .- 2 * scanrate * t
		elseif  t >= 6.0 && t <= 9.0
			biasVal = - 3 * amplitude .+ scanrate * t
		end

		return biasVal

	end

end

# ╔═╡ 80bf2ced-c5f9-4ae3-95f2-e450b85055a2
mutable struct Data
	inEquilibrium::Bool
	Data() = new()
end

# ╔═╡ a5bf9ba2-2e59-4616-8c17-d343bbc7ef71
begin
	data               = Data()
	data.inEquilibrium = true;
end

# ╔═╡ 05d39cf1-85d8-41de-b50d-94ec00f21cc2
begin
	h_total  = h

	# region numbers
	region   = 1
	regions  = [region]

	# boundary region numbers
	bregion1 = 1
	bregion2 = 2
	bregions = [bregion1, bregion2]

	coord1   = geomspace(0.0, h/2, 7e-4 * h, 1e-2 * h)
	coord2   = geomspace(h/2, h,   1e-2 * h, 7e-4 * h)
	coord    = glue(coord1, coord2)

	grid     = simplexgrid(coord)

	# specify inner regions
	cellmask!(grid, [0.0], [h], region)

	# specifiy outer regions
	bfacemask!(grid, [0.0],     [0.0],     bregion1)
	bfacemask!(grid, [h_total], [h_total], bregion2)
	gridplot(grid,legend=:lt, fignumber = 1, resolution=(600,200))
end

# ╔═╡ 5481508b-597a-4ff8-bce4-30612b2d8226
md"""
# Current implementation with VoronoiFVM (with potentials)
"""

# ╔═╡ 49ef7e60-a27b-49fa-891e-bfad5c9db78c
md"""
In the current implementation we use the quasi Fermi potentials with the electric potential as set of unknowns.
"""

# ╔═╡ d9f62141-5856-4ecd-b1c5-b6f44597e48e
md"""
### Physics environment
"""

# ╔═╡ 9da3f051-8ab0-44ee-9c17-1453935487f1
begin

	function reactionqF!(f, u, node, data)
		if data.inEquilibrium
			f[iphin] = u[iphin]
			f[iphip] = u[iphip]
		else
			f[iphin] = 0.0
			f[iphip] = 0.0
		end

		########################
		etan    = zn * ( (u[iphin] - u[ipsi]) )
		etap    = zp * ( (u[iphip] - u[ipsi]) )

		nn      = statistics(etan)
		np      = statistics(etap)

		f[ipsi] = - (zn * (nn - Nd) + zp * np)
	end


	function fluxqF!(f, u, edge, data)
		dpsi    = u[ipsi, 2] - u[ipsi, 1]
		f[ipsi] =  - εr * dpsi

		####################################
		etan1   = zn * ( (u[iphin, 1] - u[ipsi, 1]))
        etan2   = zn * ( (u[iphin, 2] - u[ipsi, 2]))

		nn1     = statistics(etan1)
		nn2     = statistics(etan2)
        if statistics(etan1).value ≈ 0.0 || statistics(etan2).value ≈ 0.0
            bp, bm  = fbernoulli_pm(zn * dpsi)
        else
            bp, bm  = fbernoulli_pm(zn * dpsi + (etan2 - etan1) -log(statistics(etan2)) + log(statistics(etan1)) )
        end
		f[iphin]  = - zn * ( bm * nn2 - bp * nn1)
		####################################
		etap1   = zp * ( (u[iphip, 1] - u[ipsi, 1]))
        etap2   = zp * ( (u[iphip, 2] - u[ipsi, 2]))

		np1     = statistics(etap1)
		np2     = statistics(etap2)
        if statistics(etap1).value ≈ 0.0 || statistics(etap2).value ≈ 0.0
            bp, bm  = fbernoulli_pm(zp * dpsi)
        else
            bp, bm  = fbernoulli_pm(zp * dpsi + (etap2 - etap1) -log(statistics(etap2)) + log(statistics(etap1)) )
        end
		f[iphip]  = - zp * ( bm * np2 - bp * np1)
	end

	function bconditionqF!(y, u, bnode, data)
		if bnode.region == 1
			barrier = barrierLeft
			Δu      = 0.0
		elseif bnode.region == 2
			barrier = barrierRight

		end

		if bnode.time == Inf
			Δu = 0.0
		else
			Δu = biasValue(bnode.time)
		end
		############################################
		boundary_dirichlet!(y, u, bnode, species=ipsi, region=1, value=(- (barrierLeft  )) )
		boundary_dirichlet!(y, u, bnode, species=ipsi, region=2, value=(- (barrierRight ) + Δu) )
		############################################

		boundary_dirichlet!(y, u, bnode, species=iphin, region=1, value = 0.0 )
		boundary_dirichlet!(y, u, bnode, species=iphip, region=1, value = 0.0 )

		boundary_dirichlet!(y, u, bnode, species=iphin, region=2, value = Δu)
		boundary_dirichlet!(y, u, bnode, species=iphip, region=2, value = Δu)
	end

	function storageqF!(f, u, node, data)

		etan     = zn * ( (u[iphin] - u[ipsi]) )
		etap     = zp * ( (u[iphip] - u[ipsi]) )

        nn       = statistics(etan)
        np       = statistics(etap)

        f[iphin] = zn * nn
        f[iphip] = zp * np
	end
end

# ╔═╡ 4e2b4ba4-ac05-4098-9f08-fd28e1cf9cb8
begin
	sysqF = VoronoiFVM.System(grid, data            = data,
									flux            = fluxqF!,
									reaction        = reactionqF!,
									bcondition      = bconditionqF!,
									storage         = storageqF!,
									unknown_storage =:sparse)
	enable_species!(sysqF, iphin, regions)
    enable_species!(sysqF, iphip, regions)
    enable_species!(sysqF, ipsi,  regions)
end

# ╔═╡ 36cc28cd-1737-41bd-8127-c784eb532443
md"""
## Calculate initial condition
"""

# ╔═╡ 73917cbf-113b-4947-91d6-6cb0aac2dafc
begin
	data.inEquilibrium = true
	solEQqF=VoronoiFVM.solve(sysqF)
end

# ╔═╡ f46d0e11-e30c-420d-87ed-242235a2fdd9
begin

	nnEQqF = statistics.(zn .* (solEQqF[iphin, :] - solEQqF[ipsi, :]))
	npEQqF = statistics.(zp .* (solEQqF[iphip, :] - solEQqF[ipsi, :]))
end

# ╔═╡ ee968c23-170b-4c66-b258-fbf4eebde2ce
md"""
let
	PyPlot.clf()

	subplot(211)
	PyPlot.plot(coord, solEQqF[ipsi, :], linewidth = 4, color = "royalblue", label ="\$\\psi\$")
	PyPlot.plot(coord, solEQqF[iphin, :], linewidth = 4, color = "green", label ="\$\\varphi_n\$")
	PyPlot.plot(coord, solEQqF[iphip, :], linewidth = 4, color = "red", label ="\$\\varphi_p\$", linestyle = ":")

	PyPlot.xlabel("space [m]", fontsize=18)
	PyPlot.ylabel("potential [V]", fontsize=18)
	PyPlot.tick_params(axis="both", labelsize=18)
	PyPlot.title("Equilibrium")
	PyPlot.legend(fancybox = true, loc = "lower center", fontsize=15)
	PyPlot.grid()
	################################################
	################################################
	subplot(212)

	PyPlot.semilogy(coord, nnEQqF, linewidth = 4, color="green", label ="\$n_{n}\$ ")
	PyPlot.semilogy(coord, npEQqF, linewidth = 4, color="red", label ="\$n_{p}\$ ")

	PyPlot.xlabel("space [m]", fontsize=18)
	PyPlot.ylabel("density [m\$^{-3}\$]", fontsize=18)
	PyPlot.tick_params(axis="both", labelsize=18)
	PyPlot.legend(fancybox = true, loc = "lower center", fontsize=15)
	PyPlot.grid()
	PyPlot.tight_layout()
	figure(1)
end
"""

# ╔═╡ 39238121-20f0-49d2-a6e9-96236d4a414e
let
	vis=GridVisualizer(Plotter=PlutoVista,size=(600,300),
	xlabel="x/m",legend=:cc)

scalarplot!(vis,coord, solEQqF[ipsi, :], linewidth = 2, color = "royalblue", label ="\$\\psi\$")
scalarplot!(vis,coord, solEQqF[iphin, :], linewidth = 2, color = "green", label ="\$\\varphi_n\$",clear=false)
scalarplot!(vis,coord, solEQqF[iphip, :], linewidth = 2, color = "red", label ="\$\\varphi_p\$", linestyle = ":",clear=false)

	reveal(vis)
end

# ╔═╡ a306246a-ae46-4b8b-8539-df9a1e279f4f
let
vis=GridVisualizer(Plotter=PlutoVista,size=(600,300),
	xlabel="x/m",legend=:cc)

scalarplot!(vis,coord, nnEQqF, linewidth = 2, color="green", label ="\$n_{n}\$ ")
scalarplot!(vis,coord, npEQqF, linewidth = 2, color="red", label ="\$n_{p}\$ ",clear=false)
reveal(vis)
end

# ╔═╡ 612906c4-13bd-4534-a184-ee892dcd6d03
md"""
## Time solve
"""

# ╔═╡ 2d17f616-5dc6-4bfa-952d-c21f06e3ed6f
begin
	data.inEquilibrium = false
	soltqF = VoronoiFVMDiffEq.solve(sysqF, inival = solEQqF, times=(0.0, endTimeSim))
	tvalqF = soltqF.t
end

# ╔═╡ 485bb7b8-9b32-480a-a4d0-aab0426ece9b
md"""
# New implementation with VoronoiFVM (with densities)
"""

# ╔═╡ cb9347f0-33df-4140-b0de-b5b0e79bae5a
md"""
Usually, it is convenient to solve the model in terms of quasi Fermi potentials. However, this way of definiting the model is not compatible with the setting in DifferentialEquations.jl. Thus, we use a DAE approach.
"""

# ╔═╡ fb971c90-ce59-416a-b25d-6bf71b38681b
begin
	function reaction!(f, u, node, data)

        f[ipsi]  = - (zn * (u[inn] - Nd) + zp * u[inp])
        ######################

        f[inn]   = 0.0
        f[inp]   = 0.0
        ########################
        etan     = zn * ( (u[iphin] - u[ipsi]) )
        etap     = zp * ( (u[iphip] - u[ipsi]) )

        nn       = statistics(etan)
        np       = statistics(etap)

        if data.inEquilibrium
            f[iphin] = u[iphin]#u[inn] - nn
            f[iphip] = u[iphip]#u[inp] - np
        else
            f[iphin] = u[inn] - nn
            f[iphip] = u[inp] - np
        end

    end

    function flux!(f, u, edge, cata)

        dpsi    = u[ipsi, 2] - u[ipsi, 1]
        f[ipsi] =  - εr * dpsi

        ####################################
        etan1   = zn * ( (u[iphin, 1] - u[ipsi, 1]))
        etan2   = zn * ( (u[iphin, 2] - u[ipsi, 2]))

        if statistics(etan1).value ≈ 0.0 || statistics(etan2).value ≈ 0.0
            bp, bm  = fbernoulli_pm(zn * dpsi)
        else
            bp, bm  = fbernoulli_pm(zn * dpsi + (etan2 - etan1) -log(statistics(etan2)) + log(statistics(etan1)) )
        end
        f[inn]  = - zn * ( bm * u[inn, 2] - bp * u[inn, 1])
        ####################################
        etap1   = zp * ( (u[iphip, 1] - u[ipsi, 1]))
        etap2   = zp * ( (u[iphip, 2] - u[ipsi, 2]))

        if statistics(etap1).value ≈ 0.0 || statistics(etap2).value ≈ 0.0
            bp, bm  = fbernoulli_pm(zp * dpsi)
        else
            bp, bm  = fbernoulli_pm(zp * dpsi + (etap2 - etap1) -log(statistics(etap2)) + log(statistics(etap1)) )
        end

        f[inp]  = - zp * ( bm * u[inp, 2] - bp * u[inp, 1])

    end

    function bcondition!(y, u, bnode, data)

        if bnode.region == 1
            barrier = barrierLeft
            Δu      = 0.0
        elseif bnode.region == 2
            barrier = barrierRight

        end

        if bnode.time == Inf
            Δu = 0.0
        else
            Δu = biasValue(bnode.time)
        end

        ############################################
        boundary_dirichlet!(y, u, bnode, species=ipsi, region=1, value=(- (barrierLeft  )) )
        boundary_dirichlet!(y, u, bnode, species=ipsi, region=2, value=(- (barrierRight ) + Δu) )
        ############################################

        etan     = zn * ( (u[iphin] - u[ipsi]) )
        etap     = zp * ( (u[iphip] - u[ipsi]) )
        nn       = statistics(etan)
        np       = statistics(etap)

        boundary_dirichlet!(y, u, bnode, species=iphin, region=1, value = 0.0)
        boundary_dirichlet!(y, u, bnode, species=iphip, region=1, value = 0.0)

        boundary_dirichlet!(y, u, bnode, species=iphin, region=2, value = Δu)
        boundary_dirichlet!(y, u, bnode, species=iphip, region=2, value = Δu)

        ############################################
        boundary_dirichlet!(y, u, bnode, species=inn, region=1, value= nn)
        boundary_dirichlet!(y, u, bnode, species=inn, region=2, value= nn)

        boundary_dirichlet!(y, u, bnode, species=inp, region=1, value= np)
        boundary_dirichlet!(y, u, bnode, species=inp, region=2, value= np)

    end

    function storage!(f, u, node, data)

        f[inn] = zn * u[inn]
        f[inp] = zp * u[inp]

    end
end

# ╔═╡ eec5fb8e-5d71-4eae-a3a8-8d4f454ef7ab
begin
	sys  = VoronoiFVM.System(grid, data             = data,
									flux            = flux!,
									reaction        = reaction!,
									bcondition      = bcondition!,
									storage         = storage!,
									unknown_storage =:sparse,
									species         = species)

	enable_species!(sys, iphin, regions)
	enable_species!(sys, iphip, regions)
	enable_species!(sys, inn,   regions)
	enable_species!(sys, inp,   regions)
	enable_species!(sys, ipsi,  regions)
end

# ╔═╡ 214ca96c-883f-4a62-8460-3c4f63f76876
begin
	inival = unknowns(sys)
	data.inEquilibrium = false
	inival[ipsi, :]  = solEQqF[ipsi, :]
	inival[iphin, :] = solEQqF[iphin, :]
	inival[iphip, :] = solEQqF[iphip, :]
	inival[inn, :]   = nnEQqF
	inival[inp, :]   = npEQqF
	solt = VoronoiFVMDiffEq.solve(sys, inival = inival, times=(0.0, endTimeSim))
	tval = solt.t
end

# ╔═╡ 57b47175-2dcb-4b3f-b2f8-ffe7f64085f2
md"""
# Implementation with Differential Equations
"""

# ╔═╡ cfeed287-79ae-4e83-80ed-bd96f9a0d2d3
begin
	desys  = VoronoiFVM.System(grid, data             = data,
									flux            = flux!,
									reaction        = reaction!,
									bcondition      = bcondition!,
									storage         = storage!,
									unknown_storage =:sparse,
									species         = species)

	enable_species!(desys, iphin, regions)
	enable_species!(desys, iphip, regions)
	enable_species!(desys, inn,   regions)
	enable_species!(desys, inp,   regions)
	enable_species!(desys, ipsi,  regions)
end

# ╔═╡ e1fa07cf-13ee-497f-ba1d-8fe347600a9c
endTimeSim

# ╔═╡ 6d3a96f5-558e-40c5-b25f-2455260dea8c
diffeqmethods=OrderedDict(
"QNDF2 (Like matlab's ode15s)" =>  QNDF2,
"Rodas5" => Rodas5,
"Rosenbrock23 (Rosenbrock)" => Rosenbrock23,
"FBDF" => FBDF,
"Implicit Euler" => ImplicitEuler
)


# ╔═╡ 50800988-e1c6-4642-ab6b-8556f274a174
md"""
method: $(@bind method Select([keys(diffeqmethods)...]))
"""

# ╔═╡ ec1a5fac-218d-4e11-823d-8480ad369976
begin
	problem = ODEProblem(desys, inival, (0.0, endTimeSim))
	tsol    = DifferentialEquations.solve(problem, diffeqmethods[method](), reltol=1.0e-3, abstol=1.0e-2)
	sol_DE  = reshape(tsol, sys)
	details(desys.history)
end

# ╔═╡ 547a6b59-558e-43d6-819b-5c648213403f
tsol.alg

# ╔═╡ 11f84a58-a135-42c4-a6a3-6739c9fcd229
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

# ╔═╡ 39badc0a-5c54-41ba-af21-eb3f9b70902f
myaside(md"""$(vis=GridVisualizer(dim=1, resolution=(380,250),Plotter=PlutoVista,legend=:lt);)""",top=0)

# ╔═╡ 72ab1020-14c7-4a40-a49e-e889c844ecb7
myaside(md"""$(vis2=GridVisualizer(dim=1, resolution=(380,250),Plotter=PlutoVista,legend=:lt);)""",top=250)

# ╔═╡ 13c09965-9880-4df7-80a0-fe0f1cdb5f08
myaside(md"""
t=$(@bind t  PlutoUI.Slider(range(0.0,endTimeSim,length=10001),default=0.0;show_value=true))
""",top=500)

# ╔═╡ 1792dcc6-e931-463a-8d39-2a99cae4359c
begin
	uqF = soltqF(t)
	scalarplot!(vis,grid, uqF[iphin,:], clear=true, color=:green,linewidth=5, markershape=:cross, linestyle=:dash,label="qF, Voronoi")
	scalarplot!(vis,grid, uqF[iphip,:], clear=false, color=:red, markershape=:cross, linestyle=:dash,label="qF, Voronoi")
	scalarplot!(vis,grid, uqF[ipsi,:], clear=false, color=:blue, markershape=:cross, linestyle=:dash, title = "t = $t, bias = $(biasValue(t))", label="qF, Voronoi")
	reveal(vis)
end

# ╔═╡ a0c0f5a5-41fc-4527-afd1-4d7487f9ba56
begin
	nntqF        = statistics.(zn .* (uqF[iphin,:] - uqF[ipsi,:]))
	nptqF        = statistics.(zp .* (uqF[iphip,:] - uqF[ipsi,:]))
	scalarplot!(vis2,grid, nntqF, clear=true, color=:green,markershape=:cross, linestyle=:dash,label="qF, Voronoi")
	scalarplot!(vis2,grid, nptqF, clear=false, color=:red, markershape=:cross, linestyle=:dash,title = "t = $t, bias = $(biasValue(t))", label="qF, Voronoi")
	reveal(vis2)
end

# ╔═╡ 05e280de-e4cd-4db0-bcb2-7095567944f9
begin
	u = solt(t)
	scalarplot!(vis,grid, u[iphin,:], clear=false, color=:darkgreen,markershape=:circle, linewidth=2, linestyle=:dash,label="dens, Voronoi")
	scalarplot!(vis,grid, u[iphip,:], clear=false, color=:darkred, markershape=:circle, linestyle=:dash,label="dens, Voronoi")
	scalarplot!(vis,grid, u[ipsi,:], clear=false, color=:darkblue, markershape=:circle, linestyle=:dash, title = "t = $t, bias = $(biasValue(t))", label="dens, Voronoi")
	reveal(vis)
end

# ╔═╡ c10a988b-8d9c-44eb-af4c-eda533953135
begin
	scalarplot!(vis2,grid, u[inn,:], clear=false, color=:darkgreen,markershape=:circle, linestyle=:dash,label="dens, Voronoi")
	scalarplot!(vis2,grid, u[inp,:], clear=false, color=:darkred, markershape=:circle, linestyle=:dash,title = "t = $t, bias = $(biasValue(t))", label="dens, Voronoi")
	reveal(vis2)
end

# ╔═╡ e6f9ea1c-49fe-4dbe-8db4-7832ee754c83
begin
	uDE = sol_DE(t)
	scalarplot!(vis,grid, uDE[iphin,:], clear=false, color=:gray,markershape=:none, linewidth=4, linestyle=:dot, label="")
	scalarplot!(vis,grid, uDE[iphip,:], clear=false)
	scalarplot!(vis,grid, uDE[ipsi,:], clear=false,label="DiffEq")
	reveal(vis)
end

# ╔═╡ 8171e174-e02d-4705-8ea7-2b1a9659df66
begin
	scalarplot!(vis2,grid, uDE[inn,:], clear=false, color=:gray,markershape=:none, linewidth=4, linestyle=:dot,  label="DiffEq")
	scalarplot!(vis2,grid, uDE[inp,:], clear=false, label="")
	reveal(vis2)
end

# ╔═╡ b4a73b4e-a322-447e-b88a-8a839660a78a
md"""
This notebook uses the documentation environment of the package and cannot be started outside of the `examples` subdirectory. To make it relocateable, make the next cell a markdown cell and restart the notebook. Please be aware that this way, version information is rebuilt from scratch.

"""

# ╔═╡ Cell order:
# ╠═5e728bfa-7571-11ed-09a3-e3fe7605ab0d
# ╟─fcdf6695-c96c-4579-ae57-53ec332a60ae
# ╠═2d35f805-34f1-4f99-b226-dad8924f5bd3
# ╠═88d0c7bd-6d8a-47e1-8064-3d1de6f46a3f
# ╠═73339ba6-71fc-4c25-9b9f-50a780d84ac2
# ╟─c2e90dcc-35a3-441f-8cee-2a1a9fb1e391
# ╠═80bf2ced-c5f9-4ae3-95f2-e450b85055a2
# ╟─a5bf9ba2-2e59-4616-8c17-d343bbc7ef71
# ╟─05d39cf1-85d8-41de-b50d-94ec00f21cc2
# ╟─5481508b-597a-4ff8-bce4-30612b2d8226
# ╟─49ef7e60-a27b-49fa-891e-bfad5c9db78c
# ╟─d9f62141-5856-4ecd-b1c5-b6f44597e48e
# ╠═9da3f051-8ab0-44ee-9c17-1453935487f1
# ╠═4e2b4ba4-ac05-4098-9f08-fd28e1cf9cb8
# ╟─36cc28cd-1737-41bd-8127-c784eb532443
# ╠═73917cbf-113b-4947-91d6-6cb0aac2dafc
# ╠═f46d0e11-e30c-420d-87ed-242235a2fdd9
# ╟─ee968c23-170b-4c66-b258-fbf4eebde2ce
# ╠═39238121-20f0-49d2-a6e9-96236d4a414e
# ╠═a306246a-ae46-4b8b-8539-df9a1e279f4f
# ╟─612906c4-13bd-4534-a184-ee892dcd6d03
# ╠═2d17f616-5dc6-4bfa-952d-c21f06e3ed6f
# ╠═1792dcc6-e931-463a-8d39-2a99cae4359c
# ╠═a0c0f5a5-41fc-4527-afd1-4d7487f9ba56
# ╟─485bb7b8-9b32-480a-a4d0-aab0426ece9b
# ╟─cb9347f0-33df-4140-b0de-b5b0e79bae5a
# ╠═fb971c90-ce59-416a-b25d-6bf71b38681b
# ╠═eec5fb8e-5d71-4eae-a3a8-8d4f454ef7ab
# ╠═214ca96c-883f-4a62-8460-3c4f63f76876
# ╠═05e280de-e4cd-4db0-bcb2-7095567944f9
# ╠═c10a988b-8d9c-44eb-af4c-eda533953135
# ╟─57b47175-2dcb-4b3f-b2f8-ffe7f64085f2
# ╠═cfeed287-79ae-4e83-80ed-bd96f9a0d2d3
# ╠═e1fa07cf-13ee-497f-ba1d-8fe347600a9c
# ╟─6d3a96f5-558e-40c5-b25f-2455260dea8c
# ╟─50800988-e1c6-4642-ab6b-8556f274a174
# ╠═ec1a5fac-218d-4e11-823d-8480ad369976
# ╠═547a6b59-558e-43d6-819b-5c648213403f
# ╠═e6f9ea1c-49fe-4dbe-8db4-7832ee754c83
# ╠═8171e174-e02d-4705-8ea7-2b1a9659df66
# ╠═11f84a58-a135-42c4-a6a3-6739c9fcd229
# ╠═39badc0a-5c54-41ba-af21-eb3f9b70902f
# ╠═72ab1020-14c7-4a40-a49e-e889c844ecb7
# ╠═13c09965-9880-4df7-80a0-fe0f1cdb5f08
# ╟─b4a73b4e-a322-447e-b88a-8a839660a78a
# ╠═e9f8097f-8dcc-4924-9eee-d8b0f49b0e4d
