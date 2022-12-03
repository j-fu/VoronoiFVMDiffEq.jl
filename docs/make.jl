using LinearAlgebra, Pkg
using Documenter, VoronoiFVMDiffEq, PlutoSliderServer, SciMLBase




function make_all()
                 
    ENV["VORONOIFVM_CHECK_ALLOCS"]="false"
    notebook_html_dir  = joinpath(@__DIR__,"src","nbhtml")

    with_notebooks=true
    if with_notebooks
        
        #
        # Run notebooks
        #
        notebooks=[
            "1D nonlinear diffusion" => "Notebook101_NonlinearDiffusion1D.jl",
            "1D wave equation" => "Notebook102_WaveEquation1D.jl",
            "1D nonlinear storage" => "Notebook103_NonlinearStorage1D.jl",
            "1D/2D brusselator" => "Notebook203_Brusselator.jl"
        ]
        
        notebookjl=last.(notebooks)
        notebookmd=[]
        export_directory(joinpath(@__DIR__,"..","examples"),
                         notebook_paths=notebookjl,
                         Export_output_dir=joinpath(notebook_html_dir),
                         Export_offer_binder=false)
        
        # generate frame markdown for each notebook
        for notebook in notebookjl
            base=split(notebook,".")[1]
            mdstring=
"""
##### [$(base).jl](@id $(base))
[Download](https://github.com/j-fu/VoronoiFVMDiffEq.jl/blob/main/examples/$(notebook))
this [Pluto.jl](https://github.com/fonsp/Pluto.jl) notebook.

```@raw html
<iframe style="height:15000px" width="100%" src="../$(base).html"> </iframe>
```
"""
            mdname=base*".md"
            push!(notebookmd,joinpath("nbhtml",mdname))
            io=open(joinpath(notebook_html_dir,mdname),"w")
            write(io,mdstring)
            close(io)
        end     
        
        notebooks=first.(notebooks).=> notebookmd
        pushfirst!(notebooks, "About the notebooks"=> "notebooks.md")
    else
        notebooks=nothing
    end
    
    makedocs(
        sitename="VoronoiFVMDiffEq.jl",
        modules = [VoronoiFVMDiffEq],
        clean = false, 
        doctest = true,
        authors = "J. Fuhrmann",
        repo="https://github.com/j-fu/VoronoiFVM.jl",
        pages=[ 
            "Home"=>"index.md",
            "Tutorial Notebooks" => notebooks
        ]
    )

    if with_notebooks
        rm(notebook_html_dir,recursive=true)
    end
    
    if !isinteractive()
        deploydocs(repo = "github.com/j-fu/VoronoiFVMDiffEq.jl.git",devbranch="main")
    end
end

make_all()

