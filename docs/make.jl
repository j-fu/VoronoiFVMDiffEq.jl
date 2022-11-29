using LinearAlgebra, Pkg
using Documenter, VoronoiFVMDiffEq, PlutoSliderServer




function make_all()
    wd=pwd()
    Pkg.activate(joinpath("..","test"))
    Pkg.instantiate()
    Pkg.activate(wd)
                 
    ENV["VORONOIFVM_CHECK_ALLOCS"]="false"
    generated_examples=[]
    notebook_html_dir  = joinpath(@__DIR__,"src","nbhtml")
    #
    # Run notebooks
    #
    notebooks=[
        "1D nonlinear diffusion" => "Notebook101_NonlinearDiffusion1D.jl",
        "1D wave equation" => "Notebook102_WaveEquation1D.jl",
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
    
        
    makedocs(
        sitename="VoronoiFVMDiffEq.jl",
        modules = [VoronoiFVM],
        clean = false, 
        doctest = true,
        authors = "J. Fuhrmann",
        repo="https://github.com/j-fu/VoronoiFVM.jl",
        pages=[ 
            "Home"=>"index.md",
            "Tutorial Notebooks" => notebooks
        ]
    )

    rm(notebook_html_dir,recursive=true)
    
    if !isinteractive()
        deploydocs(repo = "github.com/j-fu/VoronoiFVMDiffEq.jl.git",devbranch="main")
    end
end

make_all()

