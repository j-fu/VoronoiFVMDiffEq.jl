using Test
using Pluto
using Pkg

# Activate assembly loop allocation checking
# as default.
ENV["VORONOIFVM_CHECK_ALLOCS"]="true"

modname(fname)=splitext(basename(fname))[1]

function testnotebook(input)
    # run notebook and check for cell errors
    session = Pluto.ServerSession();
    session.options.server.disable_writing_notebook_files=true
    session.options.server.show_file_system=false
    session.options.server.launch_browser=false
    session.options.server.dismiss_update_notification
    session.options.evaluation.capture_stdout=false
    session.options.evaluation.workspace_use_distributed=false # this makes it fast
    
    wd=pwd()
    @time notebook = Pluto.SessionActions.open(session, input; run_async=false)
    cd(wd)
    errored=false
    for c in notebook.cells
        if c.errored
            errored=true
            @error "Error in  $(c.cell_id): $(c.output.body[:msg])\n $(c.code)"
        end
    end
    !errored
end

function run_all_tests()
    wd=pwd()
    Pkg.instantiate()
    Pkg.activate(wd)

    ENV["VORONOIFVM_CHECK_ALLOCS"]="true"
    notebooks=["Notebook101_NonlinearDiffusion1D.jl",
               "Notebook102_WaveEquation1D.jl",
               "Notebook103_NonlinearStorage1D.jl",
               "Notebook203_Brusselator.jl"]
    
    @testset "notebooks" begin
        for notebook in notebooks
            #            include(joinpath(@__DIR__,"..","pluto-examples",notebook))
            @info "notebook $(notebook):"
            @test testnotebook(joinpath(@__DIR__,"..","examples",notebook))
            @info "notebook $(notebook) ok"
        end
    end
end

run_all_tests()
