using PDMP, Base.Test

@testset "gaussian"    begin include("gaussian_test.jl")   end
@testset "logreg"      begin include("logreg_test.jl")     end
@testset "pmf"         begin include("pmf_test.jl")        end

@testset "geometry"    begin include("geometry_test.jl")   end
@testset "ippsampler"  begin include("ippsampler_test.jl") end
@testset "kernels"     begin include("kernels_test.jl")    end
@testset "path"        begin include("path_test.jl")       end
@testset "simulate"    begin include("simulate_test.jl")   end

@testset "l_event"     begin include("local_event_test.jl")        end
@testset "l_fgraph"    begin include("local_factorgraph_test.jl")  end
@testset "l_sim"       begin include("local_simulate_test.jl")     end
