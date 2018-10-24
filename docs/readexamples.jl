baseURL = "https://github.com/alan-turing-institute/PDSampler.jl/blob/master/"
tests   = readdir("../test/")
for name in tests
    if occursin(Regex("^ex_"), name) && name[end-2:end]==".jl"
        #########################
        pathin  = "../test/$name"
        URL     = baseURL * "test/$name"
        mdname  = "$(name[1:end-3]).md"
        pathout = "src/examples/$mdname"
        open(pathin) do fi
            open(pathout, "w") do fo
                skip    = true
                iscode  = false
                llclose = false
                for ln in eachline(fi)
                    # Chomp markers
                    if occursin(Regex("^#@startexample"),ln)
                        write(fo, "#"*match(Regex("^#@startexample(.*)"),ln)[1]*"\n\n")
                        write(fo, "(*the code for this example can be found "*
                                  "[here]($URL), note that the doc rendered "* "here was automatically generated, if you "*
                                  "want to fix it, please do it in the "*
                                  "julia code directly*)\n\n")
                        skip = false
                    elseif occursin(Regex("^#@endexample"),ln)
                        break
                    elseif occursin(Regex("^#="),ln)
                        if iscode
                            # close the block
                            write(fo, "\n```\n")
                            iscode = false
                        end
                    elseif occursin(Regex("^=#"),ln)
                        iscode  = true
                        llclose = true
                    else
                        if iscode && llclose
                            write(fo, "\n```julia\n")
                            llclose = false
                        end
                        skip ? nothing : write(fo, ln*"\n")
                    end
                end
            end
        end
    end
end
