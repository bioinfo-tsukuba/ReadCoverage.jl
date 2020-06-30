export plot_absolute_coverage, plot_relative_coverage, plot_read_coverage

"""
    plot_absolute_coverage(abcov::Array{Float64,1}; out_path="")

Plot the absolute gene body coverage.
"""
function plot_absolute_coverage(abcov::Array{Float64,1}; out_path="")
    x = collect(1:length(abcov)) .* 100
    pyplot()
    p = Plots.plot(x, abcov, xlabel = "Distance from the \$3^\\prime\$ end [bp] (\$5^\\prime\$ -> \$3^\\prime\$)", label = "", title = "", xflip=true);
    if out_path != ""
        Plots.savefig(out_path)
    else
        p
    end
end


"""
    plot_relative_coverage(relcov::Array{Float64,1}; out_path="")

Plot the relative gene body coverage.
"""
function plot_relative_coverage(relcov::Array{Float64,1}; out_path="")
    x = collect(1:length(relcov)) .* 100
    pyplot()
    p = Plots.plot(1:100,  (relcov.-minimum(relcov))./(maximum(relcov)-minimum(relcov)), 
            xlabel = "Gene body percentile (\$5^\\prime\$ -> \$3^\\prime\$)", label = "", title = "");
    if out_path != ""
        Plots.savefig(p, out_path)
    else
        p
    end
end


"""
    plot_read_coverage(cov::Array{Int64,1}; out_path="")

Plot read coverage.
"""
function plot_read_coverage(cov::Array{Int64,1}; out_path="")
    x = collect(1:length(cov)) .* 100
    pyplot()
    p = Plots.plot(x, cov, xlabel = "Genomic coordinate [bp]", label = "", title = "");
    if out_path != ""
        Plots.savefig(p, out_path)
    else
        p
    end
end
