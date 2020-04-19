module ReadCoverage

using GenomicFeatures
using BioAlignments
using DelimitedFiles
using Printf
using Plots
using PyPlot
pyplot()

include("relative_genebodycoverage.jl")
include("readcovearge_bam.jl")
include("absolute_genebodycoverage.jl")
include("load_transcript.jl")
include("plot.jl")
include("utils.jl")

end # module