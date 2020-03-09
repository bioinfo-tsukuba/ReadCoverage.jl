

# Ref: GenomicFeatures.jl https://biojulia.net/GenomicFeatures.jl/latest/io/bed/
# Ref strand!! https://github.com/BioJulia/GenomicFeatures.jl/blob/8fc34ff680f5e742e25d2bd3d4722cb33fbe3cd5/src/strand.jl
# Future TODO: Test
# Future TODO: How to make subset of 'expressed' transcripts (Give list of transcript id?)
# Ref: Plot https://qiita.com/yuifu/items/b2031d0367da878fddae#bioalignments
# Future TODO: Consider Stranded RNA-seq

using GenomicFeatures
using BioAlignments
using DelimitedFiles
using Printf

include("relative_genebodycoverage.jl")
include("readcovearge_bam.jl")
include("absolute_genebodycoverage.jl")
include("load_transcript.jl")


