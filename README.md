# ReadCoverage.jl

ReadCoverage.jl is a fast tool to calculate absolute and relative gene body coverage of bulk/single-cell RNA-seq data.

## Installation

```julia
pkg> add https://github.com/bioinfo-tsukuba/ReadCoverage.jl
```

## Usage

```julia
using ReadCoverage

path_bam = "examples/data/RamDA_72h_A09.uniq.q40.bam"
path_bed12 = "examples/data/gencode.vM15.primary_assembly.annotation.protein_coding.head.bed"
out_prefix = "examples/out/RamDA_72h_A09.uniq.q40"
```

### Read coverage

```julia
chrom = "chr19"
leftpos = 3205000
rightpos = 3207000
cov = readcoverage_bam(path_bam, chrom, leftpos, rightpos)
cov = readcoverage_bam(path_bam, chrom, leftpos, rightpos, output_prefix=out_prefix)
plot_read_coverage(cov)
```

### Absolute genebody coverage

```julia
abcov = absolute_genebodycoverage(path_bam, path_bed12)
abcov = absolute_genebodycoverage(path_bam, path_bed12, output_prefix=out_prefix, bin_size=100)
plot_absolute_coverage(abcov)
```


### Relative genebody coverage

```julia
relcov = relative_genebodycoverage(path_bam, path_bed12);

max_depth = 0
relcov = relative_genebodycoverage(path_bam, path_bed12, output_prefix=out_prefix, max_depth=max_depth);
plot_relative_coverage(relcov)
```

### Variability of Depth of Coverage (VDoC) score

```julia
array_path_bam = [
    "examples/data/RamDA_72h_A09.uniq.q40.bam",
    "examples/data/RamDA_72h_A10.uniq.q40.bam",
    "examples/data/RamDA_72h_A11.uniq.q40.bam"
]
path_bed12 = "examples/data/gencode.vM15.primary_assembly.annotation.protein_coding.head.bed"

array_transcript_name, array_transcript_length, array_vdoc = calc_vdoc(array_path_bam, path_bed12)
```


## CLI (Commad line interface)

```julia
pkg> add ArgParse
```

```bash
julia cli/run.jl relcov <BAM> <BED12> <output_prefix>
```

## Docker

The Docker image is available on [Docker Hub](https://hub.docker.com/r/yuifu/readcoverage.jl).

You can run ReadCoverage.jl via CLI.

### `relcov`

```bash
docker run yuifu/readcoverage.jl:0.1.2 \
  julia /opt/run.jl relcov <BAM> <BED12> <output_prefix>
```

### `vdoccov`

```julia
docker run yuifu/readcoverage.jl:0.1.2 \
  julia /opt/run.jl vdoc <BED12> <output_prefix> <BAM1> <BAM2> ...
```
