FROM julia:1.8.5

LABEL maintainer="Haruka Ozaki <harukao.cb@gmail.com>" license="MIT"

ENV JULIA_DEPOT_PATH=/opt/julia

RUN apt update && apt install -y g++ procps

RUN cd /opt && \
	curl -L https://github.com/bioinfo-tsukuba/ReadCoverage.jl/archive/refs/tags/v0.1.3.tar.gz --output v0.1.3.tar.gz && \
	tar xvzf v0.1.3.tar.gz && \
	rm v0.1.3.tar.gz && \
	cp /opt/ReadCoverage.jl-0.1.3/cli/run.jl /opt/ && \
	julia --project=/opt/ReadCoverage.jl-0.1.3 -e 'using Pkg; Pkg.instantiate(); Pkg.precompile(); Pkg.add("ArgParse"); using ArgParse; Pkg.add("PackageCompiler"); using PackageCompiler; create_sysimage([:ArgParse, :GenomicFeatures, :XAM, :BED, :BioAlignments, :DelimitedFiles, :Printf, :Plots, :Statistics, :GZip, :ReadCoverage], sysimage_path="/opt/ReadCoverage.so", precompile_execution_file="/opt/run.jl")'

# RUN find /opt/julia/compiled/v1.8 -type f -exec chmod 666 {} + && chmod -R 777 /opt/julia/packages/ReadCoverage

ENTRYPOINT ["julia", "--project=/opt/ReadCoverage.jl-0.1.3", "-J", "/opt/ReadCoverage.so", "/opt/run.jl"]

