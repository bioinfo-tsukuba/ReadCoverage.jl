FROM julia:1.3.1

LABEL maintainer="Haruka Ozaki <harukao.cb@gmail.com>" license="MIT"

ENV JULIA_DEPOT_PATH=/opt/julia

RUN julia -e 'using Pkg; Pkg.add(PackageSpec(url="https://github.com/bioinfo-tsukuba/ReadCoverage.jl")); using ReadCoverage; Pkg.add("ArgParse"); using ArgParse'

RUN julia -e 'using Pkg; Pkg.add("PyPlot"); using PyPlot'

RUN mkdir -p /opt && curl https://raw.githubusercontent.com/bioinfo-tsukuba/ReadCoverage.jl/master/cli/run.jl --output /opt/run.jl && julia /opt/run.jl -h

RUN find /opt/julia/compiled/v1.3 -type f -exec chmod 666 {} + && chmod -R 777 /opt/julia/packages/ReadCoverage

RUN apt update && apt install -y procps