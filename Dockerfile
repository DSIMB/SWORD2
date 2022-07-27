### BUILD-STAGE: conda environment
##################################

FROM condaforge/mambaforge:4.12.0-2 AS mamba_build

RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc g++ libc-dev\
    && rm -rf /var/lib/apt/lists/*

# Install required python dependencies
COPY environment.yml .
RUN mamba env create -f environment.yml && \
    conda clean -afy

# Install conda-pack:
RUN mamba install -c conda-forge conda-pack -y && \
    conda clean -afy

# Use conda-pack to create a standalone sword2 environment in /venv:
RUN conda-pack -j -1 -n sword2 -o /tmp/sword2_env.tar && \
    mkdir /venv && \
    cd /venv && \
    tar xf /tmp/sword2_env.tar && \
    rm /tmp/sword2_env.tar

# Finish unpacking the environment after unarchiving.
# Cleans up absolute prefixes in any remaining files
RUN /venv/bin/conda-unpack

### Install and run 
###################

FROM ubuntu:20.04

LABEL program="SWORD2"
LABEL description="SWift and Optimized Recognition of protein Domains"
LABEL version="1"
LABEL maintainer="gabriel.cretin@u-paris.fr"

WORKDIR /sword2

# Keep only necessary files from previous stage: conda env
COPY --from=mamba_build /venv /venv

# Copy sources to build the program
COPY bin/ bin/
COPY data/ data/
COPY SWORD2.py SWORD2.py

# Use `bash --login`:
SHELL ["/bin/bash", "--login", "-c"]

# Activate the conda environment by setting env paths
ENV PATH="/venv/bin:$PATH"
ENV CONDA_PREFIX="/venv"

ENTRYPOINT ["./SWORD2.py"]
CMD ["--help"]
