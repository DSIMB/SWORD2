### BUILD-STAGE: conda environment
##################################

FROM continuumio/miniconda3 AS conda_build

RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc g++ libc-dev\
    && rm -rf /var/lib/apt/lists/*

# Install required python dependencies
COPY environment.yml .
RUN conda update -n base -c defaults conda && \
    conda env create -f environment.yml

# Install conda-pack:
RUN conda install -c conda-forge conda-pack -y

# Use conda-pack to create a standalone sword2 environment in /venv:
RUN conda-pack -j -1 -n sword2 -o /tmp/sword2_env.tar && \
    mkdir /venv && \
    cd /venv && \
    tar xf /tmp/sword2_env.tar && \
    rm /tmp/sword2_env.tar

# Finish unpacking the environment after unarchiving.
# Cleans up absolute prefixes in any remaining files
RUN /venv/bin/conda-unpack


### COMPILE-STAGE: Install and compile programs 
###############################################
FROM ubuntu:20.04 AS compile

COPY --from=conda_build /venv /venv

WORKDIR /sword2

# Copy sources to build the program
COPY bin/ bin/
COPY data/ data/
COPY SWORD2.py SWORD2.py

### RUNTIME-STAGE: Use the slimest image possible
#################################################

FROM ubuntu:20.04 as runtime

LABEL program="SWORD2"
LABEL description="SWift and Optimized Recognition of protein Domains"
LABEL version="1"
LABEL maintainer="gabriel.cretin@u-paris.fr"

WORKDIR /sword2

# Keep only necessary files from previous stages: conda env & sword2
COPY --from=conda_build /venv /venv
COPY --from=compile /sword2 /sword2

# Use `bash --login`:
SHELL ["/bin/bash", "--login", "-c"]

# Activate the conda environment by setting env paths
ENV PATH="/venv/bin:$PATH"
ENV CONDA_PREFIX="/venv"

ENTRYPOINT ["./SWORD2.py"]
CMD ["--help"]
