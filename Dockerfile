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

FROM ubuntu:22.04

RUN apt-get update && apt-get install -y --no-install-recommends \
    make gcc g++ libc-dev libc6 gosu\
    && rm -rf /var/lib/apt/lists/*

LABEL program="SWORD2"
LABEL description="SWift and Optimized Recognition of protein Domains"
LABEL version="2.0.0"
LABEL maintainer="gabriel.cretin@u-paris.fr"

# Keep only necessary files from previous stage: conda env
COPY --from=mamba_build /venv /venv

# Set the working directory to /app
WORKDIR /app

# Copy sources to build the program
COPY install.sh install.sh
COPY bin/ bin/
COPY SWORD2.py SWORD2.py
COPY docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh

# Make the entrypoint script executable
RUN chmod +x /usr/local/bin/docker-entrypoint.sh

# Make sure the Conda environment's binaries are in the PATH
ENV PATH=/venv/bin:$PATH
ENV CONDA_PREFIX="/venv"

# Activate the Conda environment and run the install.sh script
# This step compiles all C/C++ dependencies
RUN bash install.sh

# Change ownership of the /app directory to root initially
RUN chown -R root:root /app

# Switch back to root to allow the entrypoint script to manage user creation
USER root

# Set the entrypoint to the entrypoint script
ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]

# Define a volume for the output directory
VOLUME /output