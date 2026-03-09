### BUILD-STAGE: Compile Rust SWORD2 binary
##########################################
FROM rust:1.76-bookworm AS rust_build

# Set the working directory to /app
WORKDIR /app

# Copy the rust source code
COPY sword2-rs/ sword2-rs/

# Build the release binary
RUN cd sword2-rs && cargo build --release

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

# Set the working directory to /app
WORKDIR /app

# Copy the compiled Rust binary from the build stage
COPY --from=rust_build /app/sword2-rs/target/release/sword2 /usr/local/bin/sword2

# Copy sources to build the internal tools and data files needed
COPY install.sh install.sh
COPY bin/ bin/
COPY docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh

# Make the entrypoint script executable
RUN chmod +x /usr/local/bin/docker-entrypoint.sh

# Run the install.sh script
# This step compiles all C/C++ dependencies in bin/
RUN bash install.sh

# Change ownership of the /app directory to root initially
RUN chown -R root:root /app

# Switch back to root to allow the entrypoint script to manage user creation
USER root

# Set the entrypoint to the entrypoint script
ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]

# Define a volume for the output directory
VOLUME /output