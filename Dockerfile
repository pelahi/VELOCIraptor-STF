FROM ubuntu:18.04

# DYNAMIC BRANCH CAN BE OVERWRITTEN DURING RUNTIME
ENV BRANCH=feature/threadpool-struct

# NEEDS TO BE UPDATED FOR SPECIFIC VERSIONS
RUN apt update && \
    apt install -y g++ libomp-dev libgsl-dev libhdf5-serial-dev git cmake

WORKDIR /home/ubuntu/VELOCIraptor-STF

# INITIALISE PROJECT DIRECTORY
# RUN git clone https://github.com/pelahi/VELOCIraptor-STF.git && \
#     cd VELOCIraptor-STF && git checkout ${BRANCH} && git submodule update --init --recursive
COPY . .

# BUILD BINARY
RUN mkdir build && cd build && cmake .. && make all
