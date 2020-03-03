FROM ubuntu:18.04

# DYNAMIC BRANCH CAN BE OVERWRITTEN DURING RUNTIME
ENV BRANCH=""

# NEEDS TO BE UPDATED FOR SPECIFIC VERSIONS
RUN apt update && \
    apt install -y g++ libomp-dev libgsl-dev libhdf5-serial-dev git cmake

WORKDIR /home/ubuntu/

# COPY PROJECT DIRECTORY
# RUN git clone https://github.com/pelahi/VELOCIraptor-STF.git 
COPY . /home/ubuntu/VELOCIraptor-STF

WORKDIR /home/ubuntu/VELOCIraptor-STF

#RUN git checkout ${BRANCH} 
#&& git submodule update --init --recursive

# BUILD BINARY
RUN mkdir build && cd build && cmake .. && make all
