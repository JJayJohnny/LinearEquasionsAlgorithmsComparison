# Get the base Ubuntu image from Docker Hub
FROM ubuntu:latest

# Update apps on the base image
RUN apt-get -y update && apt-get install -y

# Install vcpkg
RUN apt-get -y install build-essential pkg-config git curl unzip zip tar
RUN cd /opt && git clone https://github.com/microsoft/vcpkg && ./vcpkg/bootstrap-vcpkg.sh
ENV PATH="${PATH}:/opt/vcpkg"

#install python with required libraries
RUN apt-get -y install python3 python3-dev python3-pip
RUN pip3 install matplotlib
#RUN pip3 install numpy
#RUN apt-get -y install python3-matplotlib

# Install matplotlib library
RUN vcpkg install matplotlib-cpp

# Specify the working directory
WORKDIR /usr/src/metody_numeryczne

#RUN ln -s  /usr/local/lib/python3.10/dist-packages/numpy/core/include/numpy /usr/include/numpy

#############
# Run the program
CMD ["sh", "entrypoint.sh"]