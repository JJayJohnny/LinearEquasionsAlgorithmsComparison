#!/bin/sh
# Create build directory
mkdir -p build
# Clean up the build directory
rm -rf build/*
# Use GCC to compile the program
g++ -I/opt/vcpkg/installed/x64-linux/include -I/usr/include/python3.10 -DWITHOUT_NUMPY -std=c++17 -O2 -o build/metody_numeryczne src/*.cpp -lpython3.10
# Run the application
sh -c build/metody_numeryczne
