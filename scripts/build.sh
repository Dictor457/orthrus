#!/usr/bin/env bash
set -e
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build
