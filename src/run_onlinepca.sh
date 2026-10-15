#!/bin/bash
mkdir -p /tmp/julia_depot
export JULIA_DEPOT_PATH="/tmp/julia_depot:/julia_packages"
export JULIA_HISTORY="/dev/null"
/usr/local/julia/bin/julia src/run_onlinepca.jl "$@"
