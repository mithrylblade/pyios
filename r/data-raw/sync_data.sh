#!/bin/sh
# Copy the master CSVs from the Python package into the R package.
set -e
cd "$(dirname "$0")/.."
cp ../src/pyios/data/reference_equations.csv ../src/pyios/data/reference_metadata.csv inst/extdata/
echo "R package data synced from src/pyios/data"
