#!/usr/bin/env bash

# URL of the STL file
stl_url="https://openlb.net/stl/kit_campus.stl"

# Output file name
output_file="kit_campus.stl"

# Download the STL file using wget
wget -O "$output_file" "$stl_url"

# Check if the download was successful
if [ $? -eq 0 ]; then
    echo "STL file downloaded successfully to: $output_file"
else
    echo "Error: Unable to download STL file."
fi
