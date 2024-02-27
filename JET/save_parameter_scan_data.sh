#!/bin/bash

# Read command-line arguments
shot_number=$1
dBB=$2
assimilation=$3

# Create folder name
folder_name="dBB_${dBB}_assim_${assimilation}"

# Create new folder
mkdir "$folder_name"

# Move specified files and folders into the newly created folder
mv output "$folder_name/"
mv output.h5 "$folder_name/"
mv __pycache__ "$folder_name/"
mv simulation_settings.log "$folder_name/"

# Ensure the target directory exists
target_dir="../JETresults/parameter_scans/${shot_number}"
mkdir -p "$target_dir"

# Move the newly created folder to the desired location within the shot number directory
mv "$folder_name" "$target_dir/"

echo "Operation completed. Folder ${folder_name} moved to ${target_dir}/"

