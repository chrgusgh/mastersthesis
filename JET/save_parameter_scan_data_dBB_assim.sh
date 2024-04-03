#!/bin/bash

# Read command-line arguments
shot_number=$1
dBB=$2
assim=$3

# Correctly create folder name using B_factor and Arfrac
folder_name="dBB_${dBB}_assim_${assim}"

# Create new folder
mkdir "$folder_name"

# Move specified files and folders into the newly created folder
mv output "$folder_name/"
mv output.h5 "$folder_name/"
mv __pycache__ "$folder_name/"
mv simulation_settings.log "$folder_name/"
mv I_RE.txt "$folder_name/"
mv I_Ohm.txt "$folder_name/"
mv I_hot.txt "$folder_name/"
mv I_tot.txt "$folder_name/"
mv t.txt "$folder_name/"
mv tau_CQ.txt "$folder_name/"

# Ensure the target directory exists
target_dir="../../../../../mnt/DISK4/christiang/resultat/parameter_scans/${shot_number}"
mkdir -p "$target_dir"

# Move the newly created folder to the desired location within the shot number directory
mv "$folder_name" "$target_dir/"

echo "Operation completed. Folder ${folder_name} moved to ${target_dir}/"

