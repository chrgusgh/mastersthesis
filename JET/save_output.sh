#!/bin/bash

# Check if at least one argument was provided
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <FolderName> [\"Comment\"]"
    echo "FolderName must be one of: 85021, 85445, 85450, 85451, 85453, 85943"
    exit 1
fi

folder_name="$1"
user_comment="${2:-No comment provided.}"

# Validate the folder name
case "$folder_name" in
    85021|85445|85450|85451|85453|85943)
        ;;
    *)
        echo "Invalid FolderName: $folder_name"
        echo "FolderName must be one of: 85021, 85445, 85450, 85451, 85453, 85943"
        exit 1
        ;;
esac

# Extract dBB_cold, assimilation, and t_TQ values from simulation_settings.log
dBB_cold=$(grep "dBB_cold" simulation_settings.log | cut -d " " -f3)
assimilation=$(grep "assimilation" simulation_settings.log | cut -d " " -f2)
t_TQ=$(grep "t_TQ" simulation_settings.log | cut -d " " -f2)

# Create a custom folder name based on the extracted values
custom_folder_name="dBB_${dBB_cold}_assim_${assimilation}_tTQ_${t_TQ}"

# Define the source directory and create it if it doesn't exist
output_dir="./output"
mkdir -p "$output_dir"

# Move the output.h5 file, the __pycache__ folder, and simulation_settings.log to the output directory
mv -v "./output.h5" "$output_dir/"
mv -v "./__pycache__" "$output_dir/"
mv -v "./simulation_settings.log" "$output_dir/"

# Save the comment in a new file within the output directory
echo "$user_comment" > "$output_dir/move_comment.txt"

# Use the custom folder name to rename the output directory
new_output_dir_name="output_${custom_folder_name}"

# Rename the output directory
mv -v "$output_dir" "$new_output_dir_name"

# Define the final destination directory within JETresults and ensure it exists
jet_results_dir="../JETresults/$folder_name/"
mkdir -p "$jet_results_dir"

# Move the newly named output directory to the specified JETresults subdirectory
mv -v "$new_output_dir_name" "$jet_results_dir/"

echo "Files have been moved successfully to $jet_results_dir/$new_output_dir_name."

