#!/bin/bash

# Check if at least one argument was provided
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <FolderName> [\"Comment\"] [\"CustomOutputFolderName\"]"
    echo "FolderName must be one of: 85021, 85445, 85450, 85451, 85453, 85943"
    exit 1
fi

# Define the folder name from the first argument
folder_name="$1"

# Use the second argument as the comment, if provided
user_comment="${2:-No comment provided.}"

# Use the third argument as the custom output folder name, if provided
custom_output_folder_name="${3}"

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

# Define the source directory and create it if it doesn't exist
output_dir="./output"
mkdir -p "$output_dir"

# Move the output.h5 file and the __pycache__ folder to the output directory
mv -v "./output.h5" "$output_dir/"
mv -v "./__pycache__" "$output_dir/"
mv -v "./simulation_settings.log" "$output_dir/"

# Save the comment in a new file within the output directory
echo "$user_comment" > "$output_dir/move_comment.txt"

# Check if a custom output folder name was provided
if [ -n "$custom_output_folder_name" ]; then
    new_output_dir_name="output_$custom_output_folder_name"
else
    # Get current date and time to rename the output directory if no custom name provided
    current_datetime=$(date "+%Y-%m-%d_%H-%M-%S")
    new_output_dir_name="output_$current_datetime"
fi

# Rename the output directory
mv -v "$output_dir" "$new_output_dir_name"

# Define the final destination directory within JETresults and ensure it exists
jet_results_dir="../../../../../mnt/DISK4/christiang/resultat/$folder_name"
mkdir -p "$jet_results_dir"

# Move the newly named output directory to the specified JETresults subdirectory
mv -v "$new_output_dir_name" "$jet_results_dir/"

echo "Files have been moved successfully to $jet_results_dir/$new_output_dir_name."

