#!/bin/bash

# Directory containing PDB files
input_dir="./PDB_splitted"

# Output directory for MOL2 files
output_dir="./mol2_splitted"

# Loop through all PDB files in the input directory
for pdb_file in "$input_dir"/*.pdb; do
  # Extract the base filename without extension
  base_name=$(basename "$pdb_file" .pdb)
  
  # Define the output MOL2 file path
  mol2_file="$output_dir/${base_name}.mol2"
  
  # Convert PDB to MOL2 using obabel
  obabel "$pdb_file" -O "$mol2_file"
  
  # Check if the conversion was successful
  if [ $? -eq 0 ]; then
    echo "Converted: $pdb_file -> $mol2_file"
  else
    echo "Failed to convert: $pdb_file"
  fi
done
