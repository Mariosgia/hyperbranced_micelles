#!/bin/bash

# Directory path
hm=$(pwd)
dirs=$( ls -d phi_runs/chosen*/phi*/stretch*/ )
for dir in $dirs; do
  if [ -d "$dir" ]; then
    echo "Processing directory: $dir"
    output_file="${dir}/runs_details.dat"
    # Create or truncate the output file
    > "${output_file}"


    str="# 1.Asp_ratio  2. Edif 3.Con_ene"
    echo "$str" >> "${output_file}"
    # Array to store folder names and corresponding values in the 7th column
    declare -A folder_values

    # Iterate through folders in the specified directory
    for folder in "${dir}/con"*; do
        if [ -d "${folder}" ]; then
            # Find details.dat file in each folder
            details_file="${folder}/data/details.dat"
            input_file="${folder}/input_general"
            echo ${folder}        
            # Check if details.dat file exists
            if [ -f "${details_file}" ]; then
                # Extract the last line of details.dat and append to box_details.dat
                last_line=$(tail -n 1 "${details_file}")
                
	        edif=$(echo "${last_line}" | awk '{print $18}')
                con_ene=$(echo "${last_line}" | awk '{print $19}')
                asp_ratio=$(sed -n '/asp_ratio/{s/d0//;s/ .*//;p;q;}' $input_file) #Extract mu_offset
                echo "$asp_ratio  $edif $con_ene" >> "${output_file}"
            else
	        echo "not found details file"
                          
            fi
        fi
    done

echo "Extraction and insertion complete. Results stored in ${output_file}"


  fi
done
exit





