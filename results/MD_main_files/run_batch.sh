#!/bin/bash
count_running_jobs() {
    ps_count=$(ps aux | grep -c 'python3')
    ps_count=$(( $ps_count-1))
    echo "$ps_count"
}

# Function to find the next available folder name
get_next_folder_name() {
    i=1
    while [ -d "batch_$i" ]; do
        ((i++))
    done
    echo "batch_$i"
}

MAX_CPUS=$(nproc --all)
linear=$1 #Linear molecules
core=$2 #Core molecules
mu=$3 #Average hyperbranched molecules
dw=$4 #Polydispersity

X_VALUE=$(($MAX_CPUS * 2))  # Number of runs
echo $X_VALUE

SOURCE_PY_FILES="run_files/*.py"  
TARGET_PY_FILE="run.py"  
INITIAL_STATE_GENERATOR="initial_state_generator.py"



# Create X folders
for ((count=1; count<=$X_VALUE; count++)); do

    check=1
    while [ "${check}" -eq "1" ]; do
       if [ "$(count_running_jobs)" -le "$MAX_CPUS" ]; then
           folder_name=$(get_next_folder_name)
	   while [ -d  $folder_name ]; do
              sleep 1
              folder_name=$(get_next_folder_name)
           done
           mkdir "$folder_name"
           
	   cp $SOURCE_PY_FILES "$folder_name"
           cd "$folder_name"
           python3 "$INITIAL_STATE_GENERATOR" ${linear} ${core}
           nohup python3 "$TARGET_PY_FILE" ${mu} ${dw} > "output.log" 2>&1 &
           cd ..
           echo "New simulation job launched. $folder_name"
           check=0
       else
           #echo "Maximum number of simulation jobs reached. Sleeping for 5 seconds..."
           sleep "5"
           check=1
       fi
     done 	

done
