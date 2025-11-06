#!/bin/bash
# Script to run different values of size_x until the one with the lowest free energy is found 

#0 if false 1 if true
lessthan() {
   local var1="$1"
   local var2="$2"

   local var=$(awk 'BEGIN {OFMT = "%0.5f"; a = '$var1' ; b = '$var2' ; if (a < b) print 1 ; else print 0}')  
   echo $var

}

#0 if false 1 if true
contains() {
    local var1="$1"   
    local arr=($2)
   
    local leh=${#arr[@]}

    for (( j=0; j<${leh}; j++ )); 
    do
        var2="${arr[j]}"
        local tst=$(awk 'BEGIN {OFMT = "%0.5f"; a = '$var1' ; b = '$var2' ; if (a == b) print 1 ; else print 0}')   
        if [[ $tst -eq 1 ]] 
        then
            break
        fi
    done
    echo $tst
}


############################################
########      BEGINNING OF MAIN      #######
############################################


qread=1 # -1 Produce from program 1 feed in densities
constraint=3
kapcon=1
ensemble=1 #Canonical ensemble 


dasp=0.2
aspmax=1.0
aspmin=6.0



dir=$( ls -d phi_runs/chosen*/phi*/ )

# Check if the directory exists
if [ -n "$dir" ]; then
  cd "$dir"
else
  echo "No directory found in chosen case."
  exit
fi



if [ -f "./data/cnf_densities.out" ]; then
    echo "cnf out exists"
    else
     echo "No cnf_densities.out file found. Aborting"
     exit
fi


mkdir "stretch_runs"

phibar=$(tail  -n 1 data/details.dat  | awk '{ print $8}' )

asps=($(seq $dasp $aspmin $aspmax))


for j in ${asps[@]}; do

  folder="asp_ratio_$j"

  asp=$j

  test=$(awk 'BEGIN {OFMT = "%0.5f"; a = '$j' ; b = 1 ; if (a == b) print 1 ; else print 0}')

  if [ $test -eq 1 ]; then
    cons=0
  else
    cons=$constraint
  fi
  


  echo "asp : $asp "

  mkdir "stretch_runs/$folder"
  cd "stretch_runs/$folder"

  cp "../../data/cnf_densities.out" "./cnf_densities.in"
  cp -r "../../hyperbranched_data" .
  rm  hyperbranched_data/weights.dat

  if [ -e "../../analysis/new_weights/weights.dat" ]; then 
     cp ../../analysis/new_weights/weights.dat hyperbranched_data
  else 
     echo "New weights in analysis don't exist" 
     exit
  fi

  cp "../../input_general" "../../input_interactions" .
  cp  "../../main_hyper.exe"  .
  cp  "../../run_parallel.sbatch"  .


  sed -i "/qread/c\\${qread}     qread" input_general

  sed -i "/polymer_ensemble/c\\${ensemble}     polymer_ensemble" input_general 
  sed -i "/constraint/c\\${cons}     constraint" input_general 
  sed -i "/kapcon/c\\${kapcon}d0     kapcon" input_general 
  sed -i "/asp_ratio/c\\${asp}d0     asp_ratio" input_general 
  sed -i "/phibar/c\\${phibar}d0     phibar" input_general 





  sbatch ./run_parallel.sbatch

  cd ../..


done




