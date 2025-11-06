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
constraint=1
kapcon=1
con_rr=1.5
con_shift=10 

accuracy_scf=0.00000000000001
#lambda=0.01

dasph=0.1
asphmax=3.0
asphmin=1.0



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


mkdir "asph_runs"



con_asph=($(seq $asphmin $dasph $asphmax))

for j in ${con_asph[@]}; do

  folder="apsh_$j"

  apsh=$j

  test=$(awk 'BEGIN {OFMT = "%0.5f"; a = '$j' ; b = 1 ; if (a == b) print 1 ; else print 0}')

  if [ $test -eq 1 ]; then
    cons=0
  else
    cons=$constraint
  fi
  


  echo "apsh : $apsh"
  mkdir "asph_runs/$folder"
  cd "asph_runs/$folder"

  cp "../../data/cnf_densities.out" "./cnf_densities.in"
  cp -r "../../hyperbranched_data" .
  cp "../../input_general" "../../input_interactions" .
  cp  "../../main_hyper.exe"  .
  cp  "../../run_parallel.sbatch"  .


  sed -i "/qread/c\\${qread}     qread" input_general
  sed -i "/accuracy_scf/c\\${accuracy_scf}     accuracy_scf" input_general

  sed -i "/constraint/c\\${cons}     constraint" input_general 
  sed -i "/kapcon/c\\${kapcon}d0     kapcon" input_general 
  sed -i "/con_rr/c\\${con_rr}d0     con_rr" input_general 
  sed -i "/con_shift/c\\${shift}d0     con_shift" input_general 
  sed -i "/asp_ratio/c\\${apsh}d0     asp_ratio" input_general 




  sbatch ./run_parallel.sbatch

  cd ../..


done




