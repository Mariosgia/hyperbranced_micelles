#!/bin/bash

#Script that runs the analysis.exe for the chosen folder in 
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

wd=$(pwd)
qread=1
filee="phi_runs"
cd $filee
echo $filee
folders=(phi_runs/chosen*)



for j in ${folders[@]}; do


     folder=$j
     cd $folder
     
    

     if [ -f "./data/cnf_densities.out" ]; then
      echo "cnf exist in $folder"

      mkdir "analysis"
      cd "analysis"
      cp ${wd}/analysis.exe .
      cp "../data/cnf_densities.out" "cnf_densities.in"
      cp -r "../hyperbranched_data/" .
      cp "../input*" .
      cp ${wd}/analysis_submit.sbatch .

      sed -i "/qread/c\\${qread}     qread" input_general

      sbatch ./analysis_submit.sbatch
     else
	 echo "No cnf_densities.out file found in $folder"
     fi


   
     cd ..


done

cd ..

