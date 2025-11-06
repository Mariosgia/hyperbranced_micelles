#!/bin/bash
# Script to run different values of phibars 

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

size=15
grid=64
every_out=1000
every_show=10
accuracy_scf=0.000000000001
lambda=0.05
qread=1 # -1 Produce from program 1 feed in densities
um_rr=1.5
um_phi=0.2
umbrella=0
constraint=0
threads=128

dphi=0.00005
phst=0.00175
phen=0.00245

polymer_ensemble=2



phis=($(seq $phst $dphi $phen))

fol="phi_runs"




mkdir $fol



cp   run_parallel.sbatch main_hyper.exe input_interactions  input_general cnf_densities.in $fol
cp -r hyperbranched_data $fol
cd $fol




for (( i=0; i<${#phis[@]}; i++ )); do
    

     ph=${phis[i]}
     mkdir "phibar_${ph}"

     nam="phibar_${ph}"
     cp   run_parallel.sbatch main_hyper.exe input_interactions  input_general cnf_densities.in $nam
    
     cp -r hyperbranched_data $nam
     cd $nam

     sed -i "/gridx/c\\${grid}     gridx" input_general
     sed -i "/gridy/c\\${grid}     gridy" input_general
     sed -i "/gridz/c\\${grid}     gridz" input_general
     sed -i "/threads/c\\${threads}     threads" input_general
     sed -i "/lambda/c\\${lambda}     lambda" input_general
  
     sed -i "/sizex/c\\${size}d0     sizex" input_general
     sed -i "/sizey/c\\${size}d0     sizey" input_general
     sed -i "/sizez/c\\${size}d0     sizez" input_general
     sed -i "/phibar/c\\${ph}d0     phibar" input_general

     
     sed -i "/every_out/c\\${every_out}     every_out" input_general
     sed -i "/every_show/c\\${every_show}     every_show" input_general
     sed -i "/accuracy_scf/c\\${accuracy_scf}d0     accuracy_scf" input_general
     sed -i "/qread/c\\${qread}     qread" input_general
   
     sed -i "/umbrella/c\\${umbrella}     umbrella" input_general
     sed -i "/um_rr/c\\${um_rr}d0     um_rr" input_general
     sed -i "/um_phi/c\\${um_phi}d0     um_phi" input_general
     
     sed -i "/constraint/c\\${constraint}     constraint" input_general
     sed -i "/polymer_ensemble/c\\${polymer_ensemble}     polymer_ensemble" input_general

     echo "phibar = $ph "
     
    
     sbatch ./run_parallel.sbatch     

     cd ..


done

cd ..



