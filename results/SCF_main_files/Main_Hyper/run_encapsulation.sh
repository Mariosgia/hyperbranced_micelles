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

#############################################
####    Script assumes last molecule in   ###
#### hyperbranched_data/topology_data.dat ###
####       is a homopolymer (drug)        ###
#############################################
size=15
grid=64
every_out=1000
every_show=10
accuracy_scf=0.000000000001
lambda=0.05
qread=1 # -1 Produce from program 1 feed in densities
um_rr=0.0
um_phi=0.0
umbrella=0
constraint=0
threads=2

phihomopoly=0.000001 #Homopolymer phibar

num_mol=$(grep "Number_Molecules" hyperbranched_data/topology_data.dat | cut -d',' -f2 | tr -d ' ')

num_pol=$(($num_mol-1)) #Last molecule is model drug

sum_len_all=$(sed -n '4,$p' hyperbranched_data/topology_data.dat | grep -v "Mol-" | awk 'NF {print $3}' | awk '{sum += $1} END {print sum}')
len_drug=$(tail -n 1 hyperbranched_data/topology_data.dat | awk -F',' '{print $2}' | tr -d ' ')
sum_len_pol=$(($sum_len_all-$len_drug)) #Sum of len of polymers excluding drug

average_size_pol=$( awk -v var1=$sum_len_pol -v var2=$num_pol 'BEGIN {print var1 / (var2)}')




echo " sum len pol : $sum_len_pol and len of drug : $len_drug  average size pol : $average_size_pol" 
 
dphi=0.00005
phst=0.00050
phen=0.00200

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


#     
     wh=$( awk -v v1=$phihomopoly -v v2=$sum_len_pol -v v3=$len_drug  -v v4=$num_pol -v v5=$ph  'BEGIN {printf"%.15f\n", v1*v2/(v4*(v3*v5-v3*v1+v1*v2/v4))}') #weight of drug 
     wp=$( awk -v v1=$wh -v v2=$num_pol 'BEGIN {printf"%.15f\n", (1.0-v1)/v2 }') #weight of polymers 
     echo "phibar = $ph  wh = $wh wp = $wp "


     # Creating weights
    echo "#Weights n_i/n (ratio of number of chains) of each polymer found in the topology file in the same order" > hyperbranched_data/weights.dat

    # Append $wp to the file $num_pol times
    for ((j=1; j<=num_pol; j++)); do
      echo "$wp" >> hyperbranched_data/weights.dat
    done

    # Append $wd at the last line
    echo "$wh" >> hyperbranched_data/weights.dat


    
     sbatch ./run_parallel.sbatch     

     cd ..


done

cd ..



