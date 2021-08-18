#!/usr/bin/env bash

notloaded=true
mod="julia"; for i in `module list`; do [[ $i == *$mod* ]] && notloaded=false; done
[[ notloaded ]] && module load julia; 

sigJs=( 01 02 03 ); 
sigYs=( 05 10 20); 
sigTs=( 05 10 20 );
lengs=( 4 );
sings=( true );
for sing in "${sings[@]}";
do
  if [ $sing == true ]
  then
    singstr='s'
  else
    singstr='up'
  fi
  for ind in {1..10};
  do
    for leng in "${lengs[@]}";
    do
      for i in "${sigJs[@]}"; 
        do for j in "${sigYs[@]}"; 
          do for k in "${sigTs[@]}"; 
            do 
        filename=jdata$ind/mathematica_${leng}_${singstr}_β0.01_γ0.10_σJ0.${i}_σγ0.${j}_στ0.${k}_01000_0.03;
        #echo $filename
        if [ ! -s $filename ]
        then
          sbatch --time=2:00:00 --mem=2G --ntasks=1 --nodes=1 --mail-user=nfoulk@umd.edu --mail-type=FAIL -J ""$leng""$singstr""$i""_""$j""_""$k"".""$ind"" runner.jl --singlet=$sing --length=$leng --sigmas=0.$i,0.$j,0.$k --nReals=1000 --dirspec=$ind
        fi 
          done; 
        done;
      done;
    done;
  done;
done
