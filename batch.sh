sigJs=( 01 02 03 ); 
sigYs=( 05 10 20 ); 
sigTs=( 05 10 20 );
for i in "${sigJs[@]}"; 
  do for j in "${sigYs[@]}"; 
    do for k in "${sigTs[@]}"; 
      do sbatch --time=4:00:00 --mem=2G --ntasks=1 --nodes=1 --mail-user=nfoulk@umd.edu --mail-type=FAIL -J ""6""$i""_""$j""_""$k"" runner.jl 0.$i 0.$j 0.$k
    done; 
  done;
done; 
