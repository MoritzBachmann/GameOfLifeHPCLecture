#!/bin/bash
#SBATCH --job-name=Perf_Test  # Job name
#SBATCH --ntasks=1            # Run on a single CPU
#SBATCH --time=00:10:00       # Time limit hrs:min:sec
#SBATCH --partition 16        # Partition 16

program=gameoflife
step_count=1000 #TODO
field_size_list=(100 200 300) #TODO
# Need to be same size (example will create 1x1 2x1 2x2 decomposition)
thread_decomposition_x=(1 2 2) #TODO
thread_decomposition_y=(1 1 2) #TODO

for field_size in ${field_size_list[@]}; do
  echo "$field_size"
  for i in ${!thread_decomposition_x[@]}; do
    thread_count=$((${thread_decomposition_x[$i]} * ${thread_decomposition_y[$i]}))
    OMP_NUM_THREADS=${thread_count} srun $program $field_size $field_size $step_count ${x_num[$i]} ${y_num[$i]}
  done
done