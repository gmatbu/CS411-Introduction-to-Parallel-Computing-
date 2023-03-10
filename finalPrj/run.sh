#!/bin/sh

rm *.out
rm *.csv

sbatch -N1 --time=1:00:00 get_runtimes.sh
sbatch -N1 --time=1:00:00 changing_dampening.sh
sbatch -N1 --time=1:00:00 changing_walk_size.sh
sbatch -N1 --time=1:00:00 changing_walk_size_and_threads.sh
