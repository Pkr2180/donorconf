#!/bin/sh
cd "d:/pending work/bioinformatics - top q1/donorconf/experiments"
python real_b_adapt.py --name b_blood --sets blood_A blood_B blood_C --scheme blood_fine --permute-roles --n-ref 30 --n-cal 30 --n-tgt 30 --reps 30 > ../results/real/b_blood.log 2>&1
python real_b_adapt.py --name b_oral --sets gingiva oropharynx oral_scc --scheme oral_coarse --n-ref 20 --n-cal 10 --n-tgt 12 --min-cells 40 --min-donors 4 --reps 30 > ../results/real/b_oral.log 2>&1
python real_b_adapt.py --name b_oral_multisite --sets gingiva oral_normal_multisite oral_scc --scheme oral_coarse --n-ref 20 --n-cal 10 --n-tgt 12 --min-cells 40 --min-donors 4 --reps 30 > ../results/real/b_oral_multisite.log 2>&1
python real_benchmark.py --n-triples 24 --reps-per-triple 2 > ../results/real/c_benchmark.log 2>&1
echo DONE > ../results/real/queue2.done
