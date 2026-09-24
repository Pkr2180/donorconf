#!/bin/sh
cd "d:/pending work/bioinformatics - top q1/donorconf/experiments"
while ! grep -aq "written to" ../results/real/r1.log; do sleep 10; done
python real_crossdataset.py --name r3_oral --sets gingiva oropharynx oral_scc --scheme oral_coarse --n-ref 20 --n-cal 10 --n-tgt 10 --min-cells 40 --min-donors 4 --reps 60 > ../results/real/r3.log 2>&1
python real_crossdataset.py --name r3b_oral --sets gingiva oral_normal_multisite oral_scc --scheme oral_coarse --n-ref 20 --n-cal 10 --n-tgt 10 --min-cells 40 --min-donors 4 --reps 60 > ../results/real/r3b.log 2>&1
python real_crossdataset.py --name r2_blood --sets blood_A blood_B blood_C --scheme blood_fine --permute-roles --n-ref 30 --n-cal 30 --n-tgt 30 --min-cells 100 --reps 60 > ../results/real/r2.log 2>&1
echo QUEUE_DONE > ../results/real/queue.done
