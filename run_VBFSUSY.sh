#!/bin/bash
#set -e # exit when any command fails

deltaeta=3.0
nevents=50000

# Higgsino WinoBino
for params in Higgsino; do
    # 13 100
    for ecms in 100; do
	# 150 250 2000
	for mass in 250; do
	    min_mjj=0.5
            mmjj_step=2
            max_mmjj_TeV=4.5
            if [[ $ecms == 100 ]]; then
		# don't slice for now
		nevents=200000
		min_mjj=1
		mmjj_step=1
		max_mmjj_TeV=1
            fi
            for i_mmjj in $(seq ${min_mjj} ${mmjj_step} ${max_mmjj_TeV}); do

		mmjj=$(bc <<< "scale=0; ${i_mmjj}*1000/1")
		mmjj_max=$(bc <<< "scale=0; (${mmjj}+${mmjj_step}*1000)/1")
		if [[ $i_mmjj == $max_mmjj_TeV ]]; then
                    mmjj_max="-1"
		fi

		seed=0
		./run_VBFSUSY_standalone.sh \
		    -E ${ecms} \
		    -M ${mass} \
		    -P VBFSUSY_EWKQCD \
		    -c 4 \
		    -p ${params} \
		    -N ${nevents} \
		    -m ${mmjj} \
		    -x ${mmjj_max} \
		    -e ${deltaeta} \
		    -d ${seed}

	    done
	done
    done
done
