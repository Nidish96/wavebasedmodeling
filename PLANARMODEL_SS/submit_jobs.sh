#!/usr/bin/sh

SMP=17

for EPT in `seq 1 8`; do
    for DOF in \'x\' \'y\'; do
        for FMP in 1 9 17 25; do
            for WFRC in `seq 250 270`; do
                echo $SMP $EPT $DOF $FMP $WFRC
                sbatch --export=ALL,SMP=$SMP,EPT=$EPT,DOF=$DOF,FMP=$FMP,WFRC=$WFRC --job-name="$EPT$DOF-$FMP-$WFRC" run_setup.slurm
            done
        done
    done
done

# echo $SMP $EPT $DOF $FMP $WFRC

# sbatch --time=00:00:30 --export=ALL,SMP=$SMP,EPT=$EPT,DOF=$DOF,FMP=$FMP,WFRC=$WFRC --job-name="$EPT$DOF-$FMP$WFRC" run_setup.slurm
            
