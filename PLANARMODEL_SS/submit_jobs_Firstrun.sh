#!/usr/bin/sh

SMP=18
EPT=1
DOF=y
TTOT=1

for SWP in up down; do
    for FMP in 1 9 17 25; do
        for WFRC in `seq 250 270`; do
            # sbatch --time=02:00:25 --export=ALL,SMP=$SMP,EPT=$EPT,DOF=\'$DOF\',FMP=$FMP,WFRC=$WFRC,SWP=\'$SWP\',TTOT=$TTOT --job-name="$EPT$DOF-$FMP-$WFRC-$SWP" run_setup_FirstStudy.slurm
            if ! [ -f "./DATS/FIRST_1MAR21/RESU${SMP}_PT${EPT}${DOF}_F${FMP}000_W${WFRC}_Sweep${SWP}.mat" ]; then
                # echo "./DATS/FIRST_1MAR21/RESU${SMP}_PT${EPT}${DOF}_F${FMP}000_W${WFRC}_Sweep${SWP}.mat Exists"
            # else
                sbatch --time=02:00:25 --export=ALL,SMP=$SMP,EPT=$EPT,DOF=\'$DOF\',FMP=$FMP,WFRC=$WFRC,SWP=\'$SWP\',TTOT=$TTOT --job-name="$EPT$DOF-$FMP-$WFRC-$SWP" run_setup_FirstStudy.slurm
                echo "./DATS/FIRST_1MAR21/RESU${SMP}_PT${EPT}${DOF}_F${FMP}000_W${WFRC}_Sweep${SWP}.mat Doesn't Exist"
            fi
        done
    done
done

# echo $SMP $EPT $DOF $FMP $WFRC

# sbatch --time=00:00:30 --export=ALL,SMP=$SMP,EPT=$EPT,DOF=$DOF,FMP=$FMP,WFRC=$WFRC --job-name="$EPT$DOF-$FMP$WFRC" run_setup.slurm
            
