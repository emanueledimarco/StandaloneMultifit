#for i in $(seq 0 10); do bash -c "./CreateData.x -100 5000 16 6.25 200 200 1 1 TestBeamPhase2 0 0 0 0 1 $RANDOM&"; done

#./CreateData.x -1 10 16 6.25 0 100 1 1 TestBeamPhase2 0 0 0 0 0

for i in $(seq 0 10); do bash -c "./CreateData.x 0 100 16 6.25 200 100 1 1 TestBeamPhase2 0 0 0 0 0 $RANDOM&"; done
