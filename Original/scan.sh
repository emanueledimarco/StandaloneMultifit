NOISE=1

for ps in $(seq -6 6); do
  cat Pulse_template.h | sed "s/%PARAM/$ps/g" > Pulse.h
  rm *.x
  make
  rm inputExternal/mysample_100_0.000_0.000_16_6.25_100.00_200.00_${NOISE}.000_1.00_TestBeamPhase2_NoiseUncorrelated*.root
  for i in $(seq 0 10); do
    ./CreateData.x 0 100 16 6.25 200 100 ${NOISE} 1 TestBeamPhase2 0 0 0 0 0 $RANDOM &
  done
  wait
  hadd -f inputExternal/mysample_100_0.000_0.000_16_6.25_100.00_200.00_${NOISE}.000_1.00_TestBeamPhase2_NoiseUncorrelated_cat.root inputExternal/mysample_100_0.000_0.000_16_6.25_100.00_200.00_${NOISE}.000_1.00_TestBeamPhase2_NoiseUncorrelated_job*.root
  cp inputExternal/mysample_100_0.000_0.000_16_6.25_100.00_200.00_${NOISE}.000_1.00_TestBeamPhase2_NoiseUncorrelated_cat.root inputExternal/mysample_100_0.000_0.000_16_6.25_100.00_200.00_${NOISE}.000_1.00_TestBeamPhase2_NoiseUncorrelated_cat_globalshift_${ps}.root
  ./multifit.x inputExternal/mysample_100_0.000_0.000_16_6.25_100.00_200.00_${NOISE}.000_1.00_TestBeamPhase2_NoiseUncorrelated_cat.root $NOISE 1
  mv output_-_noisescale${NOISE}.00_fitPed1.root output/output_-_noisescale1.00_fitPed1_globalshift_${ps}.root
done
