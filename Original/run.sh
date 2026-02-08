#./CreateData.x 0 1000 16 6.25 200 10 1 1 TestBeamPhase2 0 0.5 0.5 0 1

input="inputExternal/mysample_1000_0.000_0.000_16_6.25_10.00_200.00_1.000_1.00_TestBeamPhase2_0.50_slew_0.00_FlatEnergy.root"
#input="inputExternal/mysample_11_0.000_0.000_16_6.25_10.00_200.00_1.000_1.00_TestBeamPhase2_0.50_slew_0.00_FlatEnergy.root"

./multifit.x $input 1 1

#output="output_-_noisescale1.00_fitPed1.root"

#root plot/plotPulse.C\(\"$output\",\"RecoAndSim\"\)

