#input="inputExternal/mysample_1000_0.000_0.000_16_6.25_10.00_200.00_1.000_1.00_TestBeamPhase2_0.50_slew_0.00_FlatEnergy.root"
#input="inputExternal/mysample_11_0.000_0.000_16_6.25_10.00_200.00_1.000_1.00_TestBeamPhase2_0.50_slew_0.00_FlatEnergy.root"


#input="inputExternal/mysample_10_0.000_0.000_16_6.25_200.00_200.00_1.000_1.00_TestBeamPhase2_0.50_slew_0.00_FlatEnergy_cat.root"

#input="inputExternal/mysample_10_0.000_0.000_16_6.25_10.00_200.00_1.000_1.00_TestBeamPhase2_0.00_slew_0.00_FlatEnergy.root"
#input="inputExternal/mysample_10_0.000_0.000_16_6.25_10.00_200.00_1.000_1.00_TestBeamPhase2_NoiseUncorrelated_FlatEnergy.root"

#input="inputExternal/mysample_10_0.000_0.000_16_6.25_10.00_200.00_1.000_1.00_TestBeamPhase2_NoiseUncorrelated_FlatEnergy.root"

#input="inputExternal/mysample_10_0.000_0.000_16_6.25_10.00_200.00_1.000_1.00_TestBeamPhase2_NoiseUncorrelated_FlatEnergy.root"

#input="inputExternal/mysample_10_0.000_0.000_16_6.25_200.00_200.00_1.000_1.00_TestBeamPhase2_NoiseUncorrelated_FlatEnergy.root"
#input="inputExternal/mysample_1000_0.000_0.000_16_6.25_200.00_200.00_1.000_1.00_TestBeamPhase2_NoiseUncorrelated_FlatEnergy_cat.root"

#input="inputExternal/mysample_100_4.000_0.000_16_6.25_100.00_200.00_1.000_1.00_TestBeamPhase2_NoiseUncorrelated.root"

#input="inputExternal/mysample_100_0.000_0.000_16_6.25_100.00_200.00_1.000_1.00_TestBeamPhase2_NoiseUncorrelated.root"

#input="inputExternal/mysample_100_0.000_0.000_16_6.25_100.00_200.00_1.000_1.00_TestBeamPhase2_NoiseUncorrelated_cat.root"

input=$1
./multifit.x $input 1 1

#output="output_-_noisescale1.00_fitPed1.root"

#root plot/plotPulse.C\(\"$output\",\"RecoAndSim\"\)
