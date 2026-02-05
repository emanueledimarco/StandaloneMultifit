Based on https://github.com/cms-eb-upgrade/vfe-toymc

Compile
====
./configure, then make


Create data
====

    #   #argv              1                   2                3         4       5           6                   7              8             9              10              11               12               13                  14
    #./CreateData   $temporal_shift   $number_of_events   $NSAMPLES   $NFREQ   $nPU   $signalAmplitude   $sigmaNoiseScale   $puFactor   $wf_name_string   $pu_shift   $noise-correlation   $pedestal   $slew-rate-distortion   $randomEnergy

    # examples
     ./CreateData        0                    11              16        6.25    200         10                  3                2       TestBeamPhase2        0              0.5              0.0              0                    1
    #./CreateData        0                    11              16        6.25    200         10                  1                1       TestBeamPhase2        0              0.5              0.5


    #./CreateData   0        11    10     25     0     10    0         0      CRRC43     0    0.5       0
    #./CreateData -13        11    10     25     0     10    0         0      CRRC43     0    0.5       0
    #./CreateData -13        1     10     25     0     10    0         0      CRRC43     0    0.5       0

Fit
====

    #   #argv        1              2                3             4            5
    #./multifit   inputFile   sigmaNoiseScale   fitPedestal   outputSuffix   maxEvents

    input="inputExternal/mysample_11_0.000_0.000_16_6.25_10.00_200.00_1.000_1.00_TestBeamPhase2_0.50_slew_0.00_FlatEnergy.root"

     ./multifit     $input          1                1

Plot
====

     root  plot/plotPulseInput.C\(\"$input\"\)

     output="output_-_noisescale1.00_fitPed1.root"

    #root plot/plot.C\(\"$output\"\)

    #                        fileName      treeName   nEvent
     root plot/plotPulse.C\(\"$output\",\"RecoAndSim\", 1\)

Example script
===
    source run.sh
