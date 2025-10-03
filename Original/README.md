Based on https://github.com/cms-eb-upgrade/vfe-toymc

Create data
====

    g++ -o CreateData  CreateData.cc -std=c++11 `root-config --cflags --glibs`

    ./CreateData     
                        1                     2                  3            4        5           6                   7                 8           9              10              11               12              13
    ./CreateData   temporal_shift      number_of_events       NSAMPLES       NFREQ     nPU    signalAmplitude    sigmaNoiseScale     puFactor  wf_name_string    pu_shift     noise-correlation    pedestal  slew-rate-distortion

    ./CreateData        0                   11                  16           6.25       0           10                 1                 1      IdealPSphase2        0             0.5                0

    ./CreateData -13        11    10     25     0     10    0         0      CRRC43     0    0.5       0
    ./CreateData -13        1     10     25     0     10    0         0      CRRC43     0    0.5       0
    
Fit
====

    g++ -o simple.multifitsimple.multifit.cc PulseChiSqSNNLS.cc -std=c++11 `root-config --cflags --glibs`

                   input     output    NSAMPLES   NFREQ    time-shift (13 has to go for 0)      pedestal-shift
    ./simple.multifit  input/mysample_11_-13.000_0.000_10_25.00_10.00_0.00_0.000_0.00_CRRC43_0.00.root       output/mysample_11_-13.000_0.000_10_25.00_10.00_0.00_0.000_0.00_CRRC43_0.00.root         10            25        13        0.0

    ./simple.multifit  inputExternal/mysample_1_-13.000_0.000_10_25.00_10.00_0.00_0.000_0.00_CRRC43_0.00.root       output/mysample_1_-13.000_0.000_10_25.00_10.00_0.00_0.000_0.00_CRRC43_0.00.root   10            25        13        0.0

    
    
Plot
====

    r99t  plot/plotPulseInput.C\(\"input/mysample_11_-13.000_0.000_10_25.00_10.00_0.00_0.000_0.00_CRRC43_0.00.root\"\)
    r99t  plot/plotPulseInput.C\(\"input/mysample_100_-13.000_0.000_10_25.00_10.00_80.00_0.000_0.00_CRRC43_-1.00.root\"\)
    
    r99t  plot/plot.C\(\"output/mysample_11_-13.000_0.000_10_25.00_10.00_0.00_0.000_0.00_CRRC43_0.00.root\"\)
    r99t plot/plotPulse.C\(\"output/mysample_11_-13.000_0.000_10_25.00_10.00_0.00_0.000_0.00_CRRC43_0.00.root\",2\)
    r99t plot/plotPulse.C\(\"output/mysample_11_-13.000_0.000_10_25.00_10.00_0.00_0.000_0.00_CRRC43_0.00.root\",1\)
