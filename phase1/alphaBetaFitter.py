import ROOT as rt
import math

class AlphaBetaFitter:

    def __init__(self,
                 doEB=True,
                 pedestal=0.,
                 sampling=25):

        self.doEB = doEB

        alpha = 1.250 if doEB else 1.283
        beta  = 1.600 if doEB else 1.674

        beta *= sampling

        self.fcn = rt.TF1("alphabeta",
            """[4] + [0] * (
              ( (x - [3]) > -([1]*[2]) )
              ?
              (
                  TMath::Power(
                      1 + (x - [3])/([1]*[2]),
                      [1]
                  )
                  *
                  TMath::Exp(
                      -(x - [3])/[3]
                  )
              )
              : 0
          )""",
            0,
            12 * sampling,
        )
        self.fcn.SetParNames(
            "norm",
            "alpha",
            "beta",
            "tmax",
            "pedestal"
        )

        self.fcn.FixParameter(0,1.)

        self.fcn.SetParameter(1,alpha)
        self.fcn.SetParameter(2,beta)
        self.fcn.SetParameter(3,5.5*sampling)

        if doEB:
            self.fcn.SetParLimits(1,0.8,2.5)
        else:
            self.fcn.SetParLimits(1,0.5,2.5)

        self.fcn.SetParLimits(2,0.8*sampling,3*sampling)
        self.fcn.SetParLimits(3,1.5*sampling,7*sampling)

        self.fcn.FixParameter(4,pedestal)

    def fit(self,graph):

        graph.Fit(
            self.fcn,
            "Q0"
        )

        pars = [
            self.fcn.GetParameter(i)
            for i in range(5)
        ]

        errs = [
            self.fcn.GetParError(i)
            for i in range(5)
        ]

        return pars, errs

    def getFcn(self):

        return self.fcn

