#!/usr/bin/env python
# =============================================================================
# @file   simpleFit.py
# @author C. Marin Benito (carla.marin.benito@cern.ch)
# @date   22.12.2014
# =============================================================================
"""This script fits a simple pdf to a given dataset using RooFit."""

# imports
import os
import ROOT
from uncertainties import ufloat

RooFit         = ROOT.RooFit
RooRealVar     = ROOT.RooRealVar
RooArgList     = ROOT.RooArgList
RooArgSet      = ROOT.RooArgSet
RooDataSet     = ROOT.RooDataSet
RooGaussian    = ROOT.RooGaussian
RooExponential = ROOT.RooExponential
RooAddPdf      = ROOT.RooAddPdf

# definition of functions for this script
def simpleFit(tree, cuts, mean, xmin = 4000, xmax = 7000):
    """
    This function fits the "B_M" variable of a given TTree
    with a model formed by a Gaussian and an exponential pdf.
    All shape parameters are allowed to float in the fit. The 
    initial and range values are hardcoded in the code, except
    for the initial value of the Gaussian mean and the range
    of the B_M variable to be used.
    Returns the composed model (RooAbsPdf), the residual plot object
    and its chi2 value (float)
    Definition of the arguments:
    :tree: type TTree
    the root TTree that contains the variable to be fitted
    :cuts: type str
    optional cuts to apply to the TTree before fitting
    :mean: type float
    initial value for the Gaussian mean that will be floated
    during the fit
    :xmin: type float, optional
    minimum value of the B_M range to be fitted. Default: 4000
    :xmax: type float, optional
    maximum value of the B_M range to be fitted. Default: 7000
    """
    
    # define variables and pdfs
    B_M = RooRealVar("B_M","B_M", xmin, xmax)
    
    mean  = RooRealVar("mean", "mean",  mean, mean-50, mean+50)
    sigma = RooRealVar("sigma", "sigma", 80, 10, 150)
    gauss = RooGaussian("gauss", "gauss", B_M, mean, sigma)
    
    tau = RooRealVar("tau", "tau", -0.005, -0.01, 0.)
    exp = RooExponential("exp", "exp", B_M, tau)
    
    # define coefficiencts
    nsig = RooRealVar("nsig", "nsig", 1000, 0, 20000)
    nbkg = RooRealVar("nbkg", "nbkg", 1000, 0, 20000)
    
    # build model
    suma = RooArgList()
    coeff = RooArgList()
    
    suma.add(gauss)
    suma.add(exp)
    
    coeff.add(nsig)
    coeff.add(nbkg)
    
    model = ROOT.RooAddPdf("model", "model", suma, coeff)
    
    # define dataset
    if (cuts!=""): tree = tree.CopyTree(cuts)
    ds = RooDataSet("data", "dataset with x", tree, RooArgSet(B_M))
    
    # plot dataset and fit
    massFrame = B_M.frame()
    ds.plotOn(massFrame, Name="histo_data")
    
    fitResults = model.fitTo(ds)
    model.plotOn(massFrame, RooFit.VisualizeError(fitResults, 1),
                 RooFit.Name("curve_model"))
    model.plotOn(massFrame, RooFit.Components("gauss"), RooFit.LineColor(2),
                 RooFit.VisualizeError(fitResults, 1))
    model.plotOn(massFrame, RooFit.Components("exp")  , RooFit.LineColor(3),
                 RooFit.VisualizeError(fitResults, 1))
    model.paramOn(massFrame, Layout=(.55,.95,.93),
                  Parameters=RooArgSet(nsig, nbkg, mean, sigma, tau))
    ds.Draw()

    # print results
    print "{} has been fit to {} with a chi2 = {}".format(model.GetName(),
                                                          tree.GetName(), chi2)
    print "Total number of entries is: {}".format(ds.numEntries())
    print "Number of sig entries is: {:.0f} +- {:.0f}".format(nsig.getValV(),
                                                              nsig.getError())
    print "Number of bkg entries is: {:.0f} +- {:.0f}".format(nbkg.getValV(),
                                                              nbkg.getError())
    
    # compute S/B with error propagation from uncertainties module
    nsig = ufloat(nsig.getValV(), nsig.getError())
    nbkg = ufloat(nbkg.getValV(), nbkg.getError())
    signif = nsig/nbkg
    print "S/B = {:.2f} +- {:.2f}".format(signif.nominal_value, signif.std_dev)
    
    return ds, model


if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("file", action="store", type=str)
    parser.add_argument("-t", "--tree", default="DecayTree",
                        action="store", type=str)
    parser.add_argument("-m", "--mean", default=5280.,
                        action="store", type=float)
    parser.add_argument("-n", "--xmin", default=4000.,
                        action="store", type=float)
    parser.add_argument("-x", "--xmax", default=7000.,
                        action="store", type=float)
    parser.add_argument("-c", "--cuts", default="",
                        action="store", type=str)
    args = parser.parse_args()

    # sanity check
    if not os.path.exists(args.file):
        print "File doesn't exist! Exiting..."
        exit()

    # read data
    f = ROOT.TFile(args.file)
    t = f.Get(args.tree)

    ds, model, Plot, chi2 = simpleFit(t, args.cuts, args.mean,
                                      args.xmin, args.xmax)


#EOF

