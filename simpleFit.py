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
def simpleFit(tree, cuts, xmean, xmin = 4000, xmax = 7000):
    """
    This function fits the "Lb_M" variable of a given TTree
    with a model formed by a Gaussian and an exponential pdf.
    All shape parameters are allowed to float in the fit. The 
    initial and range values are hardcoded in the code, except
    for the initial value of the Gaussian mean and the range
    of the Lb_M variable to be used.
    Returns the dataset and the composed model (RooAbsPdf)
    Definition of the arguments:
    :tree: type TTree
    the root TTree that contains the variable to be fitted
    :cuts: type str
    optional cuts to apply to the TTree before fitting
    :mean: type float
    initial value for the Gaussian mean that will be floated
    during the fit
    :xmin: type float, optional
    minimum value of the Lb_M range to be fitted. Default: 4000
    :xmax: type float, optional
    maximum value of the Lb_M range to be fitted. Default: 7000
    """
    
    # define variables and pdfs
    Lb_M = RooRealVar("Lb_M","Lb_M", xmin, xmax)
    
    mean  = RooRealVar("mean", "mean",  xmean, xmean-50, xmean+50)
    sigma = RooRealVar("sigma", "sigma", 20, 10, 50)
    gauss = RooGaussian("gauss", "gauss", Lb_M, mean, sigma)
    
    tau = RooRealVar("tau", "tau", -0.005, -0.01, 0.)
    exp = RooExponential("exp", "exp", Lb_M, tau)
    
    # define coefficiencts
    nsig = RooRealVar("nsig", "nsig", 100, 0, 2000)
    nbkg = RooRealVar("nbkg", "nbkg", 100, 0, 2000)
    
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
    ds = RooDataSet("data", "dataset with x", tree, RooArgSet(Lb_M))
    
    #  fit and save results
    fitResults = model.fitTo(ds, RooFit.Save(True))

    # plot dataset and fit results
    c = ROOT.TCanvas()
    massFrame = Lb_M.frame()
    ds.plotOn(massFrame, RooFit.Name("histo_data"))
    
    model.plotOn(massFrame)
    model.plotOn(massFrame, RooFit.Components("gauss"), RooFit.LineColor(2))
    model.plotOn(massFrame, RooFit.Components("exp")  , RooFit.LineColor(3))
    model.paramOn(massFrame, RooFit.Layout(.55,.95,.93))
    massFrame.Draw()
    c.SaveAs("fit.png")

    # print results
    print("{} has been fit to {}".format(model.GetName(),
                                         tree.GetName()))
    print("Total number of entries is: {}".format(ds.numEntries()))
    print("Number of sig entries is: {:.0f} +- {:.0f}".format(nsig.getValV(),
                                                              nsig.getError()))
    print("Number of bkg entries is: {:.0f} +- {:.0f}".format(nbkg.getValV(),
                                                              nbkg.getError()))
    
    # compute S/B with error propagation from uncertainties module
    nsigu = ufloat(nsig.getValV(), nsig.getError())
    nbkgu = ufloat(nbkg.getValV(), nbkg.getError())
    signif = nsigu/nbkgu
    print("S/B = {:.2f} +- {:.2f}".format(signif.nominal_value, signif.std_dev))
    
    return ds, model, c


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
        print("File doesn't exist! Exiting...")
        exit()

    # read data
    f = ROOT.TFile(args.file)
    t = f.Get(args.tree)

    ds, model, c = simpleFit(t, args.cuts, args.mean,
                             args.xmin, args.xmax)


#EOF

