'''
Feb 2020
Let's use a python script to make Carla happy
'''
import ROOT as r
from math import *
from ROOT import gStyle, gROOT



import ROOT
import argparse
import getopt
import numpy as np
from math import acos
from array import array
from copy import deepcopy

from ROOT import RooFit,RooLegendre,RooAbsPdf, RooPolynomial, RooArgList,RooRealVar, RooGenericPdf


RooRealVar =  ROOT.RooRealVar
gStyle.SetOptStat("")
gROOT.ProcessLine(".x ../utils/lhcbStyle.C")


ROOT.gSystem.Load("libMathMore") #to get Legendre

from math import pow


if __name__ == '__main__':
    # handle arguments
    parse = argparse.ArgumentParser()
    parse.add_argument('fname', type=str, help    ='File to add the variable')
    parse.add_argument('--parent', type=str, default='Lb', help='parent particle')
    parse.add_argument('--child1' , type=str, default='pplus', help='child particle 1')
    parse.add_argument('--child2' , type=str, default='Kminus', help='child particle 2')
    parse.add_argument('--child3' , type=str, default='mup', help='child particle 3')
    parse.add_argument('--child4' , type=str, default='mum', help='child11 particle 4')
    parse.add_argument('-t', '--tname', default='DecayTree', type=str, help='Tree name')
    parse.add_argument('-g', '--generated', default=False, action='store_true', help='Use TRUE values?') # if set computes values for TRUE and reco, else only reco
    parse.add_argument('-b', '--bmassfit', default=False, action='store_true', help='Assume Lb_M=5620?') # if set computes values for BMassFit, else standard reco
    args = parse.parse_args()

    # File to be updated
    f = ROOT.TFile(args.fname, 'UPDATE')
    t = f.Get(args.tname)


    cut_mum_PT = 0.400
    cut_mup_PT = 0.400
    cut_pplus_PT = 1.
    cut_Kminus_PT =0.250

    # binning definition in q2 [[0.1, 3], [3,6], [6, 8.86]] GeV2/c4
    qsq_bins ={
      "bin1" : [0.1, 3],
      "bin2" : [3.0, 6],
      "bin3" : [6.0, 8.86],
      "bin4" : [1.0, 6.0],
    }

    #coefficients of the Legendre Polynomials that we need to fit
    cLep0=RooRealVar ("cLep0","cLep0",0,-1,1)
    cLep1=RooRealVar ("cLep1","cLep1",0,-1,1)
    cLep2=RooRealVar ("cLep2","cLep2",0,-1,1)
    cLep3=RooRealVar ("cLep3","cLep3",0,-1,1)
    cLep4=RooRealVar ("cLep4","cLep4",0,-1,1)
    cLep5=RooRealVar ("cLep5","cLep5",0,-1,1)
    cLep6=RooRealVar ("cLep6","cLep6",0,-1,1)


    cLambda0=RooRealVar ("cLambda0","cLambda0",0,-4,4)
    cLambda1=RooRealVar ("cLambda1","cLambda1",0,-4,4)
    cLambda2=RooRealVar ("cLambda2","cLambda2",0,-4,4)
    cLambda3=RooRealVar ("cLambda3","cLambda3",0,-4,4)
    cLambda4=RooRealVar ("cLambda4","cLambda4",0,-4,4)
    cLambda5=RooRealVar ("cLambda5","cLambda5",0,-4,4)

    Lb_M23                = ROOT.RooRealVar("Lb_M23","Lb_M23", 0, 5)
    mum_PT                = ROOT.RooRealVar("mum_PT","mum_PT", 0, 100 )
    mup_PT                = ROOT.RooRealVar("mup_PT","mup_PT", 0, 100 )
    pplus_PT              = ROOT.RooRealVar("pplus_PT","pplus_PT", 0, 100 )
    Kminus_PT             = ROOT.RooRealVar("Kminus_PT"," Kminus_PT", 0, 100 )
    cosTheta_lepton       = ROOT.RooRealVar("cosTheta_lepton","cosTheta_lepton", -1, 1 )
    cosTheta_Lst       = ROOT.RooRealVar("cosTheta_Lst","cosTheta_Lst", -1, 1 )

    data_all    = ROOT.RooDataSet( "data_all","data for cosThetaL_distribution",
                  ROOT.RooArgSet(Lb_M23, mum_PT,mup_PT, pplus_PT, Kminus_PT, cosTheta_lepton, cosTheta_Lst ), ROOT.RooFit.Import(t))



    Legendre_Lep_all = "1. \
                                +cLep1*cosTheta_lepton  \
                                +cLep2*(1./2)*(3*cosTheta_lepton*cosTheta_lepton -1)\
                                +cLep3*(1./2)*(5*cosTheta_lepton*cosTheta_lepton*cosTheta_lepton -3*cosTheta_lepton)\
                                +cLep4*(1./8)*(35*cosTheta_lepton*cosTheta_lepton*cosTheta_lepton*cosTheta_lepton -30*cosTheta_lepton*cosTheta_lepton+3)\
                                +cLep5*(1./8)*(65*cosTheta_lepton*cosTheta_lepton*cosTheta_lepton*cosTheta_lepton*cosTheta_lepton \
                                             - 70*cosTheta_lepton*cosTheta_lepton*cosTheta_lepton +15*cosTheta_lepton )\
                                +cLep6*(1./16)*(231*cosTheta_lepton*cosTheta_lepton*cosTheta_lepton*cosTheta_lepton*cosTheta_lepton*cosTheta_lepton\
                                              - 315*cosTheta_lepton*cosTheta_lepton*cosTheta_lepton*cosTheta_lepton + 105*cosTheta_lepton*cosTheta_lepton -5)"


    Legendre_Lep_even = "1. \
                                +cLep2*(1./2)*(3*cosTheta_lepton*cosTheta_lepton -1)\
                                +cLep4*(1./8)*(35*cosTheta_lepton*cosTheta_lepton*cosTheta_lepton*cosTheta_lepton -30*cosTheta_lepton*cosTheta_lepton+3)\
                                +cLep6*(1./16)*(231*cosTheta_lepton*cosTheta_lepton*cosTheta_lepton*cosTheta_lepton*cosTheta_lepton*cosTheta_lepton\
                                - 315*cosTheta_lepton*cosTheta_lepton*cosTheta_lepton*cosTheta_lepton + 105*cosTheta_lepton*cosTheta_lepton -5)"






    Legendre_Lep_o2_o4 = "1. \
                               +cLep2*(1./2)*(3*cosTheta_lepton*cosTheta_lepton -1)\
                               +cLep4*(1./8)*(35*cosTheta_lepton*cosTheta_lepton*cosTheta_lepton*cosTheta_lepton -30*cosTheta_lepton*cosTheta_lepton+3)"





    Legendre_Lambda_all =  "1. \
                                +cLambda1*cosTheta_Lst  \
                                +cLambda2*(1./2)*(3*cosTheta_Lst*cosTheta_Lst -1)\
                                +cLambda3*(1./2)*(5*cosTheta_Lst*cosTheta_Lst*cosTheta_Lst -3*cosTheta_Lst)\
                                +cLambda4*(1./8)*(35*cosTheta_Lst*cosTheta_Lst*cosTheta_Lst*cosTheta_Lst -30*cosTheta_Lst*cosTheta_Lst+3)\
                                +cLambda5*(1./8)*(65*cosTheta_Lst*cosTheta_Lst*cosTheta_Lst*cosTheta_Lst*cosTheta_Lst\
                                               - 70*cosTheta_Lst*cosTheta_Lst*cosTheta_Lst +15*cosTheta_Lst )"




    Legendre_Lambda_1_2 =  "1. \
                            +cLambda1*cosTheta_Lst  \
                            +cLambda2*(1./2)*(3*cosTheta_Lst*cosTheta_Lst -1)"



    results = open("LegendreCoefficients.txt","a")
    results.write("Acceptance fit results : \n")
    for bin in qsq_bins.items(): #loop over the bins
        print ('----------------------------------------------')
        print ('----------------------------------------------')
        print ('----------------------------------------------')
        print ('---- qsq_bin_min {}  GeV2'.format (bin[1][0]))
        print ('---- qsq_bin_max {}  GeV2'.format (bin[1][1]))

        qsq_bin_min =  (bin[1][0])
        qsq_bin_max =  (bin[1][1])
        name_of_bin = bin[0]

        # ----------------------------------------------
        # define the datsets for fitting
        # ----------------------------------------------
        bin_cut = "((Lb_M23*Lb_M23 > {}) && (Lb_M23*Lb_M23 < {}))".format(qsq_bin_min, qsq_bin_max)
        pt_cut  = "((mum_PT > {}) && (mup_PT > {}) && (pplus_PT > {}) && (Kminus_PT > {}))".format(cut_mum_PT, cut_mup_PT, cut_pplus_PT, cut_Kminus_PT)
        tot_cut = "{} && {}".format(bin_cut, pt_cut)

        data_selected = data_all.reduce(tot_cut)
        data_all.Print()
        data_selected.Print()
        data_selected.get(49).Print("cosTheta_lepton")

        # ----------------------------------------------
        # define the PDFs
        # ----------------------------------------------

        Legendre_Lep = RooGenericPdf ("Legendre_Lep", Legendre_Lep_even,\
                       RooArgList (cosTheta_lepton, cLep2, cLep4,cLep6))

        Legendre_Lambda = RooGenericPdf ("Legendre_Lambda", Legendre_Lambda_all, \
                          RooArgList (cosTheta_Lst,cLambda1, cLambda2,cLambda3, cLambda4, cLambda5))
        # ----------------------------------------------
        # Make the fit
        # ----------------------------------------------

        fit_results_RooFit_Lep = Legendre_Lep.fitTo(data_selected, ROOT.RooFit.Save())
        fit_results_RooFit_Lambda = Legendre_Lambda.fitTo(data_selected, ROOT.RooFit.Save())
        #write the the numbers that we need in a text file for the toys later
        results.write("----\n")
        results.write("----\n")
        results.write(name_of_bin+" Lepton coefficients values ord2 : {} ord4 : {} ord6 :{} \n".format(cLep2.getVal(),cLep4.getVal(),cLep6.getVal()))
        results.write(name_of_bin+" Lepton coefficients errors ord2 :{}  ord4 :{}  ord6 :{} \n".format(cLep2.getError(),cLep4.getError(),cLep6.getError()))
        results.write("----\n")
        
        # ----------------------------------------------
        #plotting
        # ----------------------------------------------


        cLep =  ROOT.TCanvas ("DCanvas3Fe","Fit",750,750)
        pad1  = ROOT.TPad("pad1", "The pad 70% of the height",0.0,0.25,1.0,1.0)
        pad2  = ROOT.TPad("pad2", "The pad 30% of the height",0.0,0.0,1.0,0.25)
        pad1.SetBottomMargin(0.065)
        pad1.SetBorderMode(0)
        pad2.SetTopMargin(0.00001)
        pad2.SetBottomMargin(0.2999)
        pad2.SetBorderMode(0)
        pad2.Draw()
        pad1.Draw()
        pad1.cd()
        frameLep = cosTheta_lepton.frame()
        data_selected.plotOn(frameLep,RooFit.Binning(25))
        Legendre_Lep.plotOn(frameLep)
        frameLep.GetYaxis().SetTitle('A.U')
        frameLep.GetXaxis().SetTitle('cos #theta_{l}')
        frameLep.Draw()
        pad2.cd()
        hresid = frameLep.pullHist()
        frameLep_res  = cosTheta_lepton.frame()
        frameLep_res.addPlotable(hresid,"P")
        frameLep_res.GetXaxis().SetTitle('cos #theta_{l}')
        frameLep_res.Draw()
        cLep.SaveAs("plots/cosThetaLepton_acceptance_pulls_{}.png".format(name_of_bin))

        cLep_final =  ROOT.TCanvas ("DCanvas3Fe_final","Fit",400,400)
        cLep_final.cd()
        frameLep.Draw()
        cLep_final.SaveAs("plots/cosThetaLepton_acceptance_{}.png".format(name_of_bin))
        # -----

        cLambda = ROOT.TCanvas("cLambda", "cLambda", 750,750)
        pad1L  = ROOT.TPad("pad1L", "The pad 70% of the height",0.0,0.25,1.0,1.0)
        pad2L  = ROOT.TPad("pad2L", "The pad 30% of the height",0.0,0.0,1.0,0.25)
        pad1L.SetBottomMargin(0.065)
        pad1L.SetBorderMode(0)
        pad2L.SetTopMargin(0.00001)
        pad2L.SetBottomMargin(0.2999)
        pad2L.SetBorderMode(0)
        pad2L.Draw()
        pad1L.Draw()

        pad1L.cd()
        frameLambda = cosTheta_Lst.frame()
        data_selected.plotOn(frameLambda,RooFit.Binning(25))
        Legendre_Lambda.plotOn(frameLambda)
        frameLambda.GetYaxis().SetTitle('A.U')
        frameLambda.GetXaxis().SetTitle('cos #theta_{#Lambda}')
        frameLambda.Draw()

        pad2L.cd()
        hresid_Lambda = frameLambda.pullHist()
        frameLambda_res  = cosTheta_Lst.frame()
        frameLambda_res.addPlotable(hresid_Lambda,"P")
        frameLambda_res.GetXaxis().SetTitle('cos #theta_{#Lambda}')
        frameLambda_res.Draw()
        cLambda.SaveAs("plots/cosThetaLambda_acceptance_pulls_{}.png".format(name_of_bin))

        cLambda_final =  ROOT.TCanvas ("DCanvas3Fe_final","Fit",400,400)
        cLambda_final.cd()
        frameLambda.Draw()
        cLambda_final.SaveAs("plots/cosThetaLambda_acceptance_{}.png".format(name_of_bin))

    print ('-------------------------')
    print ('-----')
    print ('cut_mum_PT    {} GeV'.format( cut_mum_PT))
    print ('cut_mup_PT    {} GeV'.format( cut_mup_PT))
    print ('cut_pplus_PT  {} GeV'.format( cut_pplus_PT ))
    print ('cut_Kminus_PT {} GeV'.format( cut_Kminus_PT))
    print ('-----')
