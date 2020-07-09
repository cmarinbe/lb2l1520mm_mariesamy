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
#gStyle.SetOptStat("emr")
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


    t.SetBranchStatus("*",0)
    t.SetBranchStatus("Phi",1)
    t.SetBranchStatus("cosTheta_proton",1)
    t.SetBranchStatus("cosTheta_lepton",1)

    t.SetBranchStatus("mum_PT",1)
    t.SetBranchStatus("mup_PT",1)
    t.SetBranchStatus("pplus_PT",1)
    t.SetBranchStatus("Kminus_PT",1)


    t.SetBranchStatus("Lb_M23",1)

    #
    #
    phi_acc = r.TH1F("phi_acc", "phi_acc",50, -3.16, 3.16)
    cosTheta_lepton_acc = r.TH1F("cosTheta_lepton_acc", "cosTheta_lepton_acc",20, -1., 1.)
    cosTheta_proton_acc = r.TH1F("cosTheta_proton_acc", "cosTheta_proton_acc",20, -1., 1.)
    qsq = r.TH1F("qsq","qsq", 50, 0, 17)
    #mass muon 106 MeV 0.011236 = 0.106*0.106 GeV/c2
    qsq_sel = r.TH1F("qsq_sel","qsq_sel", 50, 0.011236, 17)
    qsq_acc = r.TH1F("qsq_acc","qsq_acc", 50, 0.011236, 17)
    qsq_all = r.TH1F("qsq_all","qsq_all", 50, 0.011236, 17)


    lep_had_corr = r.TH2F("lep_had_corr", "lep_had_corr",20,-1,1, 20, -1,1)


    # cuts from RpK mu : 800 MeV, proton : 1000 MeV,  kaon 250 MeV

    cut_mum_PT = 0.400
    cut_mup_PT = 0.400
    cut_pplus_PT = 1.
    cut_Kminus_PT =0.300

    '''
    cut_mum_PT = 0
    cut_mup_PT = 0
    cut_pplus_PT = 0
    cut_Kminus_PT =0
    '''
    # binning definition in q2 [[0.1, 3], [3,6], [6, 8.86]] GeV2/c4
    qsq_bin_low_min = 0.1
    qsq_bin_low_max = 3

    qsq_bin_mid_min = 3
    qsq_bin_mid_max = 6

    qsq_bin_high_min = 6.
    qsq_bin_high_max = 8.86

    qsq_bin_central_min = 1.
    qsq_bin_central_max = 6

    qsq_bin_min = qsq_bin_central_min
    qsq_bin_max = qsq_bin_central_max

    for event in t  :
      if (t.Lb_M23*t.Lb_M23 > qsq_bin_min and  t.Lb_M23*t.Lb_M23 < qsq_bin_max  ) :
            qsq_all.Fill(t.Lb_M23*t.Lb_M23)
            if (t.mum_PT >  cut_mum_PT and  \
                t.mup_PT >  cut_mup_PT and  \
                t.pplus_PT > cut_pplus_PT and  \
                t.Kminus_PT > cut_Kminus_PT) :
                phi_acc.Fill(t.Phi)
                cosTheta_lepton_acc.Fill(t.cosTheta_lepton)
                cosTheta_proton_acc.Fill(t.cosTheta_proton)
                qsq_sel.Fill(t.Lb_M23*t.Lb_M23)
                lep_had_corr.Fill(t.cosTheta_lepton,t.cosTheta_proton )


    #eff_cosTheta_lepton = r.TF1("eff_cosTheta_lepton","[0]*(1 + [1]*x*x + [2]*x*x*x*x)")


    eff_cosTheta_lepton = r.TF1("eff_cosTheta_lepton"," [0]*(1. \
                                                       +[1]*x \
                                                       +[2]*(1./2)*(3*x*x -1)\
                                                       +[3]*(1./2)*(5*x*x*x -3*x)   \
                                                       +[4]*(1./8)*(35*x*x*x*x -30*x*x +3) \
                                                         )\
                                                        ")



    eff_cosTheta_proton = r.TF1("eff_cosTheta_proton","1. + [0]*x +[1]*x*x + [2]*x*x*x + [3]*x*x*x*x")


    #write the fit results in a root file
    results_file= ROOT.TFile.Open("fit_results.root", "update")

    #make a few plots
    lep_had_corr_canvas = r.TCanvas("lep_had_corr_canvas", "lep_had_corr_canvas", 400,400)
    lep_had_corr_canvas.cd()
    lep_had_corr.GetYaxis().SetTitle("cos #theta_{#Lambda}")
    lep_had_corr.GetXaxis().SetTitle("cos #theta_{#ell}")
    lep_had_corr.Draw("zcol")
    r.gPad.SaveAs("plots/lep_had_corr.png")


    qsq_all_canvas = r.TCanvas("qsq_all_canvas", "qsq_all_canvas", 400,400)
    qsq_all_canvas.cd()
    qsq_all.Draw("E")
    qsq_all.GetXaxis().SetTitle("q^{2} [GeV^{2}/c^{4}]")
    r.gPad.SaveAs("plots/qsq_all.png")

    qsq_sel_canvas = r.TCanvas("qsq_sel_canvas", "qsq_sel_canvas", 400,400)
    qsq_sel_canvas.cd()
    qsq_sel.Draw("E")
    qsq_sel.GetXaxis().SetTitle("q^{2} [GeV^{2}/c^{4}]")
    r.gPad.SaveAs("plots/qsq_sel.pdf")




    qsq_acc_canvas = r.TCanvas("qsq_acc_canvas", "qsq_acc_canvas", 400,400)
    qsq_acc_canvas.cd()
    qsq_acc.Divide(qsq_sel, qsq_all)
    qsq_acc.GetXaxis().SetTitle("q^{2} [GeV^{2}/c^{4}]")
    qsq_acc.Draw()
    r.gPad.SaveAs("plots/qsq_acc.png")



    phi_canvas = r.TCanvas("phi_canvas", "phi_canvas",400,400)
    phi_canvas.cd()
    phi_acc.Scale(1./phi_acc.Integral())
    phi_acc.Draw("E")
    phi_acc.GetXaxis().SetTitle("#phi")
    r.gPad.SaveAs("plots/phi.png")

    print('Fit parameters of cos Theta p/Lambda')

    proton_canvas = r.TCanvas("proton_canvas", "proton_canvas",400,400)
    proton_canvas.cd()
    cosTheta_proton_acc.Scale(1./cosTheta_proton_acc.Integral())
    cosTheta_proton_acc.GetYaxis().SetRangeUser(0,0.1)
    cosTheta_proton_acc.Draw("E")
    res_cosTheta = cosTheta_proton_acc.Fit("eff_cosTheta_proton","S")
    res_cosTheta.Write()
    eff_cosTheta_proton.SetLineColor(600)
    eff_cosTheta_proton.SetLineWidth(3)
    cosTheta_proton_acc.GetXaxis().SetTitle("cos #theta_{#Lambda}")

    r.gPad.SaveAs("plots/cosTheta_proton.png")


    print('Fit parameters of cos lepton')

    lepton_canvas = r.TCanvas("lepton_canvas", "lepton_canvas",400,400)
    lepton_canvas.cd()
    cosTheta_lepton_acc.Scale(1./cosTheta_lepton_acc.Integral())
    cosTheta_lepton_acc.GetYaxis().SetRangeUser(0,0.1)
    cosTheta_lepton_acc.Draw("E")
    res_cosTheta = cosTheta_lepton_acc.Fit("eff_cosTheta_lepton", "S")
    res_cosTheta.Write()
    eff_cosTheta_lepton.SetLineColor(2)
    eff_cosTheta_lepton.SetLineWidth(3)
    cosTheta_lepton_acc.GetXaxis().SetTitle("cos #theta_{l}")
    r.gPad.SaveAs("plots/cosTheta_lepton.png")
    r.gPad.SaveAs("plots/cosTheta_lepton.pdf")





    print ('-------------------------')
    print ('qsq_bin_min   {} GeV2'.format(qsq_bin_min))
    print ('qsq_bin_max   {} GeV2'.format( qsq_bin_max ))
    print ('-----')
    print ('cut_mum_PT    {} MeV'.format( cut_mum_PT))
    print ('cut_mup_PT    {} MeV'.format( cut_mup_PT))
    print ('cut_pplus_PT  {} MeV'.format( cut_pplus_PT ))
    print ('cut_Kminus_PT {} MeV'.format( cut_Kminus_PT))
    print ('-----')


    results_file.Close()
