  #!/usr/bin/env python
# =============================================================================
# @file   add_angles.py
# @author C. Marin Benito (carla.marin.benito@cern.ch)
# @author Y.Amhis  (yasmine.sara.amhis@cern.ch)
# @date   25.02.19
# Modification : 13.12.2019 to adapt to pKmm
# =============================================================================
"""Add helicity angles to ntuple.
   Computation based on TupleToolHelicity by Olivier"""
# =============================================================================
# By setting appropriate flags, computes everything for standard variables
# as well as TRUE variables in the MC, or assuming Lb=5620.

import ROOT
import argparse
import getopt
import numpy as np
from math import acos
from array import array
from copy import deepcopy


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

    # Variables and branches that will be added
    variables = [   'Lb_M', #0
                    'Lb_M01', #1 pK invariant mass
                    'Phi', #2
                    'cosTheta_dilepton',#3
                    'cosTheta_Resonance', #4
                    'cosTheta_Lst',#5  #this one is called theta_Lambda in arXiv:1903.00448
                    'p_P_in_pK', #6
                    'dl_P_in_Lb',#7
                    'cosTheta_lepton', #8
                    'Lb_M23',#9 dilepton invariant mass
                    'lepton_P_in_dilepton'] #10



    branches = ['%s' %v for v in variables]
    br_names = branches
    # add branches
    br_var = [array('f',[0.]) for br in br_names]
    br_bra = [t.Branch(br, var, br+'/F')
              for br, var in zip(br_names, br_var)]
    print('adding branches: ', br_names)
    # do the work
    v_p = ROOT.TLorentzVector() #parent
    v_resonance = ROOT.TLorentzVector() #resonance
    v_c1 = ROOT.TLorentzVector() #child1
    v_c2 = ROOT.TLorentzVector() #child2
    v_c3 = ROOT.TLorentzVector() #child3
    v_c4 = ROOT.TLorentzVector() #child4
    v_dilepton = ROOT.TLorentzVector() #dilepton


    PX, PY, PZ, PE = ('_PX', '_PY', '_PZ', '_E') # in RapidSim E is just E not PE like in dtt
    # Define prefixes in case BMassFit variables should be used. for pplus, kminus prefix2 is needed as well.
    PREFIX1 = ''
    PREFIX2 = ''
    n_entries = t.GetEntries()
    n_print = int(n_entries/10)+1
    for i in range(n_entries):
    #for i in range(1):

        t.GetEntry(i)
        if (i%n_print == 0): print ('Entry {} - {}%'.format(i,int((i/n_entries)*100)))
        # Get parent, vector and child momentum
        v_c1.SetPxPyPzE(getattr(t, PREFIX1+PREFIX2+args.child1+PX),
                       getattr(t, PREFIX1+PREFIX2+args.child1+PY),
                       getattr(t, PREFIX1+PREFIX2+args.child1+PZ),
                       getattr(t, PREFIX1+PREFIX2+args.child1+PE))
        v_c2.SetPxPyPzE(getattr(t, PREFIX1+PREFIX2+args.child2+PX),
                        getattr(t, PREFIX1+PREFIX2+args.child2+PY),
                        getattr(t, PREFIX1+PREFIX2+args.child2+PZ),
                        getattr(t, PREFIX1+PREFIX2+args.child2+PE))
        v_c3.SetPxPyPzE(getattr(t, PREFIX1+PREFIX2+args.child3+PX),
                        getattr(t, PREFIX1+PREFIX2+args.child3+PY),
                        getattr(t, PREFIX1+PREFIX2+args.child3+PZ),
                        getattr(t, PREFIX1+PREFIX2+args.child3+PE))
        v_c4.SetPxPyPzE(getattr(t, PREFIX1+PREFIX2+args.child4+PX),
                        getattr(t, PREFIX1+PREFIX2+args.child4+PY),
                        getattr(t, PREFIX1+PREFIX2+args.child4+PZ),
                        getattr(t, PREFIX1+PREFIX2+args.child4+PE))

        v_p.SetPxPyPzE(getattr(t, args.parent+PX),
                           getattr(t, args.parent+PY),
                           getattr(t, args.parent+PZ),
                           getattr(t, args.parent+PE))
        v_resonance  = v_c1 + v_c2   # A+B
        v_dilepton   = v_c3 + v_c4   # C+D

        if (i==0):
            print("## Original momenta for parent, resonance, dilepton, c1, c2, c3, c4:")
            for v in [v_p, v_resonance, v_dilepton, v_c1, v_c2, v_c3, v_c4]:
                v.Print()


        #----------------------------------------------
        #    get boosts to parent rest frames
        #    ----------------------------------------------
        boost_p = -v_p.BoostVector() # boost from lab to parent RF

        v_dilepton.Boost(boost_p) # boost dilepton to parent RF
        u_dl = v_dilepton.Vect().Unit() # get unit vector

        v_resonance.Boost(boost_p) # boost resonance to parent RF
        u_r = v_resonance.Vect().Unit() # get unit vector

        boost_rp = -v_resonance.BoostVector() # boost vector to res RF in Lb RF
        v_c1.Boost(boost_p) # boost child to parent RF
        v_c1.Boost(boost_rp) # boost to resonance RF in Lb RF
        u_c1 = v_c1.Vect().Unit() # get unit vector
        # for cross-check only
        v_c2.Boost(boost_p) # boost child to parent RF
        v_c2.Boost(boost_rp) # boost to resonance RF in Lb RF

        boost_dlp = -v_dilepton.BoostVector() # boost vector to dilepton RF in Lb RF
        v_c3.Boost(boost_p) # boost child to parent RF
        v_c3.Boost(boost_dlp) # boost to dilepton RF in Lb RF
        u_c3 = v_c3.Vect().Unit() # get unit vector
        # for cross-check only
        v_c4.Boost(boost_p) # boost child to parent RF
        v_c4.Boost(boost_dlp) # boost to dilepton RF in Lb RF


        u_p = v_p.Vect().Unit() # get unit vector // to parent direction
        u_resonance = v_resonance.Vect().Unit() # get unit vector // to resonance direction


        #Get normal vectors
        n_p = u_p.Cross(u_resonance) # normal to parent decay plane
        n_v = u_resonance.Cross(u_c1) # normal to resonance decay plane

        # ----------------------------------------------------
        #    compute cosTheta_dl:
        #    angle between dilepton momentum and quantization axis,
        #    defined parallel to z in EvtGen, in Lb rest frame
        #    ALSO: get dilepton momentum in Lb frame
        # ----------------------------------------------------
        if args.generated:
            z0_TRUE  = ROOT.TLorentzVector(0., 0., 1., 0.) # z axis in Lab Frame
            z0_TRUE.Boost(boost_p_TRUE) # boost to parent RF
            u_n_TRUE = z0_TRUE.Vect().Unit() # get unit vector
            v_dilepton_TRUE.Boost(boost_p_TRUE) # boost photon to parent RF
            u_dl_TRUE = v_dilepton_TRUE.Vect().Unit() # get unit vector
            br_var[len(variables)+3][0] = u_dl_TRUE.Dot(u_n_TRUE) # get cos angle between the two CHECK HERE !!!!
            br_var[len(variables)+7][0] = v_dilepton_TRUE.P() # get dilepton momentum in Lb frame

        z0  = ROOT.TLorentzVector(0., 0., 1., 0.) # z axis in Lab Frame
        z0.Boost(boost_p) # boost to parent RF
        u_n = z0.Vect().Unit() # get unit vector
        u_dl = v_dilepton.Vect().Unit() # get unit vector
        br_var[3][0] = u_dl.Dot(u_n) # get cos angle between the two
        br_var[7][0] = v_dilepton.P() # get dilepton  momentum in Lb frame CHECK HERE !!!!


        #-------------------------------------------
        # same for cosTheta_L ie : The Lambda*
        #-------------------------------------------

        if args.generated:
            v_resonance_TRUE.Boost(boost_p_TRUE) # boost resonance to parent RF
            u_r_TRUE = v_resonance_TRUE.Vect().Unit() # get unit vector

        br_var[4][0] = u_r.Dot(u_n) # get cos angle between the two
        if (i == 0):
            if args.generated:
                print ("Polarization axis and p_dl and p_r in Lb rest frame TRUE:")
                u_n_TRUE.Print(); u_dl_TRUE.Print(); u_r_TRUE.Print()
                print ("cosTheta_dl and cosTheta_Lst TRUE (should add up to pi)")
                print (acos(br_var[len(variables)+3][0]), acos(br_var[len(variables)+4][0]))
            print ("Polarization axis and p_dilepton and p_Lst in Lb rest frame:")
            u_n.Print(); u_dl.Print(); u_r.Print()
            print ("cosTheta_dilepton and cosTheta_Lst (should add up to pi)")
            print (acos(br_var[3][0]), acos(br_var[4][0]))


        # -----------------------------------------------
        #     compute cosTheta_p:
        #     angle between p momentum in L(X) rest frame and
        #     L(X) direction in Lb rest frame
        #     ALSO: get proton momentum in Lb frame
        #     -----------------------------------------------

        if args.generated:
            boost_rp_TRUE = -v_resonance_TRUE.BoostVector() # boost vector to res RF in Lb RF
            v_c1_TRUE.Boost(boost_p_TRUE) # boost child to parent RF
            v_c1_TRUE.Boost(boost_rp_TRUE) # boost to resonance RF in Lb RF
            u_c1_TRUE = v_c1_TRUE.Vect().Unit() # get unit vector
            br_var[len(variables)+5][0] = u_c1_TRUE.Dot(u_r_TRUE) # get cos angle between the two
            br_var[len(variables)+6][0] = v_c1_TRUE.P() # get proton momentum in resonance frame

        u_c1 = v_c1.Vect().Unit() # get unit vector
        br_var[5][0] = u_c1.Dot(u_r) # get cos angle between the two
        br_var[6][0] = v_c1.P() # get proton momentum in resonance frame

         #-----------------------------------------------
         #   compute cosTheta_lepton:
         #    angle between lepton momentum in Dilepton rest frame and
         #    Dilepton direction in Lb rest frame
         #    ALSO: get lepton momentum in Lb frame
         #-----------------------------------------------

        if args.generated:
            boost_dlp_TRUE = -v_dilepton_TRUE.BoostVector() # boost vector to dilepton RF in Lb RF
            v_c3_TRUE.Boost(boost_p_TRUE) # boost child to parent RF
            v_c3_TRUE.Boost(boost_dlp_TRUE) # boost to dilepton RF in Lb RF
            u_c3_TRUE = v_c3_TRUE.Vect().Unit() # get unit vector

            br_var[len(variables)+8][0] = u_c3_TRUE.Dot(u_dl_TRUE) # get cos angle between the two CHECK HERE !!!!
            br_var[len(variables)+10][0] = v_c3_TRUE.P() # get proton momentum in resonance frame
        boost_dlp = -v_dilepton.BoostVector() # boost vector to dilepton RF in Lb RF
        u_c3 = v_c3.Vect().Unit() # get unit vector
        br_var[8][0] = u_c3.Dot(u_dl) # get cos angle between the two
        br_var[10][0] = v_c3.P() # get lepton momentum in dilepton frame

        #---------------------
        # compute phi
        #---------------------
        if (i==0):
            print(" # Resonance momentum in parent RF:")
            v_resonance.Print()
            print(" # Dilepton momentum in parent RF:")
            v_dilepton.Print()
            print(" # Child1 momentum in resonance RF:")
            v_c1.Print()
            print(" # Child2 momentum in resonance RF:")
            v_c2.Print()
            print(" # Child3 momentum in dilepton RF:")
            v_c3.Print()
            print(" # Child4 momentum in dilepton RF:")
            v_c4.Print()

        decPlane_12 = (u_r.Cross(u_c1)).Unit()
        decPlane_34 = (u_dl.Cross(u_c3)).Unit()

        cosphi = (decPlane_12.Dot(decPlane_34))
        sinphi = ( decPlane_34.Cross( decPlane_12 ) ).Dot( u_dl )
        phi    = acos( cosphi )
        # Resolve ambiguity
        if sinphi > 0.0:
        
            br_var[2][0] = phi
        else:
            br_var[2][0] = -phi

        if (i==0):
            print("Normal to decay plane 12:")
            decPlane_12.Print()
            print("Normal to decay plane 34:")
            decPlane_34.Print()
            print("Phi:", br_var[2][0])

        # ------------------
        #    fill tree branches
        # ------------------
        for br in br_bra:
            br.Fill()

    # Close the file
    f.Write('',ROOT.TFile.kOverwrite)
    f.Close()
# EOF
