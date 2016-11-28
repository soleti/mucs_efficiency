#!/usr/bin/env python3.4

import os, sys
import math
from ROOT import TPad, gDirectory, TChain, TFile, TH1F, TH2F, TH3F, gStyle, TCanvas, kRed, TProfile, gPad, TPaveText
from glob import glob
from geomUtil import *
import random

gStyle.SetOptStat(0)
gStyle.SetPalette(87)
gStyle.SetPaintTextFormat(".2f")
gStyle.SetNumberContours(999)

algo = "pandoraCosmic"

fidvol = 20
x_start = 0
x_end = 256.35
y_start = -116.5
y_end = 116.5
z_start = 0
z_end = 1036.8
top = [-40,40,399.45,398.45,560,660]
bottom = [-80,0,117.5,116.5,560,660]

bin_ang = 12
bin_len = 8

h_eff = TH1F("h_eff",";Efficiency;N. Entries / 0.03",35,0,1.5)      


files = [glob("../root_files/energy_sampling/MuCSCryMCGen/*/ana_hist.root")]
        
triggered = 0
reco = 0
p_bottom = [0,y_start,0]
p_no_bottom = [0,1,0]
p_anode = [x_end,y_start,0]
p_no_anode = [1,0,0]
p_no = [0,1,0]
p_co_top = [top[0], top[2], top[4]]
p_co_bottom = [bottom[0], bottom[2], bottom[4]]
tot = files
geant_trackIds = []

for n,f in enumerate(tot):
    print(n)
    chain = TChain("analysistree/anatree")
    for filename in f:
        chain.Add(filename)
    entries = chain.GetEntries()
    trk=0

    for entry in range(entries):

        ientry = chain.LoadTree(entry)
        nb = chain.GetEntry(entry)

        for i in range(chain.geant_list_size):
            p0 = [chain.StartPointx_drifted[i], chain.StartPointy_drifted[i], chain.StartPointz_drifted[i]]
            p2 = [chain.StartPointx_drifted[i]+chain.Px[i], chain.StartPointy_drifted[i]+chain.Py[i], chain.StartPointz_drifted[i]+chain.Pz[i]]

            theta_data, phi_data = math.degrees(chain.theta[i]), math.degrees(chain.phi[i])
            theta_ok = theta_data > 60 and theta_data < 75
            phi_ok = phi_data > -90 and phi_data < -75
            
            if chain.inTPCActive[i] and (abs(chain.pdg[i]) == 13) and chain.pathlen[i] > 140 and chain.pathlen[i] < 200 and theta_ok and phi_ok:
                triggered+=1
                geant_trackIds.append(chain.TrackId[i])
                                                    

        for j in range(chain.ntracks_pandoraCosmic):
            
            x1 = chain.trkstartx_pandoraCosmic[j]
            y1 = chain.trkstarty_pandoraCosmic[j]
            z1 = chain.trkstartz_pandoraCosmic[j]

            p0 = [x1,y1,z1]
            
            best = chain.trkpidbestplane_pandoraCosmic[j]
            trkid = chain.trkidtruth_pandoraCosmic[3*j+best]
            if trkid in geant_trackIds:
                index = geant_trackIds.index(trkid)
                p0_geant = [chain.StartPointx_tpcAV[index], chain.StartPointy_tpcAV[index], chain.StartPointz_tpcAV[index]]
                dist = math.sqrt(sum([(a-b)**2 for a,b in zip(p0,p0_geant)]))                
                theta_data, phi_data = math.degrees(chain.theta[index]),math.degrees(chain.phi[index])
                
                if dist < 35 and theta_data > 60 and theta_data < 75 and phi_data > -90 and phi_data < -75 and chain.pathlen[index] > 140 and chain.pathlen[index] < 200:
                    reco+=1

elements = (triggered-reco+56)*[0] + reco*[1]

print(elements)
for i in range(1000):
    print(sum(random.sample(elements,25))/25)
    h_eff.Fill(sum(random.sample(elements,25))/25)
    

print(reco,triggered)
eff = reco/(triggered+56)
error = math.sqrt((eff*(1-eff))/triggered) 
print("Integrated efficiency",eff,error)


c_eff = TCanvas("c_eff")
h_eff.Draw("ep")
h_eff.SetMarkerStyle(20)
h_eff.Fit("gaus")
c_eff.Update()   
input()