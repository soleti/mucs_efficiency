#!/usr/bin/env python3.4

import math
from ROOT import TFile, TChain, TH1F, TH3F, gDirectory, TCanvas, kRed
from geomUtil import isect_line_plane_v3
from glob import glob

bin_ang = 12
fidvol = 20
bin_len = 8
        
f = glob("../root_files/mucs_cry_overlay_upstream/*/ana_hist.root")

entry = TChain("analysistree/anatree")
for filename in f:
    entry.Add(filename)
    
y_end = 116.5
        
p_no = [0,1,0]
p_co_tpc = [0,y_end,0]
    
in_tpc = 0
mucs_reco = 0
reco = 0
total_mucs_reco = 0

entries = entry.GetEntries()
print(entries)
for event in range(entries):

    if event % 100 == 0: print(event)
    ientry = entry.LoadTree(event)
    nb = entry.GetEntry(event) 
    geant_trackIds = [entry.TrackId[i] for i in range(entry.geant_list_size)]
        
    for i in range(entry.geant_list_size):
        if entry.inTPCActive[i] and abs(entry.pdg[i]) == 13:
            p0 = [entry.StartPointx[i],entry.StartPointy[i],entry.StartPointz[i]]
            p1 = [entry.StartPointx[i]+entry.Px[i],entry.StartPointy[i]+entry.Py[i],entry.StartPointz[i]+entry.Pz[i]]

            p_tpc = isect_line_plane_v3(p0, p1, p_co_tpc, p_no)
            mucs_event = entry.StartPointy[i] > 391 and entry.StartPointy[i] < 392

            if mucs_event and p_tpc:
                min_dist = 99999
                min_id = 0
                in_tpc += 1
                found = 0
                x = entry.Px[i]
                y = entry.Py[i]
                z = entry.Pz[i]
                r = entry.P[i]
                theta = math.degrees(math.acos(z/r))
                phi = math.degrees(math.atan2(y,x))
                l = (entry.StartPointx_tpcAV[i]-entry.EndPointx_tpcAV[i])**2+(entry.StartPointy_tpcAV[i]-entry.EndPointy_tpcAV[i])**2+(entry.StartPointz_tpcAV[i]-entry.EndPointz_tpcAV[i])**2
                l = math.sqrt(l)
                h_theta_phi_l_geant.Fill(theta, phi, l)


                for n in range(entry.ntracks_pandoraCosmic):
                    
                    best = entry.trkpidbestplane_pandoraCosmic[n]
                    trkid = entry.trkidtruth_pandoraCosmic[3*n+best]
                    p_reco = [entry.trkstartx_pandoraCosmic[n],entry.trkstarty_pandoraCosmic[n],entry.trkstartz_pandoraCosmic[n]]
                    dist = math.sqrt(sum([(x1-x2)**2 for x1,x2 in zip(p_tpc,p_reco)]))
                    if trkid == entry.TrackId[i]:
                        found = 1
                    if dist < min_dist:
                        min_dist = dist
                        min_id = trkid
                
                
                if not found:
                    print(f[int(event/25)], event%25, event, entry.event)

                if min_dist < 32:
                    reco += 1                    
                    if min_id == entry.TrackId[i]:
                        mucs_reco += 1 

                    break

                
eff = reco/in_tpc
eff_err = math.sqrt((eff*(1-eff))/in_tpc)
print("Efficiency: %.2f +- %.2f"%(round(eff*100,2), round(eff_err*100,2)))

purity = mucs_reco/reco
purity_err = math.sqrt((purity*(1-purity))/reco)
print("Purity: %.2f +- %.2f"%(round(purity*100,2), round(purity_err*100,2)))


input()
