#!/usr/bin/env python3.4

import math
from ROOT import TFile, TChain
from geomUtil import isect_line_plane_v3
from glob import glob
        
f = glob("../mucs_cry_overlay/*/ana_hist.root")    
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
for event in range(1000):
    
    if event % 100 == 0: print(event)
    ientry = entry.LoadTree(event)
    nb = entry.GetEntry(event) 
    geant_trackIds = [entry.TrackId[i] for i in range(entry.geant_list_size)]
        
    for i in range(entry.geant_list_size):
        if entry.inTPCActive[i] and abs(entry.pdg[i]) == 13:
            p0 = [entry.StartPointx[i],entry.StartPointy[i],entry.StartPointz[i]]
            p1 = [entry.StartPointx[i]+entry.Px[i],entry.StartPointy[i]+entry.Py[i],entry.StartPointz[i]+entry.Pz[i]]

            p_tpc = isect_line_plane_v3(p0, p1, p_co_tpc, p_no)
            mucs_event = entry.StartPointy[i] > 397 and entry.StartPointy[i] < 398

            if mucs_event and p_tpc:
                min_dist = 99999
                min_id = 0
                in_tpc += 1
                
                for n in range(entry.ntracks_pandoraCosmic):
                    
                    best = entry.trkpidbestplane_pandoraCosmic[n]
                    trkid = entry.trkidtruth_pandoraCosmic[3*n+best]
                    p_reco = [entry.trkstartx_pandoraCosmic[n],entry.trkstarty_pandoraCosmic[n],entry.trkstartz_pandoraCosmic[n]]
                    dist = math.sqrt(sum([(x1-x2)**2 for x1,x2 in zip(p_tpc,p_reco)]))
                    if dist < min_dist:
                        min_dist = dist
                        min_id = trkid
                        trk_n = n
                             
                if min_dist < 35:
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

#reco_eff = total_mucs_reco/in_tpc
#reco_eff_err = math.sqrt((reco_eff*(1-reco_eff))/in_tpc)
#print("Total MuCS reco. efficiency: %.2f +- %.2f"%(round(reco_eff*100,2), round(reco_eff_err*100,2)))

input()
