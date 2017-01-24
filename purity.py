#!/usr/bin/env python3.4

import math
from ROOT import TFile, TChain, TH1F, TH3F, gDirectory, TCanvas, kRed
from geomUtil import isect_line_plane_v3
from glob import glob

bin_ang = 12
fidvol = 20
bin_len = 8

f = glob("../root_files/mucs_cry_overlay/*/ana_hist.root")

h_res_phi = TH1F("h_res_phi","",100,-10,10)
h_res_theta = TH1F("h_res_theta","",100,-10,10)

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
decay = 0
ok = 0
entries = entry.GetEntries()
print(entries)
for event in range(entries):

    if event % 100 == 0: print(event)
    ientry = entry.LoadTree(event)
    nb = entry.GetEntry(event)
    geant_trackIds = [entry.TrackId[i] for i in range(entry.geant_list_size)]

    for i in range(entry.geant_list_size):
        mucs_event = entry.StartPointy[i] > 391 and entry.StartPointy[i] < 399
        if mucs_event and abs(entry.pdg[i]) == 13 and entry.process_primary[i] == 1:
            if entry.EndPointy[i] > y_end:
                print(entry.EndPointy[i],entry.Eng[i],"doesn't reach TPC")
                decay+=1
            else:
                ok+=1
            break
        if entry.inTPCActive[i] and abs(entry.pdg[i]) == 13:
            p0 = [entry.StartPointx[i],entry.StartPointy[i],entry.StartPointz[i]]
            p1 = [entry.StartPointx[i]+entry.Px[i],entry.StartPointy[i]+entry.Py[i],entry.StartPointz[i]+entry.Pz[i]]

            p_tpc = isect_line_plane_v3(p0, p1, p_co_tpc, p_no)
            mucs_event = entry.StartPointy[i] > 391 and entry.StartPointy[i] < 399

            x = entry.Px[i]
            y = entry.Py[i]
            z = entry.Pz[i]
            r = entry.P[i]
            theta = math.degrees(math.acos(z/r))
            phi = math.degrees(math.atan2(y,x))
            l = (entry.StartPointx_tpcAV[i]-entry.EndPointx_tpcAV[i])**2+(entry.StartPointy_tpcAV[i]-entry.EndPointy_tpcAV[i])**2+(entry.StartPointz_tpcAV[i]-entry.EndPointz_tpcAV[i])**2
            l = math.sqrt(l)
            if mucs_event and p_tpc and theta < 75:
                min_dist = 99999
                min_id = 0
                found = 0
                in_tpc += 1
                i_min = 0
                for n in range(entry.ntracks_pandoraCosmic):

                    best = entry.trkpidbestplane_pandoraCosmic[n]
                    trkid = entry.trkidtruth_pandoraCosmic[3*n+best]
                    p_reco = [entry.trkstartx_pandoraCosmic[n],entry.trkstarty_pandoraCosmic[n],entry.trkstartz_pandoraCosmic[n]]
                    dist = math.sqrt(sum([(x1-x2)**2 for x1,x2 in zip(p_tpc,p_reco)]))

                    if trkid == entry.TrackId[i]:
                        found = 1

                    if dist < min_dist:
                        i_min = n
                        min_dist = dist
                        min_id = trkid


                if not found:
                    print(f[int(event/25)], event%25, event, entry.event)
                reco_theta = math.degrees(entry.trktheta_pandoraCosmic[i_min])
                reco_phi = math.degrees(entry.trkphi_pandoraCosmic[i_min])

                if min_dist < 32 and abs(theta-reco_theta) < 2 and abs(phi-reco_phi) < 2:
                    reco += 1

                    h_res_theta.Fill(theta-reco_theta)
                    h_res_phi.Fill(phi-reco_phi)
                    if min_id == entry.TrackId[i]:
                        mucs_reco += 1

                    break
print(in_tpc)
print(decay,ok)
eff = reco/in_tpc
eff_err = math.sqrt((eff*(1-eff))/in_tpc)
print("Efficiency: %.2f +- %.2f"%(round(eff*100,2), round(eff_err*100,2)))

purity = mucs_reco/reco
purity_err = math.sqrt((purity*(1-purity))/reco)
print("Purity: %.2f +- %.2f"%(round(purity*100,2), round(purity_err*100,2)))
acceptance = purity*eff/0.9865
print(acceptance)
print(0.9865/eff)

c_theta = TCanvas("c_theta","")
h_res_theta.Draw()
c_theta.Update()
c_phi = TCanvas("c_phi","")
h_res_phi.Draw()
c_phi.Update()
input()
