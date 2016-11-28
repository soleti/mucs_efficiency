#!/usr/bin/env python3.4

import math
from ROOT import TFile, TChain, TH1F, TH3F, gDirectory, TCanvas, kRed
from geomUtil import isect_line_plane_v3
from glob import glob

bin_ang = 12
fidvol = 20
bin_len = 8

h_l_geant = TH1F("h_l_geant",";L [cm]; N. Entries / 80 cm",bin_len,fidvol,500)
h_l_reco = TH1F("h_l_reco",";L [cm]; N. Entries / 80 cm",bin_len,fidvol,500)

h_theta_phi_l_reco = TH3F("h_theta_phi_l_reco",";#theta [#circ]; #phi [#circ]; L [cm]",bin_ang,0,180,bin_ang,-180,0,bin_len,fidvol,500)
h_theta_phi_l_geant = TH3F("h_theta_phi_l_geant",";#theta [#circ]; #phi [#circ]; L [cm]",bin_ang,0,180,bin_ang,-180,0,bin_len,fidvol,500)

        
f = glob("../mucs_cry_overlay_upstream/*/ana_hist.root")

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
bad_angle = 0
bad_angle_min = 0
bad_angle_geant = 0
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

                if theta > 60 and theta < 75 and phi > -60 and phi < -45:
                    bad_angle_geant+=1

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
                
                if entry.TrackId[i] in entry.trkidtruth_pandoraCosmic and theta > 60 and theta < 75 and phi > -60 and phi < -45:
                    bad_angle += 1
                
                if not found:
                    #print(entry.StartPointx[i], entry.StartPointz[i])
                    #print(entry.TrackId[i], [entry.trkidtruth_pandoraCosmic[3*entry.trkpidbestplane_pandoraCosmic[n]+best] for n in range(entry.ntracks_pandoraCosmic)])
                    print(f[int(event/25)], event%25, event, entry.event)

                if min_dist < 35:
                    reco += 1                    
                    if min_id == entry.TrackId[i]:
                        mucs_reco += 1 
                        h_theta_phi_l_reco.Fill(theta, phi, l)

                    if theta > 60 and theta < 75 and phi > -60 and phi < -45:
                        bad_angle_min += 1

                    break

#print(bad_angle,bad_angle_geant,bad_angle_min)               
#print(bad_angle/bad_angle_min)                  
eff = reco/in_tpc
eff_err = math.sqrt((eff*(1-eff))/in_tpc)
print("Efficiency: %.2f +- %.2f"%(round(eff*100,2), round(eff_err*100,2)))

purity = mucs_reco/reco
purity_err = math.sqrt((purity*(1-purity))/reco)
print("Purity: %.2f +- %.2f"%(round(purity*100,2), round(purity_err*100,2)))

#reco_eff = total_mucs_reco/in_tpc
#reco_eff_err = math.sqrt((reco_eff*(1-reco_eff))/in_tpc)
#print("Total MuCS reco. efficiency: %.2f +- %.2f"%(round(reco_eff*100,2), round(reco_eff_err*100,2)))

f_theta_phi_l = TFile("plots/data/e_theta_phi_l_pandoraCosmic.root")
e_theta_phi_l = gDirectory.Get("h_theta_phi_l_tpc")
h_reco = h_theta_phi_l_reco.Clone()
h_geant = h_theta_phi_l_geant.Clone()
for i in range(1, h_theta_phi_l_reco.GetNbinsX()+2):
    for j in range(1, h_theta_phi_l_reco.GetNbinsY()+2):
        for k in range(1, h_theta_phi_l_reco.GetNbinsZ()+2):
            if h_theta_phi_l_reco.GetBinContent(i,j,k) and h_theta_phi_l_geant.GetBinContent(i,j,k):
                if e_theta_phi_l.GetBinContent(i,j,k) > 0:
                    eff = h_theta_phi_l_reco.GetBinContent(i,j,k)/h_theta_phi_l_geant.GetBinContent(i,j,k)
                    error = math.sqrt((eff*(1-eff))/h_theta_phi_l_geant.GetBinContent(i,j,k))
                else:
                    eff = 0
                    error = 0
                    h_reco.SetBinContent(i,j,k,0)
                    h_reco.SetBinError(i,j,k,0)
                    h_geant.SetBinContent(i,j,k,0)
                    h_geant.SetBinError(i,j,k,0)
                h_theta_phi_l_reco.SetBinContent(i,j,k,eff)
                h_theta_phi_l_reco.SetBinError(i,j,k,error)

c_l = TCanvas("c_l")
for i in range(1,h_theta_phi_l_reco.GetNbinsZ()+2): 
    reco = sum([h_reco.GetBinContent(j,k,i) for j in range(1,h_reco.GetNbinsX()+2) for k in range(1,h_reco.GetNbinsY()+2)])
    geant = sum([h_geant.GetBinContent(j,k,i) for j in range(1,h_geant.GetNbinsX()+2) for k in range(1,h_geant.GetNbinsY()+2)])
    
    if reco and geant:
        eff = reco/geant
        error = math.sqrt((eff*(1-eff))/geant)
        h_l_reco.SetBinContent(i,eff)
        h_l_reco.SetBinError(i,error)
        h_l_reco.SetLineColor(kRed+1)
        
h_l_reco.Draw("hist")
h_l_reco.GetYaxis().SetRangeUser(0.001,1.05)
h_l_reco_clone = h_l_reco.Clone()
h_l_reco_clone.SetFillStyle(3002)
h_l_reco_clone.SetFillColor(kRed+1)
h_l_reco_clone.Draw("e2same")
h_l_reco_clone.SaveAs("plots/mc/l_mcc7.root")     
c_l.Update()
input()
