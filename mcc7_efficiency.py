#!/usr/bin/env python3.4

from ROOT import TPad, gPad, gDirectory, TFile, gStyle, TChain, TH1F, TH2F, TH3F, TLine, TCanvas, TPad, kRed, kGray, kBlue, TBox
import math

bin_ang = 12
fidvol = 20
bin_len = 8

def getEff(h_num, h_den):
    h_num.Divide(h_den)
    
    if type(h_num) == TH2F:
        for i in range(1,h_num.GetNbinsX()+2):
            for j in range(1,h_num.GetNbinsY()+2):
                den = h_den.GetBinContent(i, j)         
                eff = h_num.GetBinContent(i, j)
                    
                if den and eff <= 1:
                    error = math.sqrt((eff*(1-eff))/den)
                    h_num.SetBinError(i, j, error)


    if type(h_num) == TH1F:
        for i in range(1,h_num.GetNbinsX()+2):
            den = h_den.GetBinContent(i)  
            eff = h_num.GetBinContent(i)

            if den and eff <=1 :
                error = math.sqrt((eff*(1-eff))/den)
                h_num.SetBinError(i, error)

            
    return h_num
            

gStyle.SetOptStat(0)
gStyle.SetPalette(87)
gStyle.SetPaintTextFormat(".2f")
gStyle.SetNumberContours(999)

f = TFile("../root_files/mcc7_ana.root")
chain = gDirectory.Get("analysistree/anatree")
entries = chain.GetEntries()
tot_prim = 0
tot_reco = 0

h_theta_phi_geant = TH2F("h_theta_phi_geant",";#theta [#circ];#phi [#circ]",bin_ang,0,180,bin_ang,-180,0)
h_theta_phi_reco = TH2F("h_theta_phi_reco",";#theta [#circ];#phi [#circ]",bin_ang,0,180,bin_ang,-180,0)

h_theta_l_geant = TH2F("h_theta_l_geant",";#theta [#circ];L [cm]",bin_ang,0,180,bin_len,fidvol,500)
h_theta_l_reco = TH2F("h_theta_l_reco",";#theta [#circ];L [cm]",bin_ang,0,180,bin_len,fidvol,500)

h_phi_l_geant = TH2F("h_phi_l_geant",";#phi [#circ];L [cm]",bin_ang,-180,0,bin_len,fidvol,500)
h_phi_l_reco = TH2F("h_phi_l_reco",";#phi [#circ];L [cm]",bin_ang,-180,0,bin_len,fidvol,500)

h_theta_reco = TH1F("h_theta_reco",";#theta [#circ]; N. Entries / 10#circ",bin_ang,0,180)
h_theta_geant = TH1F("h_theta_geant",";#theta [#circ]; N. Entries / 10#circ",bin_ang,0,180)

h_phi_reco = TH1F("h_phi_reco",";#phi [#circ]; N. Entries / 20#circ",bin_ang,-180,0)
h_phi_geant = TH1F("h_phi_geant",";#phi [#circ]; N. Entries / 20#circ",bin_ang,-180,0)

h_l_geant = TH1F("h_l_geant",";L [cm]; N. Entries / 80 cm",bin_len,fidvol,500)
h_l_reco = TH1F("h_l_reco",";L [cm]; N. Entries / 80 cm",bin_len,fidvol,500)

h_theta_phi_l_reco = TH3F("h_theta_phi_l_reco",";#theta [#circ]; #phi [#circ]; L [cm]",bin_ang,0,180,bin_ang,-180,0,bin_len,fidvol,500)
h_theta_phi_l_geant = TH3F("h_theta_phi_l_geant",";#theta [#circ]; #phi [#circ]; L [cm]",bin_ang,0,180,bin_ang,-180,0,bin_len,fidvol,500)


print(entries)
for entry in range(200):
    
    if entry % 10 == 0: print(entry)
    
    ientry = chain.LoadTree(entry)
    nb = chain.GetEntry(entry)
    geant_trackIds = [chain.TrackId[i] for i in range(chain.geant_list_size)]
    primaries = 0
    reco = 0
    
    for i in range(chain.geant_list_size):
        if chain.process_primary[i] == 1 and chain.inTPCDrifted[i] == 1:
            x = chain.Px[i]
            y = chain.Py[i]
            z = chain.Pz[i]
            r = chain.P[i]
            theta = math.degrees(math.acos(z/r))
            phi = math.degrees(math.atan2(y,x))
            l = (chain.StartPointx_tpcAV[i]-chain.EndPointx_tpcAV[i])**2+(chain.StartPointy_tpcAV[i]-chain.EndPointy_tpcAV[i])**2+(chain.StartPointz_tpcAV[i]-chain.EndPointz_tpcAV[i])**2
            l = math.sqrt(l)
            h_theta_phi_l_geant.Fill(theta, phi, l)
            primaries += 1
            
    indeces = []

    for i in range(chain.ntracks_pandoraCosmic):
        best = chain.trkpidbestplane_pandoraCosmic[i]
        trkid = chain.trkidtruth_pandoraCosmic[3*i+best]
        
        #if not (chain.trkidtruth_pandoraCosmic[3*i] == chain.trkidtruth_pandoraCosmic[3*i+1] and chain.trkidtruth_pandoraCosmic[3*i] == chain.trkidtruth_pandoraCosmic[3*i+2]):
        #    if chain.trkidtruth_pandoraCosmic[3*i] != -1 and chain.trkidtruth_pandoraCosmic[3*i+1] != -1 and chain.trkidtruth_pandoraCosmic[3*i+2] != -1:
        #        print(chain.trkidtruth_pandoraCosmic[3*i],chain.trkidtruth_pandoraCosmic[3*i+1],chain.trkidtruth_pandoraCosmic[3*i+2])
        
        if trkid in geant_trackIds:
            index = geant_trackIds.index(trkid)
            
            x1 = chain.trkstartx_pandoraCosmic[i]
            y1 = chain.trkstarty_pandoraCosmic[i]
            z1 = chain.trkstartz_pandoraCosmic[i]
            p0 = [x1,y1,z1]
            
            p0_geant = [chain.StartPointx_drifted[index], chain.StartPointy_drifted[index], chain.StartPointz_drifted[index]]
            dist = math.sqrt(sum([(a-b)**2 for a,b in zip(p0,p0_geant)]))
            
            if index not in indeces and chain.inTPCDrifted[index] == 1 and chain.process_primary[index] == 1 and dist < 35:
                reco += 1
                x = chain.Px[index]
                y = chain.Py[index]
                z = chain.Pz[index]
                r = chain.P[index]
                theta = math.degrees(math.acos(z/r))
                phi = math.degrees(math.atan2(y,x))
                
                l = (chain.StartPointx_tpcAV[index]-chain.EndPointx_tpcAV[index])**2+(chain.StartPointy_tpcAV[index]-chain.EndPointy_tpcAV[index])**2+(chain.StartPointz_tpcAV[index]-chain.EndPointz_tpcAV[index])**2
                l = math.sqrt(l)
                
                h_theta_phi_l_reco.Fill(theta, phi, l)
                indeces.append(index)
                       
        
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
            
print("Overall efficiency: %.1f%%" % round(reco/primaries*100,1))

c_theta_phi = TCanvas("c_theta_phi")
for i in range(1,h_theta_phi_l_reco.GetNbinsX()+2):
    for j in range(1,h_theta_phi_l_reco.GetNbinsY()+2):
        reco = sum([h_reco.GetBinContent(i,j,k) for k in range(1,h_reco.GetNbinsZ()+2)])
        geant = sum([h_geant.GetBinContent(i,j,k) for k in range(1,h_geant.GetNbinsZ()+2)])
        if reco and geant:
            eff = reco/geant
            error = math.sqrt((eff*(1-eff))/geant)
            h_theta_phi_reco.SetBinContent(i,j,eff)
            h_theta_phi_reco.SetBinError(i,j,error)
            
h_theta_phi_reco.Draw("colz texte")
h_theta_phi_reco.SaveAs("plots/mc/theta_phi_mcc7.root")
c_theta_phi.Update()

c_theta_l = TCanvas("c_theta_l")
for i in range(1,h_theta_phi_l_reco.GetNbinsX()+2):
    for j in range(1,h_theta_phi_l_reco.GetNbinsZ()+2):
        reco = sum([h_reco.GetBinContent(i,k,j) for k in range(1,h_reco.GetNbinsY()+2)])
        geant = sum([h_geant.GetBinContent(i,k,j) for k in range(1,h_geant.GetNbinsY()+2)])
        if reco and geant:
            eff = reco/geant
            error = math.sqrt((eff*(1-eff))/geant)
            h_theta_l_reco.SetBinContent(i,j,eff)
            h_theta_l_reco.SetBinError(i,j,error)
h_theta_l_reco.Draw("colz texte")
h_theta_l_reco.SaveAs("plots/mc/theta_l_mcc7.root")
c_theta_l.Update()

c_phi_l = TCanvas("c_phi_l")
for i in range(1,h_theta_phi_l_reco.GetNbinsY()+2):
    for j in range(1,h_theta_phi_l_reco.GetNbinsZ()+2):
        reco = sum([h_reco.GetBinContent(k,i,j) for k in range(1,h_reco.GetNbinsX()+2)])
        geant = sum([h_geant.GetBinContent(k,i,j) for k in range(1,h_geant.GetNbinsX()+2)])
        if reco and geant:
            eff = reco/geant
            error = math.sqrt((eff*(1-eff))/geant)
            h_phi_l_reco.SetBinContent(i,j,eff)
            h_phi_l_reco.SetBinError(i,j,error)
h_phi_l_reco.Draw("colz texte")
h_phi_l_reco.SaveAs("plots/mc/phi_l_mcc7.root")
c_phi_l.Update()

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

c_theta = TCanvas("c_theta")
for i in range(1,h_theta_phi_l_reco.GetNbinsX()+2): 
    reco = sum([h_reco.GetBinContent(i,j,k) for j in range(1,h_reco.GetNbinsY()+2) for k in range(1,h_reco.GetNbinsZ()+2)])
    geant = sum([h_geant.GetBinContent(i,j,k) for j in range(1,h_geant.GetNbinsY()+2) for k in range(1,h_geant.GetNbinsZ()+2)])
    
    if reco and geant:
        eff = reco/geant
        error = math.sqrt((eff*(1-eff))/geant)
        h_theta_reco.SetBinContent(i,eff)
        h_theta_reco.SetBinError(i,error)
        h_theta_reco.SetLineColor(kRed+1)

h_theta_reco.SetLineColor(kRed+1)
h_theta_reco.Draw("hist")
h_theta_reco.GetYaxis().SetRangeUser(0.001,1.05)
h_theta_reco_clone = h_theta_reco.Clone()
h_theta_reco_clone.SetFillStyle(3002)
h_theta_reco_clone.SetFillColor(kRed+1)
h_theta_reco_clone.Draw("e2same")
h_theta_reco_clone.SaveAs("plots/mc/theta_mcc7.root")

c_theta.Update()

c_phi = TCanvas("c_phi")
for i in range(1,h_theta_phi_l_reco.GetNbinsY()+2): 
    reco = sum([h_reco.GetBinContent(j,i,k) for j in range(1,h_reco.GetNbinsX()+2) for k in range(1,h_reco.GetNbinsZ()+2)])
    geant = sum([h_geant.GetBinContent(j,i,k) for j in range(1,h_geant.GetNbinsX()+2) for k in range(1,h_geant.GetNbinsZ()+2)])
    
    if reco and geant:
        eff = reco/geant
        error = math.sqrt((eff*(1-eff))/geant)
        h_phi_reco.SetBinContent(i,eff)
        h_phi_reco.SetBinError(i,error)
        h_phi_reco.SetLineColor(kRed+1)

h_phi_reco.SetLineColor(kRed+1)
h_phi_reco.Draw("hist")
h_phi_reco.GetYaxis().SetRangeUser(0.001,1.05)
h_phi_reco_clone = h_phi_reco.Clone()
h_phi_reco_clone.SetFillStyle(3002)
h_phi_reco_clone.SetFillColor(kRed+1)
h_phi_reco_clone.Draw("e2same")
h_phi_reco_clone.SaveAs("plots/mc/phi_mcc7.root")

c_phi.Update()
input()