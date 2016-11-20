#!/usr/bin/env python3.4

import os, sys
import math
from ROOT import TPad, gDirectory, TChain, TFile, TH1F, TH2F, TH3F, gStyle, TCanvas, gPad, TPaveText, kGray

if len(sys.argv) > 1:
    data = 1 if sys.argv[1] == "data" else 0
else:
    data = 1

if len(sys.argv) > 2:
    position = sys.argv[2]
else:
    position = "central_shift"

if len(sys.argv) > 3:
    algo = sys.argv[3]
else:
    algo = "pandoraCosmic"

if data: print("Data efficiency")
else: print("Monte Carlo efficiency")
print(position, "dataset")
print(algo, "algorithm")

def getThetaPhiFromThetas(qxy,qyz):
   tqxy=math.tan(qxy)
   tqyz=math.tan(qyz)
   theta=-math.pi/2+math.acos(((1+tqxy**-2+tqyz**-2)**-0.5))
   phi= math.atan2(math.sin(theta)/math.tan(qyz),math.sin(theta)/math.tan(qxy))
   return theta,phi

def weighted_mean(values, weights):
    return sum(v*w for v,w in zip(values,weights))/sum(weights)

def extrapolate3D(hist):
    hist_ex = hist.Clone()
    hist_ex.SetName(hist.GetName()+"_ex")
    
    for i in range(hist.GetNbinsX()+2):
        for j in range(hist.GetNbinsY()+2):
            for k in range(hist.GetNbinsZ()+2):
                #New map with distance:
                if not hist.GetBinContent(i,j,k):
                    distances = []
                    eff = []
                    eff_err = []
                    for l in range(hist.GetNbinsX()+2):
                        for m in range(hist.GetNbinsY()+2):
                            for n in range(hist.GetNbinsZ()+2):
                                distances += [(l-i)**2+(j-m)**2+(k-n)**2]
                                eff += [hist.GetBinContent(l,m,n)]
                                eff_err += [hist.GetBinError(l,m,n)]
                            
                
                    minimum = min([d for i,d in enumerate(distances) if eff[i] > 0])
                    closest_eff = [eff[a] for a,b in enumerate(distances) if b == minimum and eff[a] > 0]
                    closest_eff_err = [eff_err[a] for a,b in enumerate(distances) if b == minimum and eff[a] > 0]
                    weights = [1/(err**2) for err in closest_eff_err]

                    eff_ave = weighted_mean(closest_eff, weights)
                    hist_ex.SetBinContent(i,j,k,eff_ave)
                    hist_ex.SetBinError(i,j,k,1/sum(weights))
                    
    return hist_ex
    

gStyle.SetOptStat(0)
gStyle.SetPalette(87)
gStyle.SetPaintTextFormat(".2f")
gStyle.SetNumberContours(999)


chain = TChain("events_"+algo)

#chain.Add("MuCSRun3702_Group158_pdfrecoBigExtension_MergedTree.Root")
#chain.Add("MuCSRun3702_Group158500MeV_MergedTree.Root")

#central_centre == 33
#central_shift == 30
#upstream_shift == 35
#downstream_shift == 34

if position == "central_centre":
    chain.Add("../mucs_merged/MuCSRun7263_Group180_MergedTree.Root")
    chain.Add("../mucs_merged/MuCSRun7264_Group180_MergedTree.Root")
    chain.Add("../mucs_merged/MuCSRun7265_Group180_MergedTree.Root")
if position == "downstream_shift":
    chain.Add("../mucs_merged/MuCSRun7347_Group181_MergedTree.Root")
    chain.Add("../mucs_merged/MuCSRun7348_Group181_MergedTree.Root")
if position == "central_shift":
    chain.Add("../mucs_merged/CORSIKA.root")
if position == "upstream_shift":
    chain.Add("../mucs_merged/MuCSRun7702_Group182_MergedTree.Root")
    chain.Add("../mucs_merged/MuCSRun7703_Group183_MergedTree.Root")
    
fidvol = 20
x_start = 0
x_end = 256.35
y_start = -116.5
y_end = 116.5
z_start = 0
z_end = 1036.8

bin_ang = 12
bin_len = 8
h_theta_phi_l_tpc = TH3F("h_theta_phi_l_tpc",";#theta [#circ]; #phi [#circ]; L [cm]",bin_ang,0,180,bin_ang,-180,0,bin_len,fidvol,500)
h_theta_phi_l_mucs = TH3F("h_theta_phi_l_mucs",";#theta [#circ]; #phi [#circ]; L [cm]",bin_ang,0,180,bin_ang,-180,0,bin_len,fidvol,500)

entries = chain.GetEntries()
triggered = 0
reco = 0
print(entries)

#f_mip = TFile("mip_fit.root")
#mip_cut = gDirectory.Get("f_pol")

for entry in range(entries):
    if entry % 1000 == 0: print(entry)
    ientry = chain.LoadTree(entry)
    nb = chain.GetEntry(entry)
    
    l_mucs = chain.MuCS_TPC_len
    x = chain.MuCS_Start_TPC[0]
        
    if data:
        mip = l_mucs > 0 and chain.flash_pe > 0 #and cosmic_pe/l_mucs > 6.72 - 0.08*x + 0.00035*x*x - 5.72e-07*x*x*x - 0.7
        xy_tpc = math.degrees(chain.MinD_theta_xy)
        yz_tpc = math.degrees(chain.MinD_theta_yz)
    else:
        mip = l_mucs > 0 and chain.MCflash_pe > 0 
        xy_tpc = math.degrees(chain.MCTagged_theta_xy)
        yz_tpc = math.degrees(chain.MCTagged_theta_yz)
        l_tpc = chain.MCTagged_len
        
    xy_mucs = math.degrees(chain.MuCS_theta_xy)
    yz_mucs = math.degrees(chain.MuCS_theta_yz)
 
    if chain.MuCS_NHitsX < 8 and chain.MuCS_NHitsZ < 8 and mip:
        x = chain.MuCS_Start_TPC[0]-chain.MuCS_Start[0]
        y = chain.MuCS_Start_TPC[1]-chain.MuCS_Start[1]
        z = chain.MuCS_Start_TPC[2]-chain.MuCS_Start[2]
        r = math.sqrt(x**2+y**2+z**2)
        
        theta_mucs = math.acos(z/r)
        theta_mucs = math.degrees(theta_mucs)
        
        phi_mucs = math.atan2(y,x)
        phi_mucs = math.degrees(phi_mucs)
                
        triggered += 1

        if phi_mucs < -45:
            h_theta_phi_l_mucs.Fill(theta_mucs, phi_mucs, l_mucs)
            
            if xy_tpc and (chain.MinD_dist < 44 or not data):
                reco += 1
                h_theta_phi_l_tpc.Fill(theta_mucs, phi_mucs, l_mucs)

eff = reco/triggered
err = math.sqrt((eff*(1-eff))/triggered)
print("Integrated efficiency: %.1f +- %.1f"%(round(eff*100,1), round(err*100,1)))
for i in range(1, h_theta_phi_l_tpc.GetNbinsX()+2):
    for j in range(1, h_theta_phi_l_tpc.GetNbinsY()+2):
        for k in range(1, h_theta_phi_l_tpc.GetNbinsZ()+2):
            if h_theta_phi_l_tpc.GetBinContent(i,j,k) < 5:
                h_theta_phi_l_tpc.SetBinContent(i,j,k,0)
                h_theta_phi_l_mucs.SetBinContent(i,j,k,0)

                
h_tpc = h_theta_phi_l_tpc.Clone()
h_tpc.SetName("h_tpc")
h_mucs = h_theta_phi_l_mucs.Clone()
h_mucs.SetName("h_mucs")

f = open("output/%s.txt" % position,"w")
for i in range(1, h_theta_phi_l_tpc.GetNbinsX()+2):
    for j in range(1, h_theta_phi_l_tpc.GetNbinsY()+2):
        for k in range(1, h_theta_phi_l_tpc.GetNbinsZ()+2):
            if h_theta_phi_l_mucs.GetBinContent(i,j,k) and h_theta_phi_l_tpc.GetBinContent(i,j,k):
                eff = h_theta_phi_l_tpc.GetBinContent(i,j,k)/h_theta_phi_l_mucs.GetBinContent(i,j,k)
                mucs = h_theta_phi_l_mucs.GetBinContent(i,j,k)
                error = math.sqrt((eff*(1-eff))/mucs)

                print(i,j,k,eff,error,file=f)
                h_theta_phi_l_tpc.SetBinContent(i,j,k,eff)
                h_theta_phi_l_tpc.SetBinError(i,j,k,error)

f.close()


h_theta_phi_l_tpc.GetZaxis().SetTitleOffset(1.7)
h_theta_phi_l_tpc.GetXaxis().SetTitleOffset(1.7)
h_theta_phi_l_tpc.GetYaxis().SetTitleOffset(1.7)

h_theta_phi = TH2F("h_theta_phi",";#theta [#circ];#phi [#circ]",bin_ang,0,180,bin_ang,-180,0)

h_theta_l = TH2F("h_theta_l",";#theta [#circ];L [cm]",bin_ang,0,180,bin_len,fidvol,500)

h_phi_l = TH2F("h_phi_l",";#phi [#circ];L [cm]",bin_ang,-180,0,bin_len,fidvol,500)

h_theta = TH1F("h_theta", ";#theta [#circ];N. Entries / %i#circ" % int(90/bin_ang),bin_ang,0,180)
h_phi = TH1F("h_phi", ";#phi [#circ];N. Entries / %i#circ" % int(180/bin_ang),bin_ang,-180,0)
h_l = TH1F("h_l", ";L [cm];N. Entries / %i cm" % int((500-fidvol)/bin_len), bin_len, fidvol, 500)

for i in range(1,h_theta_phi_l_tpc.GetNbinsX()+2):

    tpc = sum([h_tpc.GetBinContent(i,j,k) for j in range(1,h_tpc.GetNbinsY()+2) for k in range(1,h_tpc.GetNbinsZ()+2)])
    mucs = sum([h_mucs.GetBinContent(i,j,k) for j in range(1,h_mucs.GetNbinsY()+2) for k in range(1,h_mucs.GetNbinsZ()+2)])
    if tpc and mucs:
        eff = tpc/mucs
        error = math.sqrt((eff*(1-eff))/mucs)
        h_theta.SetBinContent(i,eff)
        h_theta.SetBinError(i,error)
        
for i in range(1,h_theta_phi_l_tpc.GetNbinsY()+2):   
    tpc = sum([h_tpc.GetBinContent(j,i,k) for j in range(1,h_tpc.GetNbinsX()+2) for k in range(1,h_tpc.GetNbinsZ()+2)])
    mucs = sum([h_mucs.GetBinContent(j,i,k) for j in range(1,h_mucs.GetNbinsX()+2) for k in range(1,h_mucs.GetNbinsZ()+2)])
    if tpc and mucs:
        eff = tpc/mucs
        error = math.sqrt((eff*(1-eff))/mucs)
        h_phi.SetBinContent(i,eff)
        h_phi.SetBinError(i,error)

for i in range(1,h_theta_phi_l_tpc.GetNbinsZ()+2):
    tpc = sum([h_tpc.GetBinContent(j,k,i) for j in range(1,h_tpc.GetNbinsX()+2) for k in range(1,h_tpc.GetNbinsY()+2)])
    mucs = sum([h_mucs.GetBinContent(j,k,i) for j in range(1,h_mucs.GetNbinsX()+2) for k in range(1,h_mucs.GetNbinsY()+2)])
    if tpc and mucs:
        eff = tpc/mucs
        error = math.sqrt((eff*(1-eff))/mucs)
        h_l.SetBinContent(i,eff)
        h_l.SetBinError(i,error)


for i in range(1,h_theta_phi_l_tpc.GetNbinsX()+2):
    for j in range(1,h_theta_phi_l_tpc.GetNbinsY()+2):
        tpc = sum([h_tpc.GetBinContent(i,j,k) for k in range(1,h_tpc.GetNbinsZ()+2)])
        mucs = sum([h_mucs.GetBinContent(i,j,k) for k in range(1,h_mucs.GetNbinsZ()+2)])
        if mucs and tpc:
            eff = tpc/mucs
            error = math.sqrt((eff*(1-eff))/mucs)
            h_theta_phi.SetBinContent(i,j,eff)
            h_theta_phi.SetBinError(i,j,error)

for i in range(1,h_theta_phi_l_tpc.GetNbinsX()+2):
    for j in range(1,h_theta_phi_l_tpc.GetNbinsZ()+2): 
        tpc = sum([h_tpc.GetBinContent(i,k,j) for k in range(1,h_tpc.GetNbinsY()+2)])
        mucs = sum([h_mucs.GetBinContent(i,k,j) for k in range(1,h_mucs.GetNbinsY()+2)])
        if tpc and mucs:
            eff = tpc/mucs
            error = math.sqrt((eff*(1-eff))/mucs)
            h_theta_l.SetBinContent(i,j,eff)
            h_theta_l.SetBinError(i,j,error)     

for i in range(1,h_theta_phi_l_tpc.GetNbinsY()+2):
    for j in range(1,h_theta_phi_l_tpc.GetNbinsZ()+2):
        tpc = sum([h_tpc.GetBinContent(k,i,j) for k in range(1,h_tpc.GetNbinsX()+2)])
        mucs = sum([h_mucs.GetBinContent(k,i,j) for k in range(1,h_mucs.GetNbinsX()+2)])
        if tpc and mucs:
            eff = tpc/mucs
            error = math.sqrt((eff*(1-eff))/mucs)
            h_phi_l.SetBinContent(i,j,eff)
            h_phi_l.SetBinError(i,j,error)           

c_theta_phi = TCanvas("c_theta_phi","theta_phi")
c_theta_phi.cd()

h_theta_phi.Draw("colz texte")
h_theta_phi.SetLineWidth(1)
h_theta_phi.SetLineColor(kGray+2)
h_theta_phi.GetZaxis().SetTitle("Efficiency")
h_theta_phi.GetZaxis().RotateTitle()
gPad.SetRightMargin(0.15)
h_theta_phi.GetZaxis().SetRangeUser(0,1)
pt = TPaveText(0.15,0.77,0.45,0.84, "ndc")
pt.AddText("MicroBooNE in progress")
pt.SetFillColor(0)
pt.SetBorderSize(0)
pt.SetShadowColor(0)
pt.Draw()
h_theta_phi.SaveAs("plots/%s/e_theta_phi_%s.root" % ("data" if data else "mc",algo))

c_theta_phi.Update()

c_theta_l = TCanvas("c_theta_l","theta_l")
c_theta_l.cd()

h_theta_l.Draw("col texte")
h_theta_l.GetZaxis().SetRangeUser(0,1)
h_theta_l.SaveAs("plots/%s/e_theta_l_%s.root" % ("data" if data else "mc",algo))

c_theta_l.Update()


c_phi_l = TCanvas("c_phi_l","phi_l")
c_phi_l.cd()

h_phi_l.Draw("col texte")
h_phi_l.GetZaxis().SetRangeUser(0,1)
h_phi_l.SaveAs("plots/%s/e_phi_l_%s.root" % ("data" if data else "mc", algo))

c_phi_l.Update()


c_theta = TCanvas("c_theta","theta")
c_theta.cd()

h_theta.Draw()
h_theta.GetYaxis().SetRangeUser(0,1)
h_theta.SetMarkerStyle(20)
h_theta.SaveAs("plots/%s/e_theta_%s.root" % ("data" if data else "mc", algo))

c_theta.Update()

c_phi = TCanvas("c_phi","phi")
c_phi.cd()

h_phi.Draw()
h_phi.GetYaxis().SetRangeUser(0,1)
h_phi.SetMarkerStyle(20)
h_phi.SaveAs("plots/%s/e_phi_%s.root" % ("data" if data else "mc", algo))

c_phi.Update()

c_l = TCanvas("c_l","l")
c_l.cd()

h_l.Draw("p")
h_l.GetYaxis().SetRangeUser(0,1)
h_l.SetMarkerStyle(20)
h_l.SaveAs("plots/%s/e_l_%s.root" % ("data" if data else "mc", algo))

c_l.Update()

gStyle.SetCanvasPreferGL(1)
c_theta_phi_l = TCanvas("c_theta_phi_l","theta_phi_l",500,500)

h_theta_phi_l_tpc.Draw("glbox")
h_theta_phi_l_tpc.GetZaxis().SetRangeUser(20,370)
h_theta_phi_l_tpc.GetXaxis().SetNdivisions(5,0,5)

h_theta_phi_l_tpc.SaveAs("plots/%s/e_theta_phi_l_%s.root" % ("data" if data else "mc", algo))

c_theta_phi_l.Update()


input()