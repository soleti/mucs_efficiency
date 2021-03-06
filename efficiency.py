#!/usr/bin/env python3.4

import os, sys
import math
from ROOT import TEfficiency
from ROOT import TPad, gDirectory, TChain, TFile, TH1F, TH2F, TH3F, gStyle, TCanvas, gPad, TPaveText, kGray, TLine, kRed, TF1, TLegend
from geomUtil import isect_line_plane_v3

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
gStyle.SetPalette(57)
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
    chain.Add("../root_files/MuCSRun7263_Group180_MergedTree.Root")
    chain.Add("../root_files/MuCSRun7264_Group180_MergedTree.Root")
    chain.Add("../root_files/MuCSRun7265_Group180_MergedTree.Root")
if position == "downstream_shift":
    chain.Add("../root_files/MuCSRun7347_Group181_MergedTree.Root")
    chain.Add("../root_files/MuCSRun7348_Group181_MergedTree.Root")
if position == "central_shift":
    chain.Add("../root_files/CORSIKA.root")
if position == "upstream_shift":
    chain.Add("../root_files/MuCSRun7702_Group182_MergedTree.Root")
    chain.Add("../root_files/MuCSRun7703_Group183_MergedTree.Root")
if position == "all":
    chain.Add("../root_files/CORSIKA.root")
    chain.Add("../root_files/MuCSRun7347_Group181_MergedTree.Root")
    chain.Add("../root_files/MuCSRun7348_Group181_MergedTree.Root")
    chain.Add("../root_files/MuCSRun7702_Group182_MergedTree.Root")
    chain.Add("../root_files/MuCSRun7703_Group183_MergedTree.Root")

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
h_dist = TH1F("h_dist",";Distance [cm]; N. Entries / 1 cm", 60, 0, 60)
entries = chain.GetEntries()

print(entries)
wrong,right=[0]*5,[0]*5
for entry in range(entries):
    if entry % 1000 == 0: print(entry)
    ientry = chain.LoadTree(entry)
    nb = chain.GetEntry(entry)

    l_mucs = chain.MuCS_TPC_len
    x = chain.MuCS_Start_TPC[0]
    cosmic_pe = chain.flash_pe

    # MIP fit
    a = 7.71739e+00
    b = -1.06721e-01
    c = 5.82904e-04
    d = -1.17978e-06
    #print(a + b*(x-5) + c*(x-5)**2 +d*(x-5)**3)
    mip = l_mucs > 0 and cosmic_pe > 0# and cosmic_pe/l_mucs > a + b*(x-5) + c*(x-5)**2 +d*(x-5)**3

    if chain.MuCS_NHitsX < 8 and chain.MuCS_NHitsZ < 8 and mip and chain.evt_number != 17453:
        x = chain.MuCS_Start_TPC[0]-chain.MuCS_Start[0]
        y = chain.MuCS_Start_TPC[1]-chain.MuCS_Start[1]
        z = chain.MuCS_Start_TPC[2]-chain.MuCS_Start[2]
        r = math.sqrt(x**2+y**2+z**2)
        pe = chain.flash_pe
        theta_mucs = math.acos(z/r)
        theta_mucs = math.degrees(theta_mucs)

        phi_mucs = math.atan2(y,x)
        phi_mucs = math.degrees(phi_mucs)

        if chain.MinD_len < 320 and phi_mucs < -45:
            p_no_tpc = [0,1,0]
            p_tpc = [0,y_end,0]
            p0 = [chain.MinD_Start[0], chain.MinD_Start[1], chain.MinD_Start[2]]
            p1 = [chain.MinD_End[0], chain.MinD_End[1], chain.MinD_End[2]]
            intersect_tpc = isect_line_plane_v3(p0,p1,p_tpc,p_no_tpc)

            length = chain.MinD_len

            through_x = intersect_tpc[0] > x_start and intersect_tpc[0] < x_end
            through_z = intersect_tpc[2] > z_start and intersect_tpc[2] < z_end
            through_tpc = through_x and through_z

            if through_tpc:
                h_theta_phi_l_mucs.Fill(theta_mucs, phi_mucs, length)
                dist = math.sqrt((chain.MuCS_Start_TPC[0]-chain.MinD_Start[0])**2+(chain.MuCS_Start_TPC[1]-chain.MinD_Start[1])**2)
                dist2 = math.sqrt((chain.MuCS_Start_TPC[0]-chain.MinD_Start[0])**2+(chain.MuCS_Start_TPC[2]-chain.MinD_Start[2])**2)
                dist3 = math.sqrt((chain.MuCS_Start_TPC[0]-chain.MinD_Start[0])**2+(chain.MuCS_Start_TPC[1]-chain.MinD_Start[1])**2+(chain.MuCS_Start_TPC[2]-chain.MinD_Start[2])**2)

                bin_mucs = int((chain.MuCS_TPC_len-20)/60)
                bin_reco = int((chain.MinD_len-20)/60)

                h_dist.Fill(dist2)
                if dist < 32:
                    h_theta_phi_l_tpc.Fill(theta_mucs, phi_mucs, length)

print("Number of events: ", h_theta_phi_l_mucs.Integral())
eff = h_theta_phi_l_tpc.Integral()/(h_theta_phi_l_mucs.Integral())
err = math.sqrt((eff*(1-eff))/h_theta_phi_l_mucs.Integral())
print("Integrated efficiency: %.1f +- %.1f" % (round(eff*100,1), round(err*100,1)))
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
sys_errors = open("output/sys_errors_3d.txt","r")
sys_errors_bins=sys_errors.readlines()

for i in range(1, h_theta_phi_l_tpc.GetNbinsX()+2):
    for j in range(1, h_theta_phi_l_tpc.GetNbinsY()+2):
        for k in range(1, h_theta_phi_l_tpc.GetNbinsZ()+2):
            bin1 = i == 6 and j == 7 and k == 5
            bin2 = i == 7 and j == 7 and k == 5

            if h_theta_phi_l_mucs.GetBinContent(i,j,k) and h_theta_phi_l_tpc.GetBinContent(i,j,k) and not bin1 and not bin2:

                eff = h_theta_phi_l_tpc.GetBinContent(i,j,k)/h_theta_phi_l_mucs.GetBinContent(i,j,k)
                mucs = h_theta_phi_l_mucs.GetBinContent(i,j,k)
                error = math.sqrt((eff*(1-eff))/mucs)

                sys_error = float(sys_errors_bins[(i-1)*13*9+(j-1)*9+(k-1)].split()[3])
                if sys_error < 0: sys_error = 0

                print(i,j,k,eff,error,file=f)

                h_theta_phi_l_tpc.SetBinContent(i,j,k,eff)
                h_theta_phi_l_tpc.SetBinError(i,j,k,error)

            else:
                h_theta_phi_l_tpc.SetBinContent(i,j,k,0)
                h_theta_phi_l_tpc.SetBinError(i,j,k,0)
                print(i,j,k,0,0,file=f)

sys_errors.close()
f.close()


h_theta_phi_l_tpc.GetZaxis().SetTitleOffset(1.7)
h_theta_phi_l_tpc.GetXaxis().SetTitleOffset(1.7)
h_theta_phi_l_tpc.GetYaxis().SetTitleOffset(1.7)

h_theta_phi = TH2F("h_theta_phi",";#theta [#circ];#phi [#circ]",bin_ang,0,180,bin_ang,-180,0)
h_theta_l = TH2F("h_theta_l",";#theta [#circ];L [cm]",bin_ang,0,180,bin_len,fidvol,500)
h_phi_l = TH2F("h_phi_l",";#phi [#circ];L [cm]",bin_ang,-180,0,bin_len,fidvol,500)

h_theta = TH1F("h_theta", ";#theta [#circ];N. Entries / %i#circ" % int(180/bin_ang),bin_ang,0,180)
h_phi = TH1F("h_phi", ";#phi [#circ];N. Entries / %i#circ" % int(180/bin_ang),bin_ang,-180,0)
h_l = TH1F("h_l", ";L [cm];N. Entries / %i cm" % int((500-fidvol)/bin_len), bin_len, fidvol, 500)
h_l_sys = TH1F("h_l_sys", ";L [cm];N. Entries / %i cm" % int((500-fidvol)/bin_len), bin_len, fidvol, 500)
h_theta_sys = TH1F("h_theta_sys", ";#theta [#circ];N. Entries / %i#circ" % int(180/bin_ang),bin_ang,0,180)
h_phi_sys = TH1F("h_phi_sys", ";#phi [#circ];N. Entries / %i#circ" % int(180/bin_ang),bin_ang,-180,0)

f_theta = open("output/theta_%s.txt" % position,"w")
f_theta_sys = open("output/sys_errors_theta.txt","r")
theta_sys=f_theta_sys.readlines()
h_theta_tpc = TH1F("h_theta_tpc","",bin_ang,0,180)
h_theta_mucs = TH1F("h_theta_mucs","",bin_ang,0,180)

dif_corr = 0.01
dif_corr_err = 0.001

for i in range(1,h_theta_phi_l_tpc.GetNbinsX()+2):

    tpc = sum([h_tpc.GetBinContent(i,j,k) for j in range(1,h_tpc.GetNbinsY()+2) for k in range(1,h_tpc.GetNbinsZ()+2)])
    mucs = sum([h_mucs.GetBinContent(i,j,k) for j in range(1,h_mucs.GetNbinsY()+2) for k in range(1,h_mucs.GetNbinsZ()+2)])

    if tpc and mucs:
        h_theta_tpc.SetBinContent(i,tpc)
        h_theta_mucs.SetBinContent(i,mucs)
        eff = tpc/mucs+dif_corr
        error = math.sqrt((eff*(1-eff))/mucs)
        sys_error = abs(float(theta_sys[i-1].split()[1]))+dif_corr_err

        print(i,eff,error,file=f_theta)
        h_theta.SetBinContent(i,eff)
        h_theta.SetBinError(i,error)
        h_theta_sys.SetBinContent(i,eff)
        h_theta_sys.SetBinError(i,math.sqrt(sys_error**2+error**2))
    else:
        print(i,0,0,file=f_theta)


f_theta.close()
f_theta_sys.close()

f_phi = open("output/phi_%s.txt" % position,"w")
f_phi_sys = open("output/sys_errors_phi.txt","r")
phi_sys = f_phi_sys.readlines()

h_phi_mucs = TH1F("h_phi_mucs","",bin_ang,-180,0)

for i in range(1,h_theta_phi_l_tpc.GetNbinsY()+2):
    tpc = sum([h_tpc.GetBinContent(j,i,k) for j in range(1,h_tpc.GetNbinsX()+2) for k in range(1,h_tpc.GetNbinsZ()+2)])
    mucs = sum([h_mucs.GetBinContent(j,i,k) for j in range(1,h_mucs.GetNbinsX()+2) for k in range(1,h_mucs.GetNbinsZ()+2)])

    if tpc and mucs:
        eff = tpc/mucs+dif_corr
        error = math.sqrt((eff*(1-eff))/mucs)
        h_phi_mucs.SetBinContent(i,mucs)
        print(i,eff,error,file=f_phi)

        sys_error = abs(float(phi_sys[i-1].split()[1]))+dif_corr_err

        h_phi.SetBinContent(i,eff)
        h_phi.SetBinError(i,error)
        h_phi_sys.SetBinContent(i,eff)
        h_phi_sys.SetBinError(i,math.sqrt(sys_error**2+error**2))
    else:
        print(i,0,0,file=f_phi)

f_l = open("output/l_%s.txt" % position,"w")
f_l_sys = open("output/sys_errors_l.txt","r")
l_sys = f_l_sys.readlines()

for i in range(1,h_theta_phi_l_tpc.GetNbinsZ()+2):
    tpc = sum([h_tpc.GetBinContent(j,k,i) for j in range(1,h_tpc.GetNbinsX()+2) for k in range(1,h_tpc.GetNbinsY()+2)])
    mucs = sum([h_mucs.GetBinContent(j,k,i) for j in range(1,h_mucs.GetNbinsX()+2) for k in range(1,h_mucs.GetNbinsY()+2)])
    if tpc and mucs:
        eff = tpc/mucs+dif_corr
        error = math.sqrt((eff*(1-eff))/mucs)

        print(i,eff,error,file=f_l)
        sys_error = abs(float(l_sys[i-1].split()[1]))+dif_corr_err

        h_l.SetBinContent(i,eff)
        h_l.SetBinError(i,error)
        h_l_sys.SetBinContent(i,eff)
        h_l_sys.SetBinError(i,math.sqrt(sys_error**2+error**2))
    else:
        print(i,0,0,file=f_l)


f_l.close()
f_l_sys.close()


f_theta_phi = open("output/theta_phi_%s.txt" % position,"w")
f_theta_phi_sys = open("output/sys_errors_theta_phi.txt","r")
theta_phi_sys = f_theta_phi_sys.readlines()

for i in range(1,h_theta_phi_l_tpc.GetNbinsX()+2):
    for j in range(1,h_theta_phi_l_tpc.GetNbinsY()+2):
        tpc = sum([h_tpc.GetBinContent(i,j,k) for k in range(1,h_tpc.GetNbinsZ()+2)])
        mucs = sum([h_mucs.GetBinContent(i,j,k) for k in range(1,h_mucs.GetNbinsZ()+2)])
        if mucs and tpc:
            eff = tpc/mucs+dif_corr
            #print(eff)
            if eff<1:
                error = math.sqrt((eff*(1-eff))/mucs)

            print(i,eff,error,file=f_theta_phi)

            sys_error = float(theta_phi_sys[(i-1)*13+(j-1)].split()[1])+dif_corr_err

            h_theta_phi.SetBinContent(i,j,eff)
            h_theta_phi.SetBinError(i,j,math.sqrt(sys_error**2+error**2))
        else:
            print(i,0,0,file=f_theta_phi)

f_theta_l = open("output/theta_l_%s.txt" % position,"w")
f_theta_l_sys = open("output/sys_errors_theta_l.txt","r")
theta_l_sys = f_theta_l_sys.readlines()

for i in range(1,h_theta_phi_l_tpc.GetNbinsX()+2):
    for j in range(1,h_theta_phi_l_tpc.GetNbinsZ()+2):
        tpc = sum([h_tpc.GetBinContent(i,k,j) for k in range(1,h_tpc.GetNbinsY()+2)])
        mucs = sum([h_mucs.GetBinContent(i,k,j) for k in range(1,h_mucs.GetNbinsY()+2)])
        if tpc and mucs:
            eff = tpc/mucs+dif_corr
            error = math.sqrt((eff*(1-eff))/mucs)

            print(i,eff,error,file=f_theta_l)

            sys_error = float(theta_l_sys[(i-1)*9+(j-1)].split()[1])+dif_corr_err
            if sys_error < 0: sys_error = 0

            h_theta_l.SetBinContent(i,j,eff)
            h_theta_l.SetBinError(i,j,math.sqrt(sys_error**2+error**2))
        else:
            print(i,0,0,file=f_theta_l)


f_phi_l = open("output/phi_l_%s.txt" % position,"w")
f_phi_l_sys = open("output/sys_errors_phi_l.txt","r")
phi_l_sys = f_phi_l_sys.readlines()

for i in range(1,h_theta_phi_l_tpc.GetNbinsY()+2):
    for j in range(1,h_theta_phi_l_tpc.GetNbinsZ()+2):
        tpc = sum([h_tpc.GetBinContent(k,i,j) for k in range(1,h_tpc.GetNbinsX()+2)])
        mucs = sum([h_mucs.GetBinContent(k,i,j) for k in range(1,h_mucs.GetNbinsX()+2)])
        if tpc and mucs:
            eff = tpc/mucs+dif_corr
            if eff<1:
                error = math.sqrt((eff*(1-eff))/mucs)

            print(i,eff,error,file=f_phi_l)
            sys_error = float(phi_l_sys[(i-1)*9+(j-1)].split()[1])+dif_corr_err
            if sys_error < 0: sys_error = 0

            h_phi_l.SetBinContent(i,j,eff)
            h_phi_l.SetBinError(i,j,math.sqrt(sys_error**2+error**2))
        else:
            print(i,0,0,file=f_phi_l)


c_theta_phi = TCanvas("c_theta_phi","theta_phi")
c_theta_phi.cd()

h_theta_phi.Draw("colz texte")
h_theta_phi.SetLineWidth(1)
h_theta_phi.SetLineColor(kGray+2)
#h_theta_phi.GetZaxis().SetTitle("Efficiency")
#h_theta_phi.GetZaxis().RotateTitle()
#gPad.SetRightMargin(0.15)
h_theta_phi.GetZaxis().SetRangeUser(0,1)

h_theta_phi.GetXaxis().SetRangeUser(60,120)
h_theta_phi.GetYaxis().SetRangeUser(-90,-45)
h_theta_phi.SetMarkerSize(2.5)
gPad.SetBottomMargin(0.17)
gPad.SetLeftMargin(0.13)
gPad.SetTopMargin(0.15)
gPad.SetRightMargin(0.15)

h_theta_phi.GetYaxis().SetTitleSize(0.07)
h_theta_phi.GetXaxis().SetTitleSize(0.07)
h_theta_phi.GetYaxis().SetTitleOffset(0.8)

pt = TPaveText(0.13,0.855,0.42,0.98, "ndc")
pt.AddText("MicroBooNE")
pt.SetFillColor(0)
pt.SetBorderSize(0)
pt.SetShadowColor(0)
h_theta_phi.Draw("colz texte")

pt.Draw()
h_theta_phi.SaveAs("plots/%s/e_theta_phi_%s.root" % ("data" if data else "mc",algo))
c_theta_phi.SaveAs("plots/%s/e_theta_phi.pdf" % ("data" if data else "mc"))

c_theta_phi.Update()

c_theta_l = TCanvas("c_theta_l","theta_l")
c_theta_l.cd()

h_theta_l.Draw("colz texte")
h_theta_l.GetZaxis().SetRangeUser(0,1)
h_theta_l.GetXaxis().SetRangeUser(60,120)
h_theta_l.GetYaxis().SetRangeUser(20,320)
h_theta_l.SetMarkerSize(2.5)
gPad.SetBottomMargin(0.17)
gPad.SetLeftMargin(0.13)
gPad.SetTopMargin(0.15)
gPad.SetRightMargin(0.15)

h_theta_l.GetYaxis().SetTitleSize(0.07)
h_theta_l.GetXaxis().SetTitleSize(0.07)
h_theta_l.GetYaxis().SetTitleOffset(0.8)

h_theta_l.Draw("colz texte")

pt.Draw()
h_theta_l.SaveAs("plots/%s/e_theta_l_%s.root" % ("data" if data else "mc",algo))
c_theta_l.SaveAs("plots/%s/e_theta_l.pdf" % ("data" if data else "mc"))


c_theta_l.Update()


c_phi_l = TCanvas("c_phi_l","phi_l")
c_phi_l.cd()

h_phi_l.Draw("colz texte")
h_phi_l.GetZaxis().SetRangeUser(0,1)
h_phi_l.GetXaxis().SetRangeUser(-90,-45)
h_phi_l.GetYaxis().SetRangeUser(20,320)
h_phi_l.SetMarkerSize(2.5)
gPad.SetBottomMargin(0.17)
gPad.SetLeftMargin(0.13)
gPad.SetTopMargin(0.15)
gPad.SetRightMargin(0.15)

h_phi_l.GetYaxis().SetTitleSize(0.07)
h_phi_l.GetXaxis().SetTitleSize(0.07)
h_phi_l.GetYaxis().SetTitleOffset(0.8)

h_phi_l.Draw("colz texte")

pt.Draw()
h_phi_l.SaveAs("plots/%s/e_phi_l_%s.root" % ("data" if data else "mc", algo))
c_phi_l.SaveAs("plots/%s/e_phi_l.pdf" % ("data" if data else "mc"))


c_phi_l.Update()


c_theta = TCanvas("c_theta","theta")
c_theta.cd()

h_theta_sys.Draw()
h_theta_sys.GetYaxis().SetRangeUser(0,1)
h_theta.SetMarkerStyle(20)
h_theta.SaveAs("plots/%s/e_theta_%s.root" % ("data" if data else "mc", algo))
h_theta_sys.SaveAs("plots/%s/e_theta_sys_%s.root" % ("data" if data else "mc", algo))
c_theta.Update()

c_phi = TCanvas("c_phi","phi")
c_phi.cd()

h_phi_sys.Draw()
h_phi_sys.GetYaxis().SetRangeUser(0,1)
h_phi.SetMarkerStyle(20)
h_phi.SaveAs("plots/%s/e_phi_%s.root" % ("data" if data else "mc", algo))
h_phi_sys.SaveAs("plots/%s/e_phi_sys_%s.root" % ("data" if data else "mc", algo))
c_phi.Update()

c_dist = TCanvas("c_dist")
h_dist_mc = TH1F("h_dist_mc","",60,0,60)
f = TF1("f","landau",0,60)
f.SetParameters(10940,3.8,1.5)
h_dist_mc.FillRandom("f",int(h_dist.GetEntries()))
h_dist.Draw("ep")
gPad.SetLeftMargin(0.13)
h_dist.SetLineColor(1)
h_dist_mc.SetLineColor(kRed+1)
h_dist_mc.Draw("same")
h_dist.SetMarkerStyle(20)
h_dist.GetYaxis().SetRangeUser(0.001,2499)
limit = TLine(32,0.001,32,2499)
limit.SetLineStyle(2)
limit.SetLineWidth(2)
limit.Draw()
leg = TLegend(0.61,0.75,0.84,0.87)
leg.AddEntry(h_dist,"Data","ep")
leg.AddEntry(h_dist_mc,"MuCS Monte Carlo","l")
leg.Draw()
pt = TPaveText(0.1,0.91,0.39,0.98, "ndc")
pt.AddText("MicroBooNE")
pt.SetFillColor(0)
pt.SetBorderSize(0)
pt.SetShadowColor(0)
pt.Draw()
c_dist.Update()

gStyle.SetCanvasPreferGL(1)
c_theta_phi_l = TCanvas("c_theta_phi_l","theta_phi_l")
h_theta_phi_l_tpc.GetZaxis().SetTitleOffset(1.7)
h_theta_phi_l_tpc.GetXaxis().SetTitleOffset(1.7)
h_theta_phi_l_tpc.GetYaxis().SetTitleOffset(1.7)
h_theta_phi_l_tpc.GetYaxis().SetNdivisions(505)
h_theta_phi_l_tpc.GetXaxis().SetRangeUser(60,120)
h_theta_phi_l_tpc.GetYaxis().SetRangeUser(-90,-45)
h_theta_phi_l_tpc.GetZaxis().SetRangeUser(20,320)
h_theta_phi_l_tpc.Draw("glbox")
h_theta_phi_l_tpc.SaveAs("plots/%s/e_theta_phi_l_%s.root" % ("data" if data else "mc", algo))

c_theta_phi_l.Update()

f_data = TFile("f_data.root","RECREATE")
h_theta_mucs.Write()
h_phi_mucs.Write()
h_mucs.Write()
f_data.Close()
input()
