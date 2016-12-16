#!/usr/bin/env python3.4

from ROOT import TPad, gPad, gDirectory, TFile, gStyle, TChain, TH1F, TH2F, TH3F, TLine, TCanvas, TPad, kRed, kGray, kBlue, TBox, TPaveText
import math
import sys
from glob import glob
from geomUtil import isect_line_plane_v3

if len(sys.argv) > 1:
    position = sys.argv[1]
else:
    position = "all"

print(position, "dataset")

bin_ang = 12
fidvol = 20
bin_len = 8

x_start = 0
x_end = 256.35
y_start = -116.5
y_end = 116.5
z_start = 0
z_end = 1036.8


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

if position == "central_shift":
    f = glob("../root_files/mucs_cry_overlay/*/ana_hist.root")
elif position == "downstream_shift":
    f = glob("../root_files/mucs_cry_overlay_downstream/*/ana_hist.root")
elif position == "upstream_shift":
    f = glob("../root_files/mucs_cry_overlay_upstream/*/ana_hist.root")
elif position == "merged":
    f = glob("../root_files/mucs_cry_overlay*/*/ana_hist.root")
elif position == "all":
    f = ["../root_files/mcc7_ana.root"]

chain = TChain("analysistree/anatree")
for filename in f:
    chain.Add(filename)

entries = chain.GetEntries()


h_theta_phi_geant = TH2F("h_theta_phi_geant",";#theta [#circ];#phi [#circ]",bin_ang,0,180,bin_ang,-180,0)
h_theta_phi_reco = TH2F("h_theta_phi_reco",";#theta [#circ];#phi [#circ]",bin_ang,0,180,bin_ang,-180,0)

h_theta_l_geant = TH2F("h_theta_l_geant",";#theta [#circ];L [cm]",bin_ang,0,180,bin_len,fidvol,500)
h_theta_l_reco = TH2F("h_theta_l_reco",";#theta [#circ];L [cm]",bin_ang,0,180,bin_len,fidvol,500)

h_phi_l_geant = TH2F("h_phi_l_geant",";#phi [#circ];L [cm]",bin_ang,-180,0,bin_len,fidvol,500)
h_phi_l_reco = TH2F("h_phi_l_reco",";#phi [#circ];L [cm]",bin_ang,-180,0,bin_len,fidvol,500)

h_theta_reco = TH1F("h_theta_reco",";#theta [#circ]; N. Entries / 15#circ",bin_ang,0,180)
h_theta_geant = TH1F("h_theta_geant",";#theta [#circ]; N. Entries / 15#circ",bin_ang,0,180)

h_phi_reco = TH1F("h_phi_reco",";#phi [#circ]; N. Entries / 15#circ",bin_ang,-180,0)
h_phi_geant = TH1F("h_phi_geant",";#phi [#circ]; N. Entries / 15#circ",bin_ang,-180,0)

h_l_geant = TH1F("h_l_geant",";L [cm]; N. Entries / 80 cm",bin_len,fidvol,500)
h_l_reco = TH1F("h_l_reco",";L [cm]; N. Entries / 80 cm",bin_len,fidvol,500)

h_theta_phi_l_reco = TH3F("h_theta_phi_l_reco",";#theta [#circ]; #phi [#circ]; L [cm]",bin_ang,0,180,bin_ang,-180,0,bin_len,fidvol,500)
h_theta_phi_l_geant = TH3F("h_theta_phi_l_geant",";#theta [#circ]; #phi [#circ]; L [cm]",bin_ang,0,180,bin_ang,-180,0,bin_len,fidvol,500)

h_dist = TH1F("h_dist",";Distance [cm]; N. Entries / 2 cm", 25, 0, 50)

anode_plane = [0,0,0]
anode_no = [1,0,0]
cathode_plane = [x_end,0,0]
cathode_no = anode_no
top_plane = [0,y_end,0]
top_no = [0,1,0]
bottom_plane = [0,y_start,0]
bottom_no = top_no
upstream_plane = [0,0,z_start]
upstream_no = [0,0,1]
downstream_plane = [0,0,z_end]
downstream_no = upstream_no

print(entries)
for entry in range(1000):

    if entry % 10 == 0: print(entry)

    ientry = chain.LoadTree(entry)
    nb = chain.GetEntry(entry)

    primaries = 0
    reco = 0

    for i in range(chain.geant_list_size):
        x_endpoint = chain.EndPointx[i]
        y_endpoint = chain.EndPointy[i]
        z_endpoint = chain.EndPointz[i]
        mip = (x_endpoint < x_start or x_endpoint > x_end) and (y_endpoint < y_start or y_endpoint > y_end) and (z_endpoint < z_start or z_endpoint > z_end)

        if position != "all":
            mucs_event = chain.StartPointy[i] > 390 and chain.StartPointy[i] < 399
        else:
            mucs_event = True

        if chain.process_primary[i] == 1 and chain.inTPCDrifted[i] == 1 and mucs_event:
            x = chain.Px[i]
            y = chain.Py[i]
            z = chain.Pz[i]
            r = chain.P[i]
            theta = math.degrees(math.acos(z/r))
            phi = math.degrees(math.atan2(y,x))
            l = (chain.StartPointx_tpcAV[i]-chain.EndPointx_tpcAV[i])**2+(chain.StartPointy_tpcAV[i]-chain.EndPointy_tpcAV[i])**2+(chain.StartPointz_tpcAV[i]-chain.EndPointz_tpcAV[i])**2
            l = math.sqrt(l)
            h_theta_phi_l_geant.Fill(theta, phi, l)

            p0_geant = [chain.StartPointx_drifted[i], chain.StartPointy_drifted[i], chain.StartPointz_drifted[i]]



            for n in range(chain.ntracks_pandoraCosmic):
                best = chain.trkpidbestplane_pandoraCosmic[n]
                trkid = chain.trkidtruth_pandoraCosmic[3*n+best]

                x1 = chain.trkstartx_pandoraCosmic[n]
                y1 = chain.trkstarty_pandoraCosmic[n]
                z1 = chain.trkstartz_pandoraCosmic[n]
                # x_dir = chain.trkstartdcosx_pandoraCosmic[n]
                # y_dir = chain.trkstartdcosy_pandoraCosmic[n]
                # z_dir = chain.trkstartdcosz_pandoraCosmic[n]
                # p_no_tpc = [0,1,0]
                # p_tpc = [0,y_end,0]
                #
                #

                p0 = [x1,y1,z1]
                #p1 = [x1+x_dir,y1+y_dir,z1+z_dir]
                # intersect_tpc = isect_line_plane_v3(p0,p1,p_tpc,p_no_tpc)
                # anode_p = isect_line_plane_v3(p0, p1, anode_plane, anode_no)
                # cathode_p = isect_line_plane_v3(p0, p1, cathode_plane, cathode_no)
                # upstream_p = isect_line_plane_v3(p0, p1, upstream_plane, upstream_no)
                # downstream_p = isect_line_plane_v3(p0, p1, downstream_plane, downstream_no)
                # top_p = isect_line_plane_v3(p0, p1, top_plane, top_no)
                # bottom_p = isect_line_plane_v3(p0, p1, bottom_plane, bottom_no)
                #
                # min_dist = 9999
                # closest_point = []
                # points = [anode_p, cathode_p, upstream_p, downstream_p, top_p, bottom_p]
                # for point in points:
                #
                #     if point:
                #         dist = math.sqrt(sum([(a-b)**2 for a,b in zip(p0,point)]))
                #         #print(point[0],point[1],point[2])
                #         if point[0] >= x_start and point[0] <= x_end and point[1] >= y_start and point[1] <= y_end and point[2] >= z_start and point[2] <= z_end:
                #             min_dist = dist
                #             closest_point = point
                #
                dist = math.sqrt(sum([(a-b)**2 for a,b in zip(p0,p0_geant)]))
                dist2 = math.sqrt((p0_geant[0]-p0[0])**2+(p0_geant[1]-p0[1])**2)
                if trkid == chain.TrackId[i]:
                    h_dist.Fill(dist)
                    if dist2 < 32:
                        h_theta_phi_l_reco.Fill(theta, phi, l)
                        break


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

overall_eff = h_reco.Integral()/h_geant.Integral()
overall_eff_err = math.sqrt((overall_eff*(1-overall_eff))/h_geant.Integral())
print("Overall efficiency: %.1f%% +- %.1f%%" % (round(overall_eff*100,1),round(overall_eff_err*100,1)))

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

pt = TPaveText(0.10,0.905,0.40,0.98, "ndc")
pt.AddText("MicroBooNE in progress")
pt.SetFillColor(0)
pt.SetBorderSize(0)
pt.SetShadowColor(0)
h_theta_phi_reco.GetXaxis().SetRangeUser(60,120)
h_theta_phi_reco.GetYaxis().SetRangeUser(-90,-45)
h_theta_phi_reco.GetZaxis().SetRangeUser(0,1)

h_theta_phi_reco.SetMarkerSize(2)
h_theta_phi_reco.Draw("colz texte")
pt.Draw()
h_theta_phi_reco.SaveAs("plots/mc/theta_phi_mcc7.root")

c_theta_phi.Update()
c_theta_phi.SaveAs("plots/mc/theta_phi_mc.pdf")

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

h_theta_l_reco.GetXaxis().SetRangeUser(60,120)
h_theta_l_reco.GetYaxis().SetRangeUser(20,320)
h_theta_l_reco.GetZaxis().SetRangeUser(0,1)
h_theta_l_reco.SetMarkerSize(2)
h_theta_l_reco.Draw("colz texte")
pt.Draw()
h_theta_l_reco.SaveAs("plots/mc/theta_l_mcc7.root")
c_theta_l.Update()
c_theta_l.SaveAs("plots/mc/theta_l_mc.pdf")

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

h_phi_l_reco.GetXaxis().SetRangeUser(-90,-45)
h_phi_l_reco.GetYaxis().SetRangeUser(20,320)
h_phi_l_reco.GetZaxis().SetRangeUser(0,1)
h_phi_l_reco.SetMarkerSize(2)
h_phi_l_reco.Draw("colz texte")
pt.Draw()
h_phi_l_reco.SaveAs("plots/mc/phi_l_mcc7.root")
c_phi_l.Update()
c_phi_l.SaveAs("plots/mc/phi_l_mc.pdf")


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

c_dist = TCanvas("c_dist")
h_dist.Draw()
c_dist.Update()

gStyle.SetCanvasPreferGL(1)

c_3d = TCanvas("c_3d","3d")
h_theta_phi_l_reco.Draw("glbox")
h_theta_phi_l_reco.GetZaxis().SetTitleOffset(1.7)
h_theta_phi_l_reco.GetXaxis().SetTitleOffset(1.7)
h_theta_phi_l_reco.GetYaxis().SetTitleOffset(1.7)
h_theta_phi_l_reco.GetYaxis().SetNdivisions(505)
h_theta_phi_l_reco.GetXaxis().SetRangeUser(60,120)
h_theta_phi_l_reco.GetYaxis().SetRangeUser(-90,-45)
h_theta_phi_l_reco.GetZaxis().SetRangeUser(20,320)

c_3d.Update()
input()
