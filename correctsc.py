#!/usr/bin/env python3.4

import os
from ROOT import TPaveText,TPolyLine3D,TF2,TGraphErrors,TBox,TLine,TProfile, TF1,TH1F, gStyle, gDirectory, TFile, TH2F, TCanvas, THStack, kBlack, kGreen, kRed,kBlue, TLegend, kAzure, kOrange, kMagenta
from array import array
import math
import sys
import random

if len(sys.argv) > 1:
    filename = str(sys.argv[1])
else:
    filename = "/Users/soleti/uboone/root_files/CORSIKA.root"

if len(sys.argv) > 2:
    output = str(sys.argv[2])
else:
    output = "plots/spacecharge"

if len(sys.argv) > 3:
    algo = str(sys.argv[3])
else:
    algo = "pandoraCosmic"

x_start = 0
x_end = 256.35
y_start = -116.5
y_end = 116.5
z_start = 0
z_end = 1036.8

def correct_theta_xy(theta_xy, start_point, f_correction):
    angle = theta_xy
    second_point = [start_point[0]+math.tan(math.radians(angle)), start_point[1]-1+y_end-f_correction(start_point[0]+math.tan(math.radians(angle)))]
    first_point = [start_point[0], start_point[1]+y_end-f_correction(start_point[0])]
    
    return -math.degrees(math.atan(second_point[0]-first_point[0])/(second_point[1]-first_point[1]))

gStyle.SetOptStat(0)
gStyle.SetOptFit(0)

#f_mip = TFile("mip_fit.root")
#mip_cut = gDirectory.Get("f_pol")

f = TFile(filename)

chain = gDirectory.Get("events_"+algo)
entries = chain.GetEntriesFast()
h_theta_res_xy_corr = TH1F("h_theta_res_xy_corr", ";#Delta#theta_{xy} [#circ];N. Entries / 0.2#circ", 100, -10, 10)
h_theta_res_yz_corr = TH1F("h_theta_res_yz_corr", ";#Delta#theta_{yz} [#circ];N. Entries / 0.2#circ", 100, -10, 10)

h_theta_res_xy = TH1F("h_theta_res_xy", ";#Delta#theta_{xy} [#circ];N. Entries / 0.2#circ", 100, -10, 10)
h_theta_res_yz = TH1F("h_theta_res_yz", ";#Delta#theta_{yz} [#circ];N. Entries / 0.2#circ", 100, -10, 10)

h_theta_res_xy_mc = TH1F("h_theta_res_xy_mc", ";#Delta#theta_{xy} [#circ];N. Entries / 0.2#circ", 100, -10, 10)
h_theta_res_yz_mc = TH1F("h_theta_res_yz_mc", ";#Delta#theta_{yz} [#circ];N. Entries / 0.2#circ", 100, -10, 10)

h_res_xy = TH2F("h_res_xy", ";x [cm]; y [cm]", 80, -8, 8, 200, x_start, x_end)


h_xy = TH2F("h_xy", ";x [cm]; y [cm]", 800, x_start-20, x_end+20, 800, y_start-40, y_end+40)
h_xy_mc = TH2F("h_xy_mc", ";x [cm]; y [cm]", 800, x_start-20, x_end+20, 800, y_start-40, y_end+40)
h_xy_corr = TH2F("h_xy_corr",";x [cm]; y [cm]", 800, x_start-20, x_end+20, 800, y_start-40, y_end+40)

h_yz = TH2F("h_yz", ";z [cm]; y [cm]", 400, z_start-40, z_end+40, 400, y_start-40, y_end+40)
h_yz_mc = TH2F("h_yz_mc", ";z [cm]; y [cm]", 400, z_start-40, z_end+40, 400, y_start-40, y_end+40)
h_yz_corr = TH2F("h_yz_corr", ";z [cm]; y [cm]", 400, z_start-40, z_end+40, 400, y_start-40, y_end+40)

h_xy_res = TH2F("h_xy_res",";#Delta x [cm];#Delta y [cm]", 500, -50, 50, 500, -50, 50)
h_xy_res_corr = TH2F("h_xy_res_corr",";#Delta x [cm];#Delta y [cm]", 500, -50, 50, 500, -50, 50)

h_l_res = TH1F("h_l_res",";#Delta Length [cm];N. Entries / 2 cm", 50, -50, 50)
h_l_res_mc = TH1F("h_l_res_mc",";#Delta Length [cm];N. Entries / 2 cm",50, -50, 50)
h_l_res_corr = TH1F("h_l_res_corr",";#Delta Length [cm];N. Entries / 2 cm",50, -50, 50)

h_l_theta = TH2F("h_l_theta", ";#Delta Length [cm]; #theta_xy [#circ]", 50, -20, 100, 100, -100, -20)


h_xys = []

h_xys_end = []
binning = 25

print(int(x_end/binning))
for i in range(int(x_end/binning)+1):
    h = TH1F("h%d"%i,"",200,0,y_end+100)
    h_end = TH1F("h_end%d"%i,"",200,y_start-100,0)
    h_xys.append(h)
    h_xys_end.append(h_end)
print(entries)
for entry in range(entries-162):
    ientry = chain.LoadTree(entry)
    nb = chain.GetEntry(entry)
    
    x = chain.Tagged_Start[0]
    y = chain.Tagged_Start[1]
    z = chain.Tagged_Start[2]
    
    x2 = chain.Tagged_End[0]
    y2 = chain.Tagged_End[1]
    z2 = chain.Tagged_End[2] 

    if x:
    #   x+=25.35
    #   x2+=25.35
        x-=2
        x2-=2
    
    cosmic_pe = chain.flash_pe
    l = chain.MuCS_TPC_len
    x_mucs = chain.MuCS_Start_TPC[0]
    l_tpc = chain.Tagged_len
    mip = l and cosmic_pe

    if chain.Tagged_Start[0] and mip and chain.MuCS_NHitsX < 8 and chain.MuCS_NHitsZ < 8:
        try:
            h_xy.Fill(x,y)
            h_xy.Fill(x2,y2)
            h_xys[int(x/binning)].Fill(y)
            h_xys_end[int(x2/binning)].Fill(y2)
        except IndexError:
            print("Out of range", x2, int(x2/binning))

max_start = array("d",[])
max_end = array("d",[])

errors_y_start = array("d",[])
errors_y_end = array("d",[])


for i,h in enumerate(h_xys):
    if h.GetEntries()<10:
        max_start.append(h.GetMean())
    else:
        max_start.append(h.GetXaxis().GetBinCenter(h.GetMaximumBin()))
    try:
        errors_y_start.append(h.GetRMS()/math.sqrt(h.GetEntries()))
    except:
        continue
        
max_start[9] = 104
errors_x = array("d", [binning/2]*len(h_xys))
index = array("d",[(i+0.5)*binning for i,h in enumerate(h_xys)])
#errors_y_start[9] = 0.1
g = TGraphErrors(len(index)-1,index,max_start,errors_x,errors_y_start)
g.GetXaxis().SetTitle("x [cm]")
g.GetYaxis().SetTitle("y [cm]")

g.SetMarkerStyle(20)
f_pol = TF1("f_pol","[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4",0,x_end)
f_pol.SetParameters(1.15930e+02, 2.54863e-02, -1.6800e-04, -3.23142e-06, 1.03878e-08)
g.Fit("f_pol","R")
f_pol.SaveAs(output+"/pol_top.root")
for h in h_xys_end:
    max_end.append(h.GetXaxis().GetBinCenter(h.GetMaximumBin()))#h.GetXaxis().GetBinCenter(h.GetMaximumBin()))
    errors_y_end.append(h.GetRMS()/math.sqrt(h.GetEntries()))

for i in range(len(max_start)):
    print(i, "TOP", max_start[i], errors_y_start[i])
    print(i, "BOT", -max_end[i], errors_y_end[i])

c_spacecharge = TCanvas("c_spacecharge","",600,600)
c_spacecharge.SetLeftMargin(0.15)
c_spacecharge.SetRightMargin(0.15)
c_spacecharge.SetTopMargin(0.15)
c_spacecharge.SetBottomMargin(0.15)

#g.Draw("AP")
g_end = TGraphErrors(len(index)-1,index,max_end,errors_x,errors_y_end)
#g_end.Draw("P")
g_end.SetMarkerStyle(20)
f_pol2 = TF1("f_pol2","[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4",0,x_end)
f_pol2.SetParameters(-1.14317e+02, -7.02387e-02, 8.15558e-04, 6.27013e-07,-7.39476e-09)
g_end.Fit("f_pol2","R")


g.GetXaxis().SetLimits(x_start-40,x_end+40)
g.GetYaxis().SetRangeUser(y_start-30,y_end+30)

f_y = TF1("f_y","0.0060304-0.0779337*x+1.39288e-06*(x**2)-2.4147e-6*(x**3)")
a,b,c,d,e = f_pol.GetParameter(0),f_pol.GetParameter(1),f_pol.GetParameter(2),f_pol.GetParameter(3),f_pol.GetParameter(4)
a2,b2,c2,d2,e2 = f_pol2.GetParameter(0),f_pol2.GetParameter(1),f_pol2.GetParameter(2),f_pol2.GetParameter(3),f_pol2.GetParameter(4)

y_value = str(f_y.GetExpFormula()).replace("x","(116.5/([0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4)*y)")
y_value_neg = str(f_y.GetExpFormula()).replace("x","-116.5/([0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4)*y")

y_max =  str(f_y.GetExpFormula()).replace("x","116.5")
y_max_neg =  str(f_y.GetExpFormula()).replace("x","-116.5")
#y_value = "1"
#y_max = "1"
#y_value_neg = "1"
#y_max_neg = "1"
print(y_value,y_max)
f3 = TF2("f_corr_pos","y+(116.5-([0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4))*"+"("+y_value+")/("+y_max+")",0,x_end,y_start+10,y_end-10)
f3.SetParameter(0,a)
f3.SetParameter(1,b)
f3.SetParameter(2,c)
f3.SetParameter(3,d)
f3.SetParameter(4,e)
f4 = TF2("f_corr_neg","y-(([0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4)+116.5)*"+"("+y_value_neg+")/("+y_max_neg+")",0,x_end,y_start+10,-1)
f4.SetParameter(0,a2)
f4.SetParameter(1,b2)
f4.SetParameter(2,c2)
f4.SetParameter(3,d2)
f4.SetParameter(4,e2)

f3.SetLineWidth(1)

f3.SaveAs(output+"/top_correction.root")
f4.SaveAs(output+"/bottom_correction.root")

c_spacecharge.Update()
average_xy = 0
n_xy = 0
average_yz = 0
n_yz = 0
for entry in range(entries-162):
    #ientry = chain.LoadTree(entry)
    
    nb = chain.GetEntry(entry)
    
    l = chain.MuCS_TPC_len
    x_mucs = chain.MuCS_Start_TPC[0]
    cosmic_pe = chain.flash_pe

    mip = l and cosmic_pe
    
    x = chain.Tagged_Start[0]
    y = chain.Tagged_Start[1]
    z = chain.Tagged_Start[2]

    x2 = chain.Tagged_End[0]
    y2 = chain.Tagged_End[1]
    z2 = chain.Tagged_End[2]

    x_mc = chain.MCTagged_Start[0]
    y_mc = chain.MCTagged_Start[1]
    z_mc = chain.MCTagged_Start[2]
    
    x2_mc = chain.MCTagged_End[0]
    y2_mc = chain.MCTagged_End[1]
    z2_mc = chain.MCTagged_End[2]
    
    y_corr = 0
    y_corr_2 = 0

    #box_cut = (f3(x,y) > 110 or f4(x,y) < -110 or x > 254) and (f3(x2,y2) > 110 or f4(x2,y2) < -110 or x2 > 254)
    box_cut = 1
    if x and mip and chain.MuCS_NHitsX < 8 and chain.MuCS_NHitsZ < 8 and box_cut:
        #x+=25.35
        #x2+=25.35
        #x_mc+=25.35
        #x2_mc+=25.35
        x-=2
        x2-=2
        max_value = f_pol(x)
        max_value_end = f_pol2(x2)
        if y > 0:
            y_corr = f3(x,y)
        else:
            y_corr = f4(x,y)
        if y2 > 0:
            y_corr_2 = f3(x2,y2)
        else:
            y_corr_2 = f4(x2,y2)
            

        h_xy_mc.Fill(x_mc,y_mc)
        h_xy_mc.Fill(x2_mc,y2_mc)

        #h_xy_corr.Fill(x,y)
        #h_xy_corr.Fill(x2,y2)
        
        h_xy_corr.Fill(x,y_corr)
        h_xy_corr.Fill(x2,y_corr_2)
        
        h_yz.Fill(z,y)
        h_yz.Fill(z2,y2)
        
        h_yz_mc.Fill(z_mc,y_mc)
        h_yz_mc.Fill(z2_mc,y2_mc)
        
        h_yz_corr.Fill(z,y_corr)
        h_yz_corr.Fill(z2,y_corr_2)

        xy_mucs = 0
        yz_mucs = 0
        
        h_xy_res.Fill(x2_mc-x2,y2_mc-y2)
        h_xy_res_corr.Fill(x2_mc-x2,y2_mc-y_corr_2)
        
        if chain.MuCS_Start_TPC[0]-chain.MuCS_Start[0] and chain.MuCS_Start_TPC[2]-chain.MuCS_Start[2]:
            xy_mucs = math.degrees(math.atan((chain.MuCS_Start_TPC[1]-chain.MuCS_Start[1])/(chain.MuCS_Start_TPC[0]-chain.MuCS_Start[0])))
            yz_mucs = math.degrees(math.atan((chain.MuCS_Start_TPC[1]-chain.MuCS_Start[1])/(chain.MuCS_Start_TPC[2]-chain.MuCS_Start[2])))
        
        
        if x2_mc-x_mc and z2_mc-z_mc:
            xy_mc = math.degrees(math.atan((y2_mc-y_mc)/(x2_mc-x_mc)))
            yz_mc = math.degrees(math.atan((y2_mc-y_mc)/(z2_mc-z_mc)))

            if xy_mucs: 
                average_xy += math.degrees(chain.MuCS_theta_xy_rms)
                n_xy += 1
                h_theta_res_xy_mc.Fill(xy_mucs-xy_mc+random.gauss(0,math.degrees(chain.MuCS_theta_xy_rms)))
                h_l_res_mc.Fill(chain.MuCS_TPC_len-chain.MCTagged_len)
            if yz_mucs: 
                average_yz += math.degrees(chain.MuCS_theta_yz_rms)
                n_yz += 1
                h_theta_res_yz_mc.Fill(yz_mucs-yz_mc+random.gauss(0,math.degrees(chain.MuCS_theta_yz_rms)))

        if x-x2 and z-z2:
            theta_xy = math.atan((y-y2)/(x-x2))
            theta_yz = math.atan((y-y2)/(z-z2))+0.034
            theta_xy_corr = math.atan((y_corr-y_corr_2)/(x-x2))
            theta_yz_corr = math.atan((y_corr-y_corr_2)/(z-z2))+0.034
            
            if yz_mucs:
                h_theta_res_yz.Fill(yz_mucs-math.degrees(theta_yz))
                h_theta_res_yz_corr.Fill(yz_mucs-math.degrees(theta_yz_corr))
                
            if xy_mucs: 
                h_l_theta.Fill(chain.MuCS_TPC_len-chain.Tagged_len,math.degrees(theta_xy))
                if x2 > 250:
                    h_l_res.Fill(chain.MuCS_TPC_len-8./math.cos(-theta_xy)-chain.Tagged_len)
                else:
                    h_l_res.Fill(chain.MuCS_TPC_len-chain.Tagged_len)
                    
                h_l_res_corr.Fill(chain.MuCS_TPC_len-chain.Tagged_len)
                
                #if chain.MuCS_TPC_len-chain.Tagged_len > 60  and chain.MuCS_TPC_len-chain.Tagged_len < 60.2:
                #    print chain.evt_number, xy_mucs, math.degrees(theta_xy)+1.7, x, y_corr      
                #    event = gDirectory.Get("tdviewpandoraCosmic_"+str(chain.evt_number))
                #    event.Draw()
                h_res_xy.Fill(xy_mucs-math.degrees(theta_xy_corr)+1.25, x_mucs)
                h_theta_res_xy.Fill(xy_mucs-math.degrees(theta_xy)+1.25)
                h_theta_res_xy_corr.Fill(xy_mucs-math.degrees(theta_xy_corr)+1.25)
                

print("xy",average_xy/n_xy)
print("yz",average_yz/n_yz)

l = TLegend(0.63,0.7,0.76,0.86)
c_res_xy = TCanvas("c_res_xy")
h_theta_res_xy_corr.Scale(h_theta_res_xy_mc.Integral()/h_theta_res_xy.Integral())
h_theta_res_xy.Scale(h_theta_res_xy_mc.Integral()/h_theta_res_xy.Integral())

h_theta_res_xy_mc.Draw("hist")
h_theta_res_xy_mc.SetLineColor(kRed+1)
h_theta_res_xy_mc.SaveAs(output+"/xy_mc_"+algo+"_res.root")

h_theta_res_xy.Draw("epsame")
h_theta_res_xy_corr.Draw("epsame")
h_theta_res_xy_corr.SaveAs(output+"/xy_tpc_"+algo+"_res.root")

h_theta_res_xy_corr.SetMarkerStyle(20)
h_theta_res_xy_corr.SetLineColor(kBlack)
h_theta_res_xy.SetMarkerStyle(4)
h_theta_res_xy.SetLineColor(kBlack)
print("no_corr",h_theta_res_xy.GetMean()-h_theta_res_xy_mc.GetMean())
print("corr",h_theta_res_xy_corr.GetMean()-h_theta_res_xy_mc.GetMean())

l.AddEntry(h_theta_res_xy_mc, "Monte Carlo", "l")
l.AddEntry(h_theta_res_xy, "Data", "ep")
l.AddEntry(h_theta_res_xy_corr, "Data - SCE corrected", "ep")
l.SetTextFont(43)
l.SetShadowColor(0)
l.SetBorderSize(0)
l.SetTextSize(14)
l.Draw()
preliminary2 = TPaveText(0.13,0.78,0.46,0.86,"NDC")
preliminary2.SetShadowColor(0)
preliminary2.SetBorderSize(0)
preliminary2.SetFillStyle(0)
preliminary2.AddText("MicroBooNE Preliminary")
preliminary2.Draw()
c_res_xy.Update()
c_res_xy.SaveAs(output+"/xy_"+algo+"_res.pdf")

c_res_yz = TCanvas("c_res_yz")
f_eff_yz = TFile("yz_mc_pandoraCosmic_eff.root")
eff_yz = gDirectory.Get("e_yz_mc")
h_theta_res_yz_corr.Scale(h_theta_res_yz_mc.Integral()/h_theta_res_yz.Integral())
h_theta_res_yz.Scale(h_theta_res_yz_mc.Integral()/h_theta_res_yz.Integral())
h_theta_res_yz_mc.SetLineColor(kRed+1)
h_theta_res_yz_mc.Draw("hist")
h_theta_res_yz_mc.SaveAs(output+"/yz_mc_"+algo+"_res.root")

h_theta_res_yz.Draw("epsame")
h_theta_res_yz_corr.Draw("epsame")

h_theta_res_yz_corr.SetMarkerStyle(20)
h_theta_res_yz_corr.SetLineColor(kBlack)
h_theta_res_yz.SetMarkerStyle(4)
h_theta_res_yz.SetLineColor(kBlack)

h_theta_res_yz_corr.SaveAs(output+"/yz_tpc_"+algo+"_res.root")
l.Draw()
print("no_corr",h_theta_res_yz.GetMean()-h_theta_res_yz_mc.GetMean())
print("corr",h_theta_res_yz_corr.GetMean()-h_theta_res_yz_mc.GetMean())

preliminary2.Draw()
c_res_yz.Update()
c_res_yz.SaveAs(output+"/yz_"+algo+"_res.pdf")


c_res_l = TCanvas("c_res_l")
h_l_res_mc.SetLineColor(kRed+1)
h_l_res_mc.Draw("hist")

h_l_res.Draw("epsame")
h_l_res_corr.Draw("epsame")

h_l_res_corr.SetMarkerStyle(20)
h_l_res_corr.SetLineColor(kBlack)
h_l_res.SetMarkerStyle(4)
h_l_res.SetLineColor(kBlack)
h_l_res_corr.SaveAs(output+"/l_tpc_"+algo+"_res.root")
h_l_res_mc.SaveAs(output+"/l_mc_"+algo+"_res.root")

l.Draw()
c_res_l.Update()
c_res_l.SaveAs(output+"/l_"+algo+"_res.pdf")

c_xy = TCanvas("c_xy")

h_xy.Draw()

h_xy_mc.Draw("same")
h_xy_corr.Draw("same")
h_xy_corr.SetMarkerColor(kBlue+1)
h_xy_mc.SetMarkerColor(kRed+1)
h_xy.GetYaxis().SetTitleOffset(1.37);
c_xy.Update()

c_spacecharge.cd()

tpc_xy = TBox(x_start, y_start, x_end, y_end)
tpc_xy.SetLineWidth(3)
tpc_xy.SetLineStyle(2)
tpc_xy.SetFillStyle(0)
tpc_xy.SetFillStyle(0)
tpc_xy.SetLineWidth(3)
#h_xy.Draw("col")
h_xy.GetYaxis().SetTitleOffset(1.3)

h_xy_corr.Draw("col")
h_xy_corr.GetYaxis().SetTitleOffset(1.3)
preliminary = TPaveText(0.48,0.74,0.81,0.84,"NDC")
preliminary.SetShadowColor(0)
preliminary.SetBorderSize(0)
preliminary.SetFillStyle(0)
preliminary.AddText("MicroBooNE Preliminary")
preliminary.Draw()
tpc_xy.Draw()

g_end.Draw("P")
g.Draw("P")
f3.SetLineStyle(1)
f3.Draw("cont2same")

c_spacecharge.Update()
c_spacecharge.SaveAs(output+"/spacecharge_"+algo+".pdf")
c_yz = TCanvas("c_yz")
h_yz.Draw()
#h_yz_mc.Draw("same")
#h_yz_corr.Draw("same")
h_yz_mc.SetMarkerColor(kRed+1)
h_yz_corr.SetMarkerColor(kBlue+1)
c_yz.Update()
#h_l_theta.Draw()
#h_xys[2].Draw()

#f3.Draw("surf1")


c_2d = TCanvas("c_2d")
h_res_xy.Draw("colz")
c_2d.Update()

c_f= TCanvas("c_f")
f_pol.Draw()
f_pol.GetXaxis().SetRangeUser(0,x_end)
f_pol.GetYaxis().SetRangeUser(y_start-15,y_end+15)
t = TPaveText(0.1,0.3,0.4,0.5, "ndc")
t.SetBorderSize(0)
t.SetShadowColor(0)
t.SetFillColor(0)
t.AddText("f_{top}(x) = %.1e %.1e x+" % (f_pol.GetParameter(0),f_pol.GetParameter(1)))
t.AddText("+%.1e x^{2} %.1e x^{3}+%.1e x^{4}" % (f_pol.GetParameter(2), f_pol.GetParameter(3),f_pol.GetParameter(4)))
t.Draw()
t2 = TPaveText(0.1,0.3,0.4,0.5, "ndc")
t2.SetBorderSize(0)
t2.SetShadowColor(0)
t2.SetFillColor(0)
t2.AddText("f_{bottom}(x) = %.1e %.1e x+" % (f_pol2.GetParameter(0),f_pol2.GetParameter(1)))
t2.AddText("+%.1e x^{2} %.1e x^{3}+%.1e x^{4}" % (f_pol2.GetParameter(2), f_pol2.GetParameter(3),f_pol2.GetParameter(4)))

t2.Draw()
f_pol2.Draw("same")
c_f.Update()

c_g = TCanvas("c_g")
f_g = TF1("f_g","abs(pol3)",y_start,y_end)

f_g.SetParameter(0, 6.034060e-03/12.9222414623)
f_g.SetParameter(1, -7.793370e-02/12.9222414623)
f_g.SetParameter(2, 1.392878e-06/12.9222414623)
f_g.SetParameter(3, -2.414689e-06/12.9222414623)
f_g.Draw()
t3 = TPaveText(0.1,0.3,0.4,0.5, "ndc")
t3.SetBorderSize(0)
t3.SetShadowColor(0)
t3.SetFillColor(0)
t3.AddText("g(y) = |%.1e %.1e y+" % (f_g.GetParameter(0),f_g.GetParameter(1)))
t3.AddText("+%.1e y^{2} %.1e y^{3}|" % (f_g.GetParameter(2), f_g.GetParameter(3)))
t3.Draw()
print(f_g(y_start))
f_g.GetYaxis().SetTitle("g(y)")
f_g.GetXaxis().SetTitle("y [cm]")

c_g.Update()

input()