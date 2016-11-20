#!/usr/bin/env python3.4

from ROOT import TGraphErrors, kBlue, kRed, TCanvas, TLegend
from numpy import array

eff_list = []
eff_err_list = []
purity_list = []
purity_err_list = []
tolerance_list = []
f = open("purity.txt")
for line in f:
    tolerance, eff, eff_err, purity, purity_err = line.split()
    tolerance_list.append(float(tolerance))
    eff_list.append(float(eff))
    eff_err_list.append(float(eff_err))
    purity_list.append(float(purity))
    purity_err_list.append(float(purity_err))

no_err = array([0]*len(eff_list))
c = TCanvas("c")
g_eff = TGraphErrors(len(eff_list), array(tolerance_list), array(eff_list), no_err, array(eff_err_list))
g_purity = TGraphErrors(len(purity_list), array(tolerance_list), array(purity_list), no_err, array(purity_err_list))

g_eff.Draw("APL")
g_eff.GetXaxis().SetTitle("Tolerance [cm]")
g_eff.GetYaxis().SetTitle("[%]")

g_eff.SetMarkerStyle(20)
g_eff.SetLineColor(kRed+1)
g_eff.SetLineWidth(2)
g_purity.Draw("PL")
g_purity.SetMarkerStyle(21)
g_purity.SetLineColor(kBlue+1)
g_purity.SetLineWidth(2)

l = TLegend(0.66,0.16,0.89,0.27)
l.SetBorderSize(0)
l.SetShadowColor(0)
l.AddEntry(g_eff,"Efficiency", "lp")
l.AddEntry(g_purity,"Purity", "lp")
l.Draw()

c.Update()
input()