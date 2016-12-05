#!/usr/bin/env python3.4

from ROOT import TGraphErrors, kBlue, kRed, kGreen, TCanvas, TLegend, TLine
from numpy import array
import math

def product_error(a,b,a_err,b_err):
    return math.sqrt((a*b)**2*((a_err/a)**2+(b_err/b)**2))

eff_list = []
eff_err_list = []
purity_list = []
purity_err_list = []
tolerance_list = []
product = []
product_err = []
f = open("output/purity.txt")
for line in f:
    tolerance, eff, eff_err, purity, purity_err = line.split()
    tolerance_list.append(float(tolerance))
    eff_list.append(float(eff)-2.66)
    eff_err_list.append(float(eff_err))
    purity_list.append(float(purity))
    purity_err_list.append(float(purity_err))
    product.append(float(purity)*(float(eff)-2.66)/100)
    product_err.append(100*product_error(float(purity)/100,float(eff)/100,float(purity_err)/100,float(eff_err)/100))

    print(float(purity)*float(eff))

no_err = array([0]*len(eff_list))
c = TCanvas("c")
g_eff = TGraphErrors(len(eff_list), array(tolerance_list), array(eff_list), no_err, array(eff_err_list))
g_purity = TGraphErrors(len(purity_list), array(tolerance_list), array(purity_list), no_err, array(purity_err_list))
g_product = TGraphErrors(len(product), array(tolerance_list), array(product), no_err, array(product_err))





g_product.Draw("APL")
g_product.GetXaxis().SetTitle("Tolerance [cm]")
g_product.GetYaxis().SetTitle("[%]")
g_product.GetYaxis().SetRangeUser(94, 100)

g_product.SetMarkerStyle(21)
g_product.SetLineColor(kGreen+1)
g_product.SetLineWidth(2)
g_eff.Draw("PL")
g_eff.SetMarkerStyle(20)
g_eff.SetLineColor(kRed+1)
g_eff.SetLineWidth(2)
g_purity.Draw("PL")
g_purity.SetMarkerStyle(21)
g_purity.SetLineColor(kBlue+1)
g_purity.SetLineWidth(2)
reco_eff_line = TLine(6,98.96-2.66,54,98.96-2.66)
reco_eff_line.SetLineColor(1)
reco_eff_line.SetLineStyle(2)
reco_eff_line.SetLineWidth(2)
reco_eff_line.Draw()

l = TLegend(0.46,0.16,0.89,0.33)
l.SetBorderSize(0)
l.SetShadowColor(0)
l.AddEntry(g_eff,"Efficiency", "lp")
l.AddEntry(g_purity,"Purity", "lp")
l.AddEntry(g_product,"Efficiency x Purity", "lp")

l.AddEntry(reco_eff_line,"Monte Carlo reco. efficiency", "l")

l.Draw()

c.Update()
input()
