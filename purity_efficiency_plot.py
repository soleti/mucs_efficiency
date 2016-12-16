#!/usr/bin/env python3.4

from ROOT import TGraphErrors, TCanvas, TLegend, TLine
from ROOT import kBlue, kRed, kGreen
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
acceptance = []
f = open("output/purity.txt")
for line in f:
    tolerance, eff, eff_err, purity, purity_err = line.split()
    tolerance_list.append(float(tolerance))
    eff_list.append(float(eff)-2.66)
    eff_err_list.append(float(eff_err))
    purity_list.append(float(purity))
    purity_err_list.append(float(purity_err))
    product_value = float(purity)*(float(eff)-2.66)/100
    product.append(product_value)
    product_err.append(100*product_error(float(purity)/100,float(eff)/100,float(purity_err)/100,float(eff_err)/100))
    acceptance_value = product_value/(98.96-2.66)
    acceptance.append(acceptance_value*100)


no_err = array([0]*len(eff_list))
c = TCanvas("c")
g_eff = TGraphErrors(len(eff_list), array(tolerance_list), array(eff_list), no_err, array(eff_err_list))
g_purity = TGraphErrors(len(purity_list), array(tolerance_list), array(purity_list), no_err, array(purity_err_list))
g_product = TGraphErrors(len(product), array(tolerance_list), array(product), no_err, array(product_err))
g_acceptance = TGraphErrors(len(acceptance),array(tolerance_list),array(acceptance),no_err, array(product_err))


#g_product.Draw("APL")


g_product.SetMarkerStyle(22)
g_product.SetLineColor(kGreen+1)
g_product.SetLineWidth(2)
g_eff.Draw("APL")
g_eff.GetXaxis().SetTitle("#it{d}_{max} [cm]")
g_eff.GetYaxis().SetTitle("[%]")
g_eff.GetYaxis().SetRangeUser(90, 100)
g_eff.SetMarkerStyle(20)
g_eff.SetLineColor(kRed+1)
g_eff.SetLineWidth(2)
g_purity.Draw("PL")
g_purity.SetMarkerStyle(21)
g_purity.SetLineColor(kBlue+1)
g_purity.SetLineWidth(2)
g_acceptance.Draw("PL")
g_acceptance.SetMarkerStyle(23)
g_acceptance.SetLineColor(kGreen+2)
g_acceptance.SetLineWidth(2)
reco_eff_line = TLine(6,98.96-2.66,54,98.96-2.66)
reco_eff_line.SetLineColor(1)
reco_eff_line.SetLineStyle(2)
reco_eff_line.SetLineWidth(2)
reco_eff_line.Draw()

l = TLegend(0.46,0.16,0.89,0.33)
l.SetBorderSize(0)
l.SetShadowColor(0)
l.AddEntry(g_eff,"Tagging efficiency #it{#epsilon}_{tag}", "lp")
l.AddEntry(g_purity,"Purity #it{P}", "lp")
#l.AddEntry(g_product,"#it{#epsilon} x #it{p}", "lp")
l.AddEntry(g_acceptance,"Acceptance #it{A}", "lp")

l.AddEntry(reco_eff_line,"Reconstruction efficiency #it{#epsilon}", "l")

l.Draw()

c.Update()
input()
