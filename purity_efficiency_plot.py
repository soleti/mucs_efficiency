#!/usr/bin/env python3.4

from ROOT import TGraphErrors, TCanvas, TLegend, TLine
from ROOT import kBlue, kRed, kGreen, kGray
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
data_eff_list = []
reco_eff_list = []
f = open("output/purity.txt")
for line in f:
    tolerance, eff, eff_err, purity, purity_err, data_eff = line.split()
    tolerance_list.append(float(tolerance))
    eff_list.append(float(eff)-2.66)
    eff_err_list.append(float(eff_err))
    purity_list.append(float(purity))
    purity_err_list.append(float(purity_err))
    product_value = float(purity)*(float(eff)-2.66)/100
    acceptance_value = product_value/(98.96-2.66)
    acceptance.append(acceptance_value*100)
    data_eff_list.append(float(data_eff))
    print(float(data_eff)*float(purity)/(acceptance_value*100), tolerance)
    reco_eff_list.append(float(data_eff)*float(purity)/(acceptance_value*100))


reco_eff = array([96.3]*len(eff_list))
no_err = array([0]*len(eff_list))

c = TCanvas("c")
g_reco_eff = TGraphErrors(len(eff_list), array(tolerance_list), reco_eff, no_err, array(eff_err_list))
g_eff = TGraphErrors(len(eff_list), array(tolerance_list), array(eff_list), no_err, array(eff_err_list))
g_purity = TGraphErrors(len(purity_list), array(tolerance_list), array(purity_list), no_err, array(purity_err_list))
g_acceptance = TGraphErrors(len(acceptance),array(tolerance_list),array(acceptance),no_err, array(purity_err_list))
g_reco_eff_data = TGraphErrors(len(reco_eff_list),array(tolerance_list),array(reco_eff_list),no_err, array(eff_err_list))
g_data_eff = TGraphErrors(len(reco_eff_list),array(tolerance_list),array(data_eff_list),no_err, array(eff_err_list))

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
g_reco_eff.Draw("PL")
g_reco_eff.SetLineColor(kGray+2)
g_reco_eff.SetLineWidth(2)
g_reco_eff.SetMarkerStyle(29)
g_reco_eff.SetMarkerSize(1.5)
g_reco_eff_data.Draw("PL")
g_data_eff.Draw("PL")
g_data_eff.SetLineColor(kRed+1)
g_data_eff.SetLineWidth(2)
g_data_eff.SetMarkerStyle(4)
g_data_eff.SetLineStyle(2)
g_reco_eff_data.SetLineColor(kGray+2)
g_reco_eff_data.SetLineWidth(2)
g_reco_eff_data.SetMarkerStyle(30)
g_reco_eff_data.SetMarkerSize(1.5)
g_reco_eff_data.SetLineStyle(2)

l = TLegend(0.46,0.16,0.85,0.425)
l.SetBorderSize(0)
l.SetShadowColor(0)

l.AddEntry(g_purity,"Purity #it{P}", "lp")
l.AddEntry(g_acceptance,"Acceptance #it{A}", "lp")

l.AddEntry(g_eff,"MC tagging efficiency #it{#epsilon}_{tag}^{MC}", "lp")
l.AddEntry(g_data_eff,"Data tagging efficiency #it{#epsilon}_{tag}^{data}", "lp")

l.AddEntry(g_reco_eff,"MC reconstruction efficiency #it{#epsilon}_{MC}", "lp")
l.AddEntry(g_reco_eff_data,"Data reconstruction efficiency #it{#epsilon}_{data}", "lp")

l.Draw()

c.Update()
input()
