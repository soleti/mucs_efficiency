#!/usr/bin/env python3.4

from ROOT import TGraphErrors, kBlue, kRed, TCanvas
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

#g_purity = TGraphErrors(len(purity_list), array(tolerance_list), array(purity_list))

g_eff.Draw("APL")
g_eff.SetMarkerStyle(20)
g_eff.SetLineColor(kRed+1)
g_eff.SetLineWidth(2)
g_purity.Draw("PL")
g_purity.SetMarkerStyle(20)
g_purity.SetLineColor(kBlue+1)
g_purity.SetLineWidth(2)

c.Update()
input()