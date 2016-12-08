#!/usr/bin/env python3.4

from ROOT import TGraphErrors, TFile, TCanvas, kGray, TH2F, TH1F, gDirectory, kRed, kGreen, gStyle, TPad, TLine, TLegend, gPad, TPaveText
from numpy import array
gStyle.SetOptStat(0)
gStyle.SetNumberContours(999)
gStyle.SetPaintTextFormat(".2f")

def histo2graph(histo):
    bins = array([histo.GetBinCenter(i) for i in range(1,histo.GetNbinsX()+2)])
    values = array([histo.GetBinContent(i) for i in range(1,histo.GetNbinsX()+2)])
    bins_err = array([histo.GetBinWidth(1)/2.]*histo.GetNbinsX())
    values_err = array([histo.GetBinError(i) for i in range(1,histo.GetNbinsX()+2)])
    graph = TGraphErrors(histo.GetNbinsX(),bins,values,bins_err,values_err)
    return graph


def draw_ratio(c_title, h_top = [], h_bottom = [], draw_top = [], draw_bottom = []):
    c = TCanvas(c_title)
    c.cd()
    pad_bottom = TPad("pad_bottom","",0, 0.02, 1, 0.3)
    pad_bottom.SetTopMargin(0)
    pad_bottom.SetBottomMargin(0.23)
    pad_bottom.Draw()
    pad_bottom.cd()
    h_bottom[0].SetStats(0)
    pad_bottom.SetGridy()

    h_bottom[0].Draw(draw_bottom[0])
    h_bottom[0].GetYaxis().SetRangeUser(0.91,1.09)
    h_bottom[0].GetYaxis().SetNdivisions(5,0,5)
    h_bottom[0].GetYaxis().SetLabelSize(0.12)
    h_bottom[0].GetYaxis().SetTitleSize(0.13)
    h_bottom[0].GetYaxis().SetTitleOffset(0.3)
    h_bottom[0].GetXaxis().SetLabelSize(0.12)
    h_bottom[0].GetXaxis().SetTitleSize(0.13)
    h_bottom[0].GetXaxis().SetTitleOffset(0.85)
    h_bottom[0].GetYaxis().SetTitle("Data/MC")
    h_bottom[0].SetMarkerStyle(20)
    h_bottom[0].SetLineColor(1)
    for i, h in enumerate(h_bottom[1:]):
        h.Draw("same"+draw_bottom[i+1])

    c.cd()
    pad_top = TPad("pad_top","",0, 0.3,1,1)
    pad_top.SetBottomMargin(0)
    pad_top.Draw()
    pad_top.cd()
    h_top[0].SetStats(0)
    h_top[0].DrawCopy(draw_top[0])
    h_top[0].GetYaxis().SetTitleOffset(0.7)
    h_top[0].SetFillStyle(3001)
    h_top[0].Draw("e2same")

    for i, h in enumerate(h_top[1:]):
        h.Draw("same"+draw_top[i+1])


    c.Update()
    return c

def draw_canvas_1d(name):
    f_mc = TFile("plots/mc/%s_mcc7.root" % name)
    h_mc = gDirectory.Get("h_%s_reco" % name)
    h_mc.SetName("h_%s_reco" % name)
    h_mc.SetLineColor(kRed+1)
    h_mc.SetFillColor(kRed+1)
    h_mc.SetFillStyle(0)
    h_mc.GetYaxis().SetTitle("Efficiency")
    h_mc.GetYaxis().SetRangeUser(0.61,1.04)
    f = TFile("plots/data/e_%s_pandoraCosmic.root" % name)
    h = gDirectory.Get("h_%s" % name)
    h.SetMarkerStyle(20)
    h.SetLineColor(1)
    f_sys = TFile("plots/data/e_%s_sys_pandoraCosmic.root" % name)
    h_sys = gDirectory.Get("h_%s_sys" % name)
    h_ratio = h_sys.Clone()
    h_ratio.Divide(h_mc)


    x_minbin = h_mc.FindFirstBinAbove()
    low = h_mc.GetXaxis().GetBinLowEdge(x_minbin)
    x_maxbin = h_mc.FindLastBinAbove()
    high = h_mc.GetXaxis().GetBinUpEdge(x_maxbin)
    h_mc.GetXaxis().SetRangeUser(low,high)
    h_ratio.GetXaxis().SetRangeUser(low,high)


    g_stat = histo2graph(h)
    g_sys = histo2graph(h_sys)

    g_sys.SetMarkerStyle(20)
    g_sys.SetLineWidth(2)
    g_stat.SetFillStyle(3001)
    g_stat.SetFillColor(kGray+2)
    g_stat.SetLineWidth(2)
    g_stat.SetLineColor(kGray+2)

    pt = TPaveText(0.10,0.905,0.40,0.98, "ndc")
    pt.AddText("MicroBooNE in progress")
    pt.SetFillColor(0)
    pt.SetBorderSize(0)
    pt.SetShadowColor(0)

    leg = TLegend(0.66,0.09,0.86,0.26)
    leg.AddEntry(g_sys,"Data - stat. #oplus sys.","ep")
    leg.AddEntry(g_stat,"Data - stat. only","f")
    leg.AddEntry(h_mc,"Monte Carlo","f")

    canvas = draw_ratio("c_%s" % name, [h_mc, g_stat, g_sys, pt, leg], [h_ratio], ["hist","2","p","",""],["ep"])
    canvas.SaveAs("plots/%s.pdf" % name)

def draw_canvas_2d(name):
    pt = TPaveText(0.10,0.905,0.40,0.98, "ndc")
    pt.AddText("MicroBooNE in progress")
    pt.SetFillColor(0)
    pt.SetBorderSize(0)
    pt.SetShadowColor(0)

    f_mc = TFile("plots/mc/%s_mcc7.root" % name)
    h_mc = gDirectory.Get("h_%s_reco" % name)
    f = TFile("plots/data/e_%s_pandoraCosmic.root" % name)
    h = gDirectory.Get("h_%s" % name)
    h.Divide(h_mc)
    h.GetZaxis().SetRangeUser(0.5,1.5)
    h.GetZaxis().SetTitle("Data/Monte Carlo")
    h.GetZaxis().RotateTitle()

    x_minbin = h.FindFirstBinAbove(0,1)
    low_x = h.GetXaxis().GetBinLowEdge(x_minbin)
    x_maxbin = h.FindLastBinAbove(0,1)
    high_x = h.GetXaxis().GetBinUpEdge(x_maxbin)

    y_minbin = h.FindFirstBinAbove(0,2)
    low_y = h.GetYaxis().GetBinLowEdge(y_minbin)
    y_maxbin = h.FindLastBinAbove(0,2)
    high_y = h.GetYaxis().GetBinUpEdge(y_maxbin)

    h.GetXaxis().SetRangeUser(low_x,high_x)
    h.GetYaxis().SetRangeUser(low_y,high_y)

    c = TCanvas("c_%s" % name)
    h.Draw("colz texte")

    pt.Draw()
    gPad.SetRightMargin(0.15)
    c.Update()
    c.SaveAs("plots/%s.pdf" % name)
    return c


draw_canvas_1d("theta")
draw_canvas_1d("phi")
draw_canvas_1d("l")

draw_canvas_2d("theta_phi")
draw_canvas_2d("theta_l")
draw_canvas_2d("phi_l")


gStyle.SetCanvasPreferGL(1)
f_theta_phi_l_mc = TFile("plots/mc/e_theta_phi_l_pandoraCosmic.root")
h_theta_phi_l_mc = gDirectory.Get("h_theta_phi_l_tpc")
f_theta_phi_l = TFile("plots/data/e_theta_phi_l_pandoraCosmic.root")
h_theta_phi_l = gDirectory.Get("h_theta_phi_l_tpc")
h_theta_phi_l.GetZaxis().SetTitleOffset(1.7)
h_theta_phi_l.GetXaxis().SetTitleOffset(1.7)
h_theta_phi_l.GetYaxis().SetTitleOffset(1.7)
h_theta_phi_l.Divide(h_theta_phi_l_mc)
h_theta_phi_l.GetZaxis().SetRangeUser(20,320)
h_theta_phi_l.GetXaxis().SetRangeUser(60,120)
h_theta_phi_l.GetYaxis().SetRangeUser(-90,-45)
c_theta_phi_l = TCanvas("c_theta_phi_l")
h_theta_phi_l.Draw("glbox")
c_theta_phi_l.Update()
c_theta_phi_l.SaveAs("plots/theta_phi_l.png")
input()
