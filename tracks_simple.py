#!/usr/bin/env python3.4

import os
from ROOT import TChain,  gStyle, gDirectory, TFile, TPad, TCanvas, TPolyLine3D, kGreen, kRed, kBlue, TLegend, kAzure, kOrange, kMagenta
from array import array
import math


c1 = TCanvas("c1")
p1 = TPad("p1","p1",0.05,0.05,0.95,0.95)
p1.Draw()
p1.cd()

red_lines = []
blue_lines = []
green_lines = []
red_lines_2D = []
blue_lines_2D = []
green_lines_2D = []

chain = TChain("events_pandoraCosmic")
chain.Add("../mucs_merged/CORSIKA.root")
chain.Add("../mucs_merged/MuCSRun7347_Group181_MergedTree.Root") #downstream
chain.Add("../mucs_merged/MuCSRun7348_Group181_MergedTree.Root") #downstream
chain.Add("../mucs_merged/MuCSRun7702_Group182_MergedTree.Root") #upstream

chain.Add("../mucs_merged/MuCSRun7703_Group183_MergedTree.Root") #upstream

entries = chain.GetEntries()

print(entries)
for entry in range(entries):
    ientry = chain.LoadTree(entry)
    nb = chain.GetEntry(entry)
    cosmic_pe = chain.flash_pe
    l_mucs = chain.MuCS_TPC_len
    if entry % 1000 == 0: print(entry)
    if chain.MuCS_NHitsX < 8 and chain.MuCS_NHitsZ < 8 and l_mucs > 0 and cosmic_pe > 0 and entry % 10 == 0:
        x = chain.MuCS_Start_TPC[0]-chain.MuCS_Start[0]
        y = chain.MuCS_Start_TPC[1]-chain.MuCS_Start[1]
        z = chain.MuCS_Start_TPC[2]-chain.MuCS_Start[2]
        r = math.sqrt(x**2+y**2+z**2)
    
        theta_mucs = math.acos(z/r)
        theta_mucs = math.degrees(theta_mucs)
    
        phi_mucs = math.atan2(y,x)
        phi_mucs = math.degrees(phi_mucs)
        
        if theta_mucs > 75 and theta_mucs < 90 and phi_mucs < -75 and phi_mucs > -90 and l_mucs > 140 and l_mucs < 260:
        
            start = [chain.MuCS_Start_TPC[0], chain.MuCS_Start_TPC[1],chain.MuCS_Start_TPC[2]]
            end_tpc = [chain.MuCS_Start[0], chain.MuCS_Start[1],chain.MuCS_Start[2]]
            direction = [a_i - b_i for a_i, b_i in zip(end_tpc, start)]
            end = [chain.MuCS_Start_TPC[0]-direction[0], chain.MuCS_Start_TPC[1]-direction[1],chain.MuCS_Start_TPC[2]-direction[2]]
            
            polyline = TPolyLine3D(2)
            polyline.SetPoint(0, start[0], start[1], start[2])
            polyline.SetPoint(1, end[0], end[1], end[2])

            #red_lines.append(polyline)
            
        if theta_mucs > 90 and theta_mucs < 105 and phi_mucs > -90 and phi_mucs < -45 and l_mucs > 200 and l_mucs < 260:
        
            start = [chain.MuCS_Start_TPC[0], chain.MuCS_Start_TPC[1],chain.MuCS_Start_TPC[2]]
            end_tpc = [chain.MuCS_Start[0], chain.MuCS_Start[1],chain.MuCS_Start[2]]
            direction = [a_i - b_i for a_i, b_i in zip(end_tpc, start)]
            end = [chain.MuCS_Start_TPC[0]-direction[0], chain.MuCS_Start_TPC[1]-direction[1],chain.MuCS_Start_TPC[2]-direction[2]]
            
            #start = [chain.Tagged_Start[0], chain.Tagged_Start[1],chain.Tagged_Start[2]]
            #end = [chain.Tagged_End[0], chain.Tagged_End[1],chain.Tagged_End[2]]
            
            polyline = TPolyLine3D(2)
            polyline.SetPoint(0, start[0],start[1],start[2])
            polyline.SetPoint(1, end[0], end[1], end[2])

            #blue_lines.append(polyline)
                     
        if theta_mucs > 60 and theta_mucs < 75 and phi_mucs < -45 and phi_mucs > -60 and ((l_mucs > 200 and l_mucs < 260) or (l_mucs > 80 and l_mucs < 140)):
        
            start = [chain.MuCS_Start_TPC[0], chain.MuCS_Start_TPC[1],chain.MuCS_Start_TPC[2]]
            end_tpc = [chain.MuCS_Start[0], chain.MuCS_Start[1],chain.MuCS_Start[2]]
            direction = [a_i - b_i for a_i, b_i in zip(end_tpc, start)]
            end = [chain.MuCS_Start_TPC[0]-direction[0], chain.MuCS_Start_TPC[1]-direction[1],chain.MuCS_Start_TPC[2]-direction[2]]
            
            polyline = TPolyLine3D(2)
            polyline.SetPoint(0, start[0],start[1],start[2])
            polyline.SetPoint(1, end[0], end[1], end[2])

            green_lines.append(polyline)
    
            
    
for line in red_lines:
    line.SetLineColor(kRed)
    line.Draw()
for line in blue_lines:
    line.SetLineColor(kBlue)
    line.Draw()
for line in green_lines:
    line.SetLineColor(kGreen)
    line.Draw()    
    
TPC = TPolyLine3D(10)
TPC.SetPoint(0,0,116.5,0)
TPC.SetPoint(1,256.35,116.5,0)
TPC.SetPoint(2,256.35,116.5,1036.8)
TPC.SetPoint(3,0,116.5,1036.8)
TPC.SetPoint(4,0,116.5,0)
TPC.SetPoint(5,0,-116.5,0)
TPC.SetPoint(6,256.35,-116.5,0)
TPC.SetPoint(7,256.35,-116.5,1036.8)
TPC.SetPoint(8,0,-116.5,1036.8)
TPC.SetPoint(9,0,-116.5,0)

TPC2 = TPolyLine3D(2)
TPC2.SetPoint(0,256.35,116.5,1036.8)
TPC2.SetPoint(1,256.35,-116.5,1036.8)

TPC3 = TPolyLine3D(2)
TPC3.SetPoint(0,0,116.5,1036.8)
TPC3.SetPoint(1,0,-116.5,1036.8)

TPC4 = TPolyLine3D(2)
TPC4.SetPoint(0,256.35,116.5,0)
TPC4.SetPoint(1,256.35,-116.5,0)

TPC.SetLineWidth(3)
TPC.SetLineColor(kRed)
TPC2.SetLineWidth(3)
TPC2.SetLineColor(kRed)
TPC3.SetLineWidth(3)
TPC3.SetLineColor(kRed)
TPC4.SetLineWidth(3)
TPC4.SetLineColor(kRed)

Top = TPolyLine3D(9)
Top.SetPoint(0,-304.6,571.7,604.8)
Top.SetPoint(1,-304.6,571.7,950.9)
Top.SetPoint(2,-131.6,571.7,950.9)
Top.SetPoint(3,-131.6,571.7,1124)
Top.SetPoint(4,387.6,571.7,1124)
Top.SetPoint(5,387.6,571.7,-87.3)
Top.SetPoint(6,-131.6,571.7,-87.3)
Top.SetPoint(7,-131.6,571.7,604.8)
Top.SetPoint(8,-304.6,571.7,604.8)
Top.SetLineWidth(3)

Under = TPolyLine3D(5)
Under.SetPoint(0,387.6,	-253.4,	210.4)
Under.SetPoint(1,-131.6,	-253.4,	210.4)
Under.SetPoint(2,-131.6,	-253.4,	556.5)
Under.SetPoint(3,387.6,	-253.4,	556.5)
Under.SetPoint(4,387.6,	-253.4,	210.4)
Under.SetLineWidth(3)



Pipe = TPolyLine3D(7)
Pipe.SetPoint(0,378.3,	252.9,	1192.2)
Pipe.SetPoint(1,378.3,	252.9,	-19.1)
Pipe.SetPoint(2,378.3,	99,	-19.1)
Pipe.SetPoint(3,380.1,	99,	-192.1)
Pipe.SetPoint(4,385.7,	-247.1,	-192.2)
Pipe.SetPoint(5,387.6,	-247.1,	1192.2)
Pipe.SetPoint(6,378.3,	252.9,	1192.2)

Pipe.SetLineWidth(3)

FT = TPolyLine3D(5)
FT.SetPoint(0,-129.7,	124.6,	-84.9)
FT.SetPoint(1,-129.7,	124.6,	1126.4)
FT.SetPoint(2,-131.6,	-221.5,	1126.4)
FT.SetPoint(3,-131.6,	-221.5,	-84.9)
FT.SetPoint(4,-129.7,	124.6,	-84.9)
FT.SetLineWidth(3)


TPC.Draw()
TPC2.Draw()
TPC3.Draw()
TPC4.Draw()
#Top.Draw()
#Under.Draw()
#Pipe.Draw()
#FT.Draw()
c1.Update()
input()