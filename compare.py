#!/usr/bin/env python3.4

# Map
# 5 7 4 60.0 -90.0 200.0
# 5 7 5 60.0 -90.0 260.0
# 5 8 4 60.0 -75.0 200.0
# 5 8 5 60.0 -75.0 260.0
# 5 9 2 60.0 -60.0 80.0
# 5 9 3 60.0 -60.0 140.0
# 5 9 4 60.0 -60.0 200.0
# 5 9 5 60.0 -60.0 260.0
# 6 7 4 75.0 -90.0 200.0
# 6 8 4 75.0 -75.0 200.0
# 6 8 5 75.0 -75.0 260.0
# 6 9 1 75.0 -60.0 20.0
# 6 9 2 75.0 -60.0 80.0
# 6 9 3 75.0 -60.0 140.0
# 6 9 4 75.0 -60.0 200.0
# 6 9 5 75.0 -60.0 260.0
# 7 7 3 90.0 -90.0 140.0
# 7 7 4 90.0 -90.0 200.0
# 7 8 4 90.0 -75.0 200.0
# 7 8 5 90.0 -75.0 260.0
# 7 9 1 90.0 -60.0 20.0
# 7 9 2 90.0 -60.0 80.0
# 7 9 3 90.0 -60.0 140.0
# 7 9 4 90.0 -60.0 200.0
# 7 9 5 90.0 -60.0 260.0
# 8 7 4 105.0 -90.0 200.0
# 8 7 5 105.0 -90.0 260.0
# 8 8 4 105.0 -75.0 200.0
# 8 8 5 105.0 -75.0 260.0
# 8 9 2 105.0 -60.0 80.0
# 8 9 3 105.0 -60.0 140.0
# 8 9 4 105.0 -60.0 200.0
# 8 9 5 105.0 -60.0 260.0

import math
from ROOT import TH1F, TCanvas, gStyle

gStyle.SetOptStat(0)

h_sig = TH1F("h_sig",";Significance [#sigma];N. Entries / 1",20,-10,10)

def measure_sig(val1,val2):
    return (val1-val2)/math.sqrt(err1**2+err2**2)
    
    
sys = 0
sys_value = 0
print("0 - Central, 1 - upstream, 2 - downstream")
with open("output/central_shift.txt") as textfile1, open("output/upstream_shift.txt") as textfile2, open("output/downstream_shift.txt") as textfile3: 
    for x, y, z in zip(textfile1, textfile2, textfile3):
        a = x.split(" ")
        a = [float(i) for i in a]
        b = y.split(" ")
        b = [float(i) for i in b]
        c = z.split(" ")
        c = [float(i) for i in c]
        
        eff = [a[3],b[3],c[3]]
        err = [a[4],b[4],c[4]]
        theta = a[0]
        phi = a[1]
        if eff[0] != 0 and eff[1] != 0 and eff[0] != 1 and eff[1] != 1:
            if (theta == 7 and phi >= 7 and phi <= 9) or (theta == 5 and phi == 9):
                sys = sys_value
            sig = (eff[0]-eff[1])/math.sqrt(err[0]**2+err[1]**2+2*sys**2)
            print(err[0],err[1])
            h_sig.Fill(sig)
            if abs(sig) > 3:
                print("%i %i %.2f %.2f %.2f %.2f %.2f %.2f %.2f" % (0,1,eff[0],eff[1],eff[2],a[0],a[1],a[2],sig))
        
        if eff[0] != 0 and eff[2] != 0 and eff[0] != 1 and eff[2] != 1:
            sig = (eff[2]-eff[0])/math.sqrt(err[0]**2+err[2]**2+2*sys**2)
            h_sig.Fill(sig)
            if abs(sig) > 3:
                print("%i %i %.2f %.2f %.2f %.2f %.2f %.2f %.2f" % (0,2,eff[0],eff[1],eff[2],a[0],a[1],a[2],sig))
                
        if eff[1] != 0 and eff[2] != 0 and eff[1] != 1 and eff[2] != 1:
            if theta == 7 and phi >= 7 and phi <= 9:
                sys = sys_value
            sig = (eff[1]-eff[2])/math.sqrt(err[1]**2+err[2]**2+2*sys**2)
            if abs(sig) > 3:
                print("%i %i %.2f %.2f %.2f %.2f %.2f %.2f %.2f" % (1,2,eff[0],eff[1],eff[2],a[0],a[1],a[2],sig))
            else:
                h_sig.Fill(sig)

                    
                    

c_sig = TCanvas("c_sig")
h_sig.Draw("ep")
h_sig.SetMarkerStyle(20)
h_sig.SetLineColor(1)
h_sig.Fit("gaus")
print(h_sig.Integral())
#h_sig.GetYaxis().SetRangeUser(0.01,h_sig.GetMaximum()*1.2)
c_sig.Update()

input()        
