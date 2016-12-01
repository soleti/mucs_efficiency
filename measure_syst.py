#!/usr/bin/env python3.4

import math



def measure_sig(val1,val2):
    return (val1-val2)/math.sqrt(err1**2+err2**2)
    
theta_angle = {7: 90, 8: 105, 5: 60, 6: 75}
phi_angle = {7: -90, 8: -75, 9: -60}
l_dict = {1: 20, 2: 80, 3: 140, 4: 200, 5: 260}
    
sys = 0
sys_value = 0.1

f_3d = open("output/sys_errors_3d.txt","w")

with open("output/central_shift.txt") as textfile1, open("output/upstream_shift.txt") as textfile2, open("output/downstream_shift.txt") as textfile3, open("output/all.txt") as all: 
    for x, y, z, m in zip(textfile1, textfile2, textfile3, all):
        a = x.split(" ")
        a = [float(i) for i in a]
        b = y.split(" ")
        b = [float(i) for i in b]
        c = z.split(" ")
        c = [float(i) for i in c]
        m = m.split(" ")
        m = [float(i) for i in m]

        eff = [a[3],b[3],c[3]]
        err = [a[4],b[4],c[4]]
        theta = a[0]
        phi = a[1]
        l = a[2]


        print("%i %i %i %f" % (theta,phi,l, max(eff)-m[3]), file=f_3d)

f_3d.close()

f_theta = open("output/sys_errors_theta.txt","w")

with open("output/theta_central_shift.txt") as textfile1, open("output/theta_upstream_shift.txt") as textfile2, open("output/theta_downstream_shift.txt") as textfile3, open("output/theta_all.txt") as all: 
    for x, y, z, m in zip(textfile1, textfile2, textfile3, all):
        a = x.split(" ")
        a = [float(i) for i in a]
        b = y.split(" ")
        b = [float(i) for i in b]
        c = z.split(" ")
        c = [float(i) for i in c]
        m = m.split(" ")
        m = [float(i) for i in m]

        eff = [a[1],b[1],c[1]]
        err = [a[2],b[2],c[2]]
        theta = a[0]

        print("%i %f" % (theta, max(eff)-m[1]), file=f_theta)

f_theta.close()

f_phi = open("output/sys_errors_phi.txt","w")

with open("output/phi_central_shift.txt") as textfile1, open("output/phi_upstream_shift.txt") as textfile2, open("output/phi_downstream_shift.txt") as textfile3, open("output/phi_all.txt") as all: 
    for x, y, z, m in zip(textfile1, textfile2, textfile3, all):
        a = x.split(" ")
        a = [float(i) for i in a]
        b = y.split(" ")
        b = [float(i) for i in b]
        c = z.split(" ")
        c = [float(i) for i in c]
        m = m.split(" ")
        m = [float(i) for i in m]

        eff = [a[1],b[1],c[1]]
        err = [a[2],b[2],c[2]]
        phi = a[0]

        print("%i %f" % (phi, max(eff)-m[1]), file=f_phi)

f_phi.close()

f_l = open("output/sys_errors_l.txt","w")

with open("output/l_central_shift.txt") as textfile1, open("output/l_upstream_shift.txt") as textfile2, open("output/l_downstream_shift.txt") as textfile3, open("output/l_all.txt") as all: 
    for x, y, z, m in zip(textfile1, textfile2, textfile3, all):
    	a = x.split(" ")
    	a = [float(i) for i in a]
    	b = y.split(" ")
    	b = [float(i) for i in b]
    	c = z.split(" ")
    	c = [float(i) for i in c]
    	m = m.split(" ")
    	m = [float(i) for i in m]

    	eff = [a[1],b[1],c[1]]
    	err = [a[2],b[2],c[2]]
    	l = a[0]
    	print("%i %f" % (l, max(eff)-m[1]), file=f_l)

f_l.close()

f_theta_phi = open("output/sys_errors_theta_phi.txt","w")

with open("output/theta_phi_central_shift.txt") as textfile1, open("output/theta_phi_upstream_shift.txt") as textfile2, open("output/theta_phi_downstream_shift.txt") as textfile3, open("output/theta_phi_all.txt") as all: 
    for x, y, z, m in zip(textfile1, textfile2, textfile3, all):
        a = x.split(" ")
        a = [float(i) for i in a]
        b = y.split(" ")
        b = [float(i) for i in b]
        c = z.split(" ")
        c = [float(i) for i in c]
        m = m.split(" ")
        m = [float(i) for i in m]

        eff = [a[1],b[1],c[1]]
        err = [a[2],b[2],c[2]]
        theta_phi = a[0]

        print("%i %f" % (theta_phi, max(eff)-m[1]), file=f_theta_phi)

f_theta_phi.close()

f_theta_l = open("output/sys_errors_theta_l.txt","w")

with open("output/theta_l_central_shift.txt") as textfile1, open("output/theta_l_upstream_shift.txt") as textfile2, open("output/theta_l_downstream_shift.txt") as textfile3, open("output/theta_l_all.txt") as all: 
    for x, y, z, m in zip(textfile1, textfile2, textfile3, all):
        a = x.split(" ")
        a = [float(i) for i in a]
        b = y.split(" ")
        b = [float(i) for i in b]
        c = z.split(" ")
        c = [float(i) for i in c]
        m = m.split(" ")
        m = [float(i) for i in m]

        eff = [a[1],b[1],c[1]]
        err = [a[2],b[2],c[2]]
        theta_l = a[0]

        print("%i %f" % (theta_l, max(eff)-m[1]), file=f_theta_l)

f_theta_l.close()

f_phi_l = open("output/sys_errors_phi_l.txt","w")

with open("output/phi_l_central_shift.txt") as textfile1, open("output/phi_l_upstream_shift.txt") as textfile2, open("output/phi_l_downstream_shift.txt") as textfile3, open("output/phi_l_all.txt") as all: 
    for x, y, z, m in zip(textfile1, textfile2, textfile3, all):
        a = x.split(" ")
        a = [float(i) for i in a]
        b = y.split(" ")
        b = [float(i) for i in b]
        c = z.split(" ")
        c = [float(i) for i in c]
        m = m.split(" ")
        m = [float(i) for i in m]

        eff = [a[1],b[1],c[1]]
        err = [a[2],b[2],c[2]]
        phi_l = a[0]

        print("%i %f" % (phi_l, max(eff)-m[1]), file=f_phi_l)

f_phi_l.close()