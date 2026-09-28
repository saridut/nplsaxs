#!/usr/bin/env python
import math
import cmath
import numpy as np
from scipy.integrate import quad, dblquad, nquad
from scipy.optimize import curve_fit
from scipy.special import jv


def sas_helical_ribbon(width, radius, pitch, N, qbeg, qend, qstep):
    n = 1 + int( (qend-qbeg)/qstep )
    q = np.linspace(qbeg, qend, num=n)
    pattern = np.zeros_like(q)
    tan_psi = pitch/(2*np.pi*radius)
    delta = width*math.sqrt(4*np.pi**2*radius**2+pitch**2)/(2*np.pi*radius)

    def func_F(phi, b, c, chi, tan_psi, qrst, qrct, qct):
        a = delta*math.sin(c)/c if c else delta
        X = a*cmath.exp(1.0j*(b*phi+c))
        t = qrst*math.cos(chi-phi) + qrct*phi*tan_psi
        F = cmath.exp(t*1.0j)*X
        #t = a2*phi + a3*math.cos(chi-phi) + c
        #ct = math.cos(t)
        #st = math.sin(t)
        #F = a1*(ct + st*1.0j)
        #F = a1*cmath.exp(t*1.0j)
        return F

    def func_intensity(theta, chi, q, radius, tan_psi, delta, N):
        ct = math.cos(theta)
        st = math.sin(theta)
        qct = q*ct
        qrct = radius*qct
        qrst = q*radius*st
        b = qrct
        c = 0.5*delta*qct
        #a1 = 2.0*math.sin(c)/qct
        #a2 = b + qrct*tan_psi
        #a3 = qrst
        F, err, info  = quad(func_F, 0, 2*np.pi*N, 
                             args=(b, c, chi, tan_psi, qrst, qrct, qct), limit=50,
                  complex_func=True, full_output=1, epsabs=1e-6, epsrel=1e-6)
        Iq = abs(F)**2*st
        return Iq


    for i in range(q.size):
        if i%1 == 0:
            print(f"i = {i}")
        #y, err  = quad(func_intensity_chi, 0.0, np.pi, 
        #               args=(q[i], radius, tan_psi, delta, N),
        #              limit=50)
        y, err, info  = nquad(func_intensity, [(0.0, np.pi/2), (0.0, np.pi)],
                       args=(q[i], radius, tan_psi, delta, N), 
                       opts={'limit':50}, full_output=True
                       )
        pattern[i] = y
    return q, pattern

if __name__ == "__main__":
    radius = 20.0
    pitch = 64 #2*np.pi*radius*np.tan(np.radians(27.0))
    width = 15.0
    N = 10
    #tan_psi = pitch/(2.0*np.pi*radius)
    #delta = width/tan_psi
    #print(f"radius = {radius}\npitch = {pitch}\nwidth={width}\ndelta={delta}")
    qbeg = 0.01
    qend = 0.2
    qstep = 0.01
    q, pattern = sas_helical_ribbon(width, radius, pitch, N, qbeg, qend, qstep)
    with open(f"sas2f_{N}.txt", 'w') as fh:
        for i in range(q.size):
            fh.write(f"{q[i]}  {pattern[i]}\n")
    import matplotlib.pyplot as plt
    figh, axh = plt.subplots(nrows=1, ncols=2, figsize=(8,3))
    axh[0].loglog(q, pattern, ls='-', marker='None', label="_nolegend_")
    axh[1].semilogx(q, pattern*q**2, ls='-', marker='o', label="_nolegend_")
    axh[0].set_xlabel(r'$q$')
    axh[0].set_ylabel(r'$I(q)$')
    axh[1].set_xlabel(r'$q$')
    axh[1].set_ylabel(r'$q \times 2I(q)$')
    plt.show()
