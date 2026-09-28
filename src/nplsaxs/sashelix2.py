#!/usr/bin/env python
import math
import cmath
import numpy as np
from scipy.integrate import quad, dblquad, nquad
from scipy.optimize import curve_fit
from scipy.special import jv


def sas_helical_ribbon(width, radius, pitch, q):
    pattern = np.zeros_like(q)
    tan_psi = pitch/(2*np.pi*radius)
    delta = width*math.sqrt(4*np.pi**2*radius**2+pitch**2)/(2*np.pi*radius)
    eps = np.finfo(np.float64).eps

    def func_intensity(theta, chi, q, radius, delta, N):
        ct = math.cos(theta)
        st = math.sin(theta)
        qct = q*ct
        qrct = radius*qct
        qrst = q*radius*st
        jn = jv(N, qrst)
        c = 0.5*delta*qct
        a = delta*math.sin(c)/c if c else delta
        Fn = 2.0*np.pi*jn*a * cmath.exp((c+N*chi)*1.0j) * 1.0j**N
        absFn = abs(Fn)
        Iq = absFn*absFn*st
        return Iq


    for i in range(q.size):
        if i%10 == 0:
            print(f"i = {i}")
        Nmax = math.ceil(pitch*q[-1]/np.pi)
        for N in range(-Nmax, Nmax+1):
            y, _  = dblquad(func_intensity, 0.0, np.pi, 0.0, np.pi/2,
                           args=(q[i], radius, delta, N)
                           )
            pattern[i] += y
    return q, pattern

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    figh, axh = plt.subplots(nrows=1, ncols=2, figsize=(12,6))

    qbeg = 0.01
    qend = 2.0
    qstep = 0.01
    qvals = np.linspace(qbeg, qend, 200)

    radius = 20.0
    pitch = 64 #2*np.pi*radius*np.tan(np.radians(27.0))
    width = 15.0

    #for radius in [15, 30, 45, 60]:
    for pitch in [40, 80, 120]:
        q, pattern = sas_helical_ribbon(width, radius, pitch, qvals)
        with open(f"sas2_r_{radius}.txt", 'w') as fh:
            for i in range(q.size):
                fh.write(f"{q[i]}  {pattern[i]}\n")
        
        axh[0].loglog(q, pattern, ls='-', marker='None',
                      label=rf"$R$={radius}nm, $P$={pitch}nm")
        axh[1].semilogx(q, pattern*q**2, ls='-', marker='None',
                        label=rf"R={radius}nm, $P$={pitch}nm")
    axh[0].legend(loc='best')
    axh[1].legend(loc='best')
    axh[0].set_xlabel(r'$q\ (\mathrm{nm}^{-1})$')
    axh[0].set_ylabel(r'$I(q)$')
    axh[1].set_xlabel(r'$q\ (\mathrm{nm}^{-1})$')
    axh[1].set_ylabel(r'$q^2 \times I(q)$')
    plt.show()
