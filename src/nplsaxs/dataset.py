import os
import sys
import math
from pathlib import Path
import numpy as np
from scipy import stats
from .nanoplatelet import Nanoplatelet


def get_rvs(dist, m, pd, rng, nrv=1, bounds=None, out=None):
    """
    Samples a univariate probability distribution, possibly truncated.

    Parameters
    ----------
    dist : `'uniform'` | `'normal'` | `'lognormal'`
        Normal of the distribution.
    m : float
        Mean or median (for lognormal)
    pd : float
        Polydispersity. Must be positive for all cases. For lognormal must be >
        1.
    rng : :py:class:`numpy.random.Generator`
        A random number generator.
    nrv : int
        Number of sampled values.
    bounds : tuple
        A tuple of the form (`lb`, `rb`), where `lb` (`ub`) represents the
        lower (upper) bound on the sampled values. `ub` can be infinity.
        Applicable only for `normal` or `lognormal` distributions.
    out : (`nrv`,) ndarray
        If given, contains the sampled values on return.

    Returns
    -------
    (`nrv`,) ndarray
        The sampled values. If `nrv`=1, will be a scalar.

    """
    if out is None:
        out = np.zeros((nrv,), dtype=np.float64)
    if out.size != nrv:
        raise ValueError(f"out must have size = {nrv}.")

    if dist == 'uniform':
        #`pd` = sigma/m
        sigma = pd*m
        low = m - sigma
        #If using numpy.random
        #high= m + sigma
        #out[:] = rng.uniform(low=low, high=high, size=nrv)
        #If using scipy.stats
        frzn = stats.uniform(loc=low, scale=2*sigma)
    elif dist == 'normal':
        sigma = pd*m
        #If using numpy.random
        #out[:] = rng.normal(loc=m, scale=sigma, size=nrv)
        #If using scipy.stats
        frzn = stats.norm(loc=m, scale=sigma)
    elif dist == 'lognormal':
        #Here `m` is the geometric mean which equals the median
        #and pd (>1) is the geometrical standard deviation.
        #If using numpy.random
        #out[:] = rng.lognormal(mean=math.log(m), sigma=math.log(pd), size=nrv)
        #If using scipy.stats
        frzn = stats.lognorm(s=math.log(pd), loc=0, scale=m)
    else:
        raise ValueError(f"Unknown value (={dist}) of dist.")
    if bounds is not None:
        cdf_low = frzn.cdf(bounds[0])
        cdf_high = frzn.cdf(bounds[1])
        if cdf_high <= cdf_low:
            raise ValueError("Cannot sample as CDF(high) <= CDF(low).")
        u = rng.uniform(low=cdf_low, high=cdf_high, size=nrv)
        out[:] = frzn.ppf(u)
    else:
        out[:] = frzn.rvs(size=nrv, random_state=rng)
    return np.squeeze(out)


def sample_rp(length, width, radius, pitch, rng, dist_pars,
                        r_vals, p_vals):
    """
    Draws a set of radius and pitch values.

    Parameters
    ----------
    length : float
        Length of the nanoplatelet.
    width : float
        Width of the nanoplatelet.
    radius : tuple
        Specification for sampling radius.
    pitch : tuple
        Specification for sampling pitch.
    rng : :py:class:`numpy.random.Generator`
        A random number generator.
    dist_pars : (4,) ndarray
        On return, contains the parameters of the distribution of radius
        and pitch.
    r_vals : (n,) ndarray
        On return, contains the sampled values of radius.
    p_vals : (n,) ndarray
        On return, contains the sampled values of pitch.

    Returns
    -------
    None

    """
    assert r_vals.size == p_vals.size
    #Sampling radius
    dist = radius[0]
    m = rng.uniform(radius[1], radius[2]) #Mean or median (for lognormal)
    #Polydispersity
    if dist == 'lognormal':
        pd = rng.uniform(1.0, radius[3])
    else:
        pd = rng.uniform(0.0, radius[3])
    dist_pars[0] = m
    dist_pars[1] = pd
    if dist == 'uniform':
        bounds = None
    else:
        bounds = (width/(2*np.pi), np.inf)
    n = r_vals.size
    get_rvs(dist, m, pd, rng, n, bounds=bounds, out=r_vals)

    #Sampling pitch
    dist = pitch[0]
    m = rng.uniform(pitch[1], pitch[2])
    if dist == 'lognormal':
        pd = rng.uniform(1.0, pitch[3])
    else:
        pd = rng.uniform(0.0, pitch[3])
    dist_pars[2] = m
    dist_pars[3] = pd
    if dist == 'uniform':
        get_rvs(dist, m, pd, rng, n, out=p_vals)
    else:
        for i in range(n):
            ri = r_vals[i]
            if 2*np.pi*ri > length:
                p_vals[i] = get_rvs(dist, m, pd, rng, 1)
            else:
                s = 4*np.pi*np.pi*ri*ri - width*width
                delta = 2*np.pi*ri*width/math.sqrt(s)
                bounds = (delta, np.inf)
                p_vals[i] = get_rvs(dist, m, pd, rng, 1, bounds=bounds)


def create(fn_out='out.npz', length=None, width=None, nlayers=None,
            radius=None, pitch=None, nsamp=1, phi=0.001, npart=128,
            calculator='AESDebye', pattern_type='x', Qbeg=0.01,
            Qend=1.0, Qstep=0.01, nthreads=-1, ncells=15,
            debyer_cmd='debyer'):
    npl = Nanoplatelet()
    npl.set_xtal_unit_cell('CdSe', 'zincblende', 6.08, 6.08, 6.08)
    npl.set_xtal_extents(length, width, None, bc='nnn', nlayers=nlayers)
    #Volume of a single particle
    npl_vol = npl.xtal_length * npl.xtal_width * npl.xtal_thickness

    rng = np.random.default_rng()
    #Distribution parameters. For each sample there are four parameters:
    #(1) median of radius (2) geometric standard derivation of radius
    #(3) median of pitch (4) geometric standard derivation of pitch
    dist_pars = np.zeros((nsamp,4), dtype=np.float64)
    r_vals = np.zeros((npart,), dtype=np.float64)
    p_vals = np.zeros((npart,), dtype=np.float64)
    Q_sample = None
    pattern_sample= None
    for isamp in range(nsamp):
        if (isamp+1)%10 == 0:
            print(f"isamp = {isamp}")
        sample_rp(length, width, radius, pitch, rng, dist_pars[isamp,:],
                  r_vals, p_vals)
        #print(f"Distribution parameters: {dist_pars[isamp,:]}")
        for ipart in range(npart):
            #print(f"  ipart = {ipart}")
            #print(f"  R={r_vals[ipart]}, P={p_vals[ipart]}")
            if p_vals[ipart] == 0.0:
                npl.set_shape('Cylinder', radius=r_vals[ipart])
            else:
                npl.set_shape('HelicalRibbon', radius=r_vals[ipart],
                              pitch=p_vals[ipart])
            npl.create(orient_along=None)
            #Calculate SAXS pattern
            Q, pattern = npl.calc_scattering_pattern(
                                fn=None, calculator=calculator,
                                pattern_type=pattern_type,
                                beg=Qbeg, end=Qend, step=Qstep,
                                nthreads=nthreads, ncells=ncells,
                                debyer_cmd=debyer_cmd)
            if Q_sample is None:
                Q_sample = Q.copy()
            if pattern_sample is None:
                pattern_sample = np.zeros((nsamp, Q_sample.size),
                                          dtype=np.float64)
            pattern_sample[isamp,:] += pattern

    pattern_sample *= (phi/(npart*npl_vol*1e-3))
    dist_pars[:,0] /= 10
    dist_pars[:,2] /= 10
    np.savez_compressed(fn_out, length=length/10, width=width/10,
                        dist_pars=dist_pars, Q=Q_sample*10,
                        pattern=pattern_sample)
