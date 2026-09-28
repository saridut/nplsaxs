#!/usr/bin/env python

from pathlib import Path
import sys
import math
import numpy as np
from scipy import stats
import pytest
#try:
#    from nplsaxs import dataset
#except ModuleNotFoundError:
srcdir = Path('src').resolve()
sys.path.insert(0, str(srcdir))
from nplsaxs import dataset, saxsds
import matplotlib.pyplot as plt

#init_params = [
#        (0, 4, 0, 1, 1, 1),
#        (2, 4, 0, 1, 1, 1),
#        (20, 4, 0, None, 1, 1, {'ngpl': None}),
#        (20, 0, 0, 1, 1, 1),
#        (20, 2, 0, 1, 1, 1),
#        (20, 4, 0, 1, None, 1, {'ngpw': None}),
#        (20, 4, -1, 1, 1, 1),
#        (20, 4, 1, 1, 1, None, {'ngpt': 1}),
#        (20, 4, 1, 1, 1, None, {'ngpt': None}),
#        ]
#
#@pytest.fixture(scope='module', params=init_params)
#def init_args(request):
#    par = request.param
#    if isinstance(par[-1], dict):
#        args = par[:-1]; kwargs = par[-1]
#    else:
#        args = par; kwargs = {}
#    return args, kwargs


#def test_init_error(init_args):
#    with pytest.raises(ValueError) as info:
#        args = init_args[0]; kwargs = init_args[1]
#        ribbon = Ribbon(*args, **kwargs)
#    assert info.type is ValueError

def test_rvs():
    rng = np.random.default_rng()
    dist = 'lognormal'
    m = 30.0
    pd = 1.51
    nrv = 1024
    out = np.zeros((nrv,), dtype=np.float64)
    #bounds = (10.0, np.inf)
    bounds = None
    dataset.get_rvs(dist, m, pd, rng, nrv, bounds, out)
    print(f"median = {np.median(out)}")
    figh, axh = plt.subplots(nrows=1, ncols=1, figsize=(8,6))
    axh.hist(out, bins=100, density=True)
    axh.set_xlim(0, out.max())
    plt.show()
    assert True


def test_sample_rp():
    length = 1000
    width = 100
    radius = ('lognormal', 50, 200, 1.2)
    pitch = ('lognormal', 100, 500, 1.2)
    npart = 512
    nbins = 50
    rng = np.random.default_rng()
    dist_pars = np.zeros((4,), dtype=np.float64)
    r_vals = np.zeros((npart,), dtype=np.float64)
    p_vals = np.zeros_like(r_vals)
    dataset.sample_rp(length, width, radius, pitch, rng, dist_pars, r_vals,
                      p_vals)
    dist_pars[0] /= 10
    dist_pars[2] /= 10
    r_vals /= 10
    p_vals /= 10
    print(f"\nDistribution parameters: {dist_pars}")
    #Fitting
    data = r_vals
    s, loc, scale = stats.lognorm.fit(data, math.log(2.0),
                                      scale=np.median(data), floc=0)
    print(f"R fit: m={scale}, pd={math.exp(s)}, loc={loc}")
    data = p_vals
    s, loc, scale = stats.lognorm.fit(data, math.log(2.0),
                                      scale=np.median(data), floc=0)
    print(f"P fit: m={scale}, pd={math.exp(s)}, loc={loc}")

    figh, axh = plt.subplots(nrows=1, ncols=1, figsize=(8,6))
    axh.hist(r_vals, bins=nbins, histtype='step', align='mid', density=True,
             label='R')
    axh.hist(p_vals, bins=nbins, histtype='step', align='mid', density=True,
             label='P')
    axh.legend(loc='best')
    #axh.set_xlim(0, out.max())
    plt.show()
    assert True


def test_create_dataset():
    fn_out = f"out.npz"
    length = 1000
    width = 100
    nlayers = 7
    radius = ('lognormal', 50, 200, 1.2)
    pitch = ('lognormal', 100, 500, 1.2)
    nsamp = 2
    phi = 0.001
    npart = 64
    calculator = 'AESDebye'
    pattern_type = 'x'
    Qbeg = 0.001
    Qend = 8.0
    Qstep = 0.002
    nthreads = -1
    ncells = 15
    debyer_cmd = 'debyer' 

    dataset.create(fn_out=fn_out, length=length, width=width, nlayers=nlayers,
                   radius=radius, pitch=pitch, nsamp=nsamp, phi=phi,
                   npart=npart, calculator=calculator,
                   pattern_type=pattern_type, Qbeg=Qbeg, Qend=Qend,
                   Qstep=Qstep, nthreads=nthreads, ncells=ncells,
                   debyer_cmd=debyer_cmd)
    assert True


#def test_saxsds():
#    saxs
