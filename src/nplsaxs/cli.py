#!/usr/bin/env python

import argparse
import sys
import glob
from pathlib import Path
import numpy as np
from .dataset import create

parser = argparse.ArgumentParser(prog='saxsds', 
                formatter_class=argparse.ArgumentDefaultsHelpFormatter
                )

subparsers = parser.add_subparsers(dest='subcommand')

sp = subparsers.add_parser('collate',
                    help='combines multiple dataset files to a'
                    ' single file.')
sp.add_argument('glob_pattern',
                   type=str,
                   help="a glob pattern for the files to combine."
                   )
sp.add_argument('--coutfile',
                   type=str,
                   default='combined.npz',
                   help="Name of the combined file (with .npz extension)."
                   )

sp = subparsers.add_parser('create',
                    help='creates a dataset.')
sp.add_argument('length',
                    type=float, 
                    help="Length of the nanoplatelet (nm).")
sp.add_argument('width',
                    type=float, 
                    help="Width of the nanoplatelet (nm).")
sp.add_argument('ML',
                    type=int,
                    help="Number of monolayers. e.g. 1 ML CdSe means one layer"
                    " of Cd atoms + one layer of Se atoms + one layer of"
                    " Cd atoms.")
sp.add_argument('--rdist',
                    type=str,
                    choices=['normal', 'lognormal', 'uniform'],
                    default='lognormal',
                    help="Distribution of radius.")
sp.add_argument('--rmin',
                    type=float,
                    default=10.0,
                    help="Lower bound of radius (nm)."
                    )
sp.add_argument('--rmax',
                    type=float,
                    default=20.0,
                    help="Upper bound of radius (nm)."
                    )
sp.add_argument('--rpd',
                    type=float,
                    default=1.1,
                    help="Maximum polydispersity in radius. Must be > 1 for"
                    " lognormal distribution, must be > 0 for all other cases."
                    )
sp.add_argument('--pdist',
                    type=str,
                    choices=['normal', 'lognormal', 'uniform'],
                    default='lognormal',
                    help="Distribution of pitch."
                    )
sp.add_argument('--pmin',
                    type=float,
                    default=20.0,
                    help="Lower bound of pitch (nm)."
                    )
sp.add_argument('--pmax',
                    type=float,
                    default=40.0,
                    help="Upper bound of pitch (nm)."
                    )
sp.add_argument('--ppd',
                    type=float,
                    default=1.1,
                    help="Maximum polydispersity in pitch. Must be > 1 for"
                    " lognormal distribution, must be > 0 for all other cases."
                    )
sp.add_argument('--outfile',
                    type=str,
                    default='ds.npz',
                    help="Name of the output file. Must have extension '.npz'.")
sp.add_argument('--nsamp',
                    type=int,
                    default=10,
                    help="Number of samples.")
sp.add_argument('--phi',
                    type=float,
                    default=0.001,
                    help="Volume fraction.")
sp.add_argument('--npart',
                    type=int,
                    default=128,
                    help="Number of nanoplatelets in each sample.")
sp.add_argument('--calculator',
                    type=str,
                    choices=['AESDebye', 'Debyer'],
                    default='AESDebye',
                    help="Program to use for calculating scattering patterns.")
sp.add_argument('--pattern_type',
                    type=str,
                    choices=['x', 'S'],
                    default='x',
                    help="'x' for X-ray scatttering, 'S' for total scattering"
                    " structure function (scattering factor).")
sp.add_argument('--Qbeg',
                    type=float,
                    default=0.01,
                    help="Start of scattering pattern (1/nm).")
sp.add_argument('--Qend',
                    type=float,
                    default=10.0,
                    help="End of scattering pattern (1/nm).")
sp.add_argument('--Qstep',
                    type=float,
                    default=0.01,
                    help="Q spacing (1/nm).")
sp.add_argument('--nthreads',
                    type=int,
                    default=-1,
                    help="Number of threads to use for AESDebye. -1 indicates"
                    " all available threads.")
sp.add_argument('--ncells',
                    type=int,
                    default=15,
                    help="Number of cells to use for AESDebye.")
sp.add_argument('--debyer_cmd',
                    type=str,
                    default='debyer',
                    help="If using Debyer, the name of the executable.")

def run():
    args = parser.parse_args()
    
    if args.subcommand == 'collate':
        fns = glob.glob(args.glob_pattern)
        nfiles = len(fns)
        if nfiles == 0:
            print(f"No files matching {args.glob_pattern} found. Exiting ...")
            raise SystemExit()
        else:
            print(f"Number of files = {nfiles}")
            #Open the first file
            with np.load(fns[0]) as fh:
                length = fh['length']
                print(f"length = {length}")
                width = fh['width']
                print(f"width = {width}")
                Q = fh['Q']
                print(f"Q vector dimension = {Q.size}")
                nsamp = fh['pattern'].shape[0]
                print(f"Number of samples per file = {nsamp}")
                pattern = np.zeros((nsamp*nfiles, Q.size), dtype=np.float64)
                dist_pars = np.zeros((nsamp*nfiles,fh['dist_pars'].shape[1]),
                                     dtype=np.float64)
                print(f"Number of distribution parameters per sample"
                      f" = {dist_pars.shape[1]}.")
            for i, fn in enumerate(fns):
                with np.load(fn) as fh:
                    ibeg = i*nsamp
                    iend = (i+1)*nsamp
                    dist_pars[ibeg:iend,:] = fh['dist_pars']
                    pattern[ibeg:iend,:] = fh['pattern']
            fn_out = args.coutfile
            np.savez_compressed(fn_out, length=length, width=width,
                            dist_pars=dist_pars, Q=Q, pattern=pattern)
    elif args.subcommand == 'create':
        length = args.length*10 #From nm to angstrom
        width = args.width*10 #From nm to angstrom
        nlayers = 2*args.ML + 1
        radius = (args.rdist, args.rmin*10, args.rmax*10, args.rpd)
        pitch = (args.pdist, args.pmin*10, args.pmax*10, args.ppd)
        nsamp = args.nsamp
        phi = args.phi
        npart = args.npart
        calculator = args.calculator
        Qbeg = args.Qbeg/10 #From 1/nm to 1/angstrom
        Qend = args.Qend/10 #From 1/nm to 1/angstrom
        Qstep = args.Qstep/10 #From 1/nm to 1/angstrom
        fn_out = args.outfile
        pattern_type = args.pattern_type
        nthreads = args.nthreads
        ncells = args.ncells
        debyer_cmd = args.debyer_cmd
        
        create(fn_out=fn_out, length=length, width=width, nlayers=nlayers,
               radius=radius, pitch=pitch, nsamp=nsamp, phi=phi,
               npart=npart, calculator=calculator,
               pattern_type=pattern_type, Qbeg=Qbeg, Qend=Qend,
               Qstep=Qstep, nthreads=nthreads, ncells=ncells,
               debyer_cmd=debyer_cmd)
