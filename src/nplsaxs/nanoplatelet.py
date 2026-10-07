"""
Class implementing a single NPL crystal.

"""

import os
import sys
import copy
import math
from pathlib import Path
import subprocess
import numpy as np
import tempfile
import ase 
import ase.build
import ase.io
import ase.geometry
import aesdebye
from rotlib import aa_rotate_vectors
from ribbon import Ribbon


class Nanoplatelet(object):
    """
    Attributes
    ----------
    bulk_uc : ase.atoms.Atoms
        Unit cell of the bulk crystal.
    xtal_bc : str
        Boundary conditions of the crystal. 
    xtal_length : float
        Length of the crystal
    xtal_width : float
        Width of the crystal
    xtal_thickness : float
        Width of the crystal
    xtal_nlayers : int
        Number of layers in the normal direction
    xtal_cut_box_params : dict
        Parameters for the cutting box
    radius : float
        Radius
    pitch : float
        Pitch
    kg : float
        Gauss curvature
    km : float
        Mean curvature
    theta : float
        Angle between the principal curvature direction and the length-wise
        direction of the crystal.
    origin : (3,) ndarray
        Origin of the nanoplatelet
    orientation : (4,) ndarray
        Quaternion giving the orientation of the nanoplatelet with respect to
        the world frame.
    origo : (3,) ndarray
        Origin of the bounding cell

    """
    def __init__(self):
        """
        Initializes an empty nanoplatelet.

        """
        self.origin = np.zeros((3,), dtype=np.float64) 
        self.origo = np.zeros((3,), dtype=np.float64) 
        #self.origo is updated after call to self.create()
        self.orientation = np.array([1, 0, 0, 0], dtype=np.float64)
        self.ribbon = None #empty ribbon
        self._atoms_ref = ase.Atoms()
        self._atoms = ase.Atoms()


    @property
    def xtal_length(self):
        return self.ribbon.length


    @property
    def xtal_width(self):
        return self.ribbon.width


    @property
    def xtal_thickness(self):
        return self.ribbon.thickness


    @property
    def radius(self):
        return self.ribbon.get_radius()


    @property
    def pitch(self):
        return self.ribbon.get_pitch()


    @property
    def kg(self):
        return self.ribbon.get_gauss_curvature()


    @property
    def km(self):
        return self.ribbon.get_mean_curvature()


    @property
    def theta(self):
        return self.ribbon.get_theta()


    def set_xtal_unit_cell(self, name, lattice, a, b=None, c=None, alpha=None,
                           covera=None, u=None):
        """
        Setter for the crystal unit cell.

        Parameters
        ----------
        name : str
            Chemical symbol or symbols (at most two atoms), e.g. 'CdSe'.
        lattice : `'sc'` | `'fcc'` | `'bcc'` | `'orthorhombic'` |
            `'zinclblende'` | `'wurtzite'`
            Crytal lattice type.
        a : float
            Lattice parameter (in nm)
        b : float, optional
            Lattice parameter (in nm)
        c : float, optional
            Lattice parameter (in nm)
        alpha : float, optional
            Angle for rhombohedral unit cells
        covera : float, optional
            c/a ratio used for hcp. Default is ideal ratio: sqrt(8/3)
        u: float, optional
            Internal coordinate for Wurtzite structure.

        """
        self.bulk_uc = ase.build.bulk(name, lattice, a, b, c, alpha=alpha,
                                      covera=covera, u=u, orthorhombic=False,
                                      cubic=True)


    def set_xtal_extents(self, length, width, thickness, bc='nnn',
                         nlayers=None, locns=0):
        """
        Sets the extents and boundary conditions of the crystal.

        Parameters
        ----------
        bc : str
            A three character string specifying the boundary conditions along
            X, Y, and Z directions, respectively, of the fixed external frame.
            Each of the characters can be either `'p'`, `'c'`, or `'n'`,
            signifying periodic, end-capped, or non-periodic boundary
            conditions.
        length : float
            Length of the crystal in nm
        width : float
            Width of the crystal in nm
        thickness : float | None
            Thickness of the crystal in nm
        nlayers : int | None
            Number of layers along the [0,0,1] direction.
        locns : float
            Fractional distance of the neutral surface from the midsurface.
            Must be <= 0.5 and >= -0.5. A negative (positive) value indicates
            that the neutral surface is below (above) the midsurface. If there
            is a single layer `locns` is ignored.

        """
        for i in range(3):
            if bc[i] not in ['p', 'c', 'n']:
                raise ValueError(
                    f"bc[{i}] (= {bc[i]}) must be one of 'p', 'c', or 'n'." )
        self.xtal_bc = bc

        Nx = 1; Ny = 1; Nz = 1
        extend = [1.0, 1.0, 1.0]

        if thickness is not None:
            if bc[2] == 'p':
                Nz = math.ceil(thickness/self.bulk_uc.cell[2,2])
            elif bc[2] == 'c':
                Nz = math.ceil(thickness/self.bulk_uc.cell[2,2])
                extend[2] = 1.01
            elif bc[2] == 'n':
                Nz = thickness/self.bulk_uc.cell[2,2]
            xtal_thickness = self.bulk_uc.cell[2] * Nz
            #Find the number of layers
            atoms = ase.build.cut(self.bulk_uc, c=(0,0,self.bulk_uc.cell[2,2]),
                                  clength=xtal_thickness, extend=extend[2])
            tags, levels = ase.geometry.get_layers(atoms, [0,0,1])
            self.nlayers = levels.size
        else:
            self.nlayers = nlayers
            if self.nlayers == 1:
                xtal_thickness = 0.0
            else:
                atoms = ase.build.cut(self.bulk_uc, c=(0,0,self.bulk_uc.cell[2,2]), 
                            nlayers=self.nlayers, tolerance=0.001)
                xtal_thickness = atoms.positions[:,2].max() \
                                - atoms.positions[:,2].min()

        if bc[0] == 'p':
            Nx = math.ceil(length/self.bulk_uc.cell[0,0])
        elif bc[0] == 'c':
            Nx = math.ceil(length/self.bulk_uc.cell[0,0])
            extend[0] = 1.01
        elif bc[0] == 'n':
            Nx = length/self.bulk_uc.cell[0,0]
        xtal_length = self.bulk_uc.cell[0,0] * Nx

        if bc[1] == 'p':
            Ny = math.ceil(width/self.bulk_uc.cell[1,1])
        elif bc[1] == 'c':
            Ny = math.ceil(width/self.bulk_uc.cell[1,1])
            extend[1] = 1.01
        elif bc[1] == 'n':
            Ny = width/self.bulk_uc.cell[1,1]
        xtal_width = self.bulk_uc.cell[1,1] * Ny

        self.xtal_cut_box_params = {'vectors': np.diag([Nx, Ny, Nz]),
                                    'extend': extend.copy(),
                                    'nlayers': nlayers}
        self.ribbon = Ribbon(xtal_length, xtal_width, xtal_thickness,
                             self.bulk_uc.cell[0,0], self.bulk_uc.cell[1,1],
                             None, ngpt=2, locns=locns)


    def set_shape(self, shape, angle=None, **params):
        """
        Sets the shape of the nanoplatelet.

        Parameters
        ----------
        shape : str
            Specifies the shape of the crystal in its current configuration.
        angle : None | float
            Pass angle = 0 for no rotation.
        params : dict
            Shape parameters.

        Returns
        -------
        None

        """
        self.ribbon.set_shape(shape, **params)

        if angle is None:
            theta = self.ribbon.get_theta()
            rotate_by = theta - np.pi/4
        else:
            rotate_by = angle - np.pi/4

        a = self.xtal_cut_box_params['vectors'][0]
        b = self.xtal_cut_box_params['vectors'][1]
        c = self.xtal_cut_box_params['vectors'][2]
        nlayers = self.xtal_cut_box_params['nlayers']
        extend = self.xtal_cut_box_params['extend']
        if np.isclose(rotate_by, 0.0, 1e-8, 1e-14):
            #No rotation necessary
            self._atoms_ref = ase.build.cut(self.bulk_uc, a=a, b=b, c=c,
                                            nlayers=nlayers, extend=extend,
                                            tolerance=0.001)
        else:
            if (self.xtal_bc[0]=='p' or self.xtal_bc[0]=='c' or
                self.xtal_bc[1]=='p' or self.xtal_bc[1]=='c'):
                raise ValueError("Crystal rotation not allowed with a periodic"
                                " boundary along X or Y direction.")
            zhat = np.array([0,0,1], dtype=np.float64)
            a = aa_rotate_vectors(a, zhat, rotate_by)
            b = aa_rotate_vectors(b, zhat, rotate_by)
            self._atoms_ref = ase.build.cut(self.bulk_uc, a=a, b=b, c=c,
                                            nlayers=nlayers, extend=extend,
                                            tolerance=0.001)
            self._atoms_ref.rotate(-rotate_by, 'z', rotate_cell=True)

        if nlayers is not None: 
            zmin = self._atoms_ref.positions[:,2].min()
            self._atoms_ref.positions[:,2] -= zmin
            zmax = self._atoms_ref.positions[:,2].max()
            self._atoms_ref.cell[2,2] = zmax

        #Ensure all atoms are within bounds
        scaled_pos = self._atoms_ref.get_scaled_positions(wrap=False)
        self._atoms_ref.cell = np.diag([self.xtal_length, self.xtal_width,
                                        self.xtal_thickness])
        self._atoms_ref.set_scaled_positions(scaled_pos)
        #Shift for ribbon. Atoms are out of cell after this.
        self._atoms_ref.positions[:,1] -= 0.5*self.xtal_width
        self._atoms_ref.positions[:,2] -= 0.5*self.xtal_thickness

        self._atoms = self._atoms_ref.copy()
        self.ribbon.set_atom_refpos(self._atoms_ref.positions, copy=False,
                                    atom_pos=self._atoms.positions)


    def create(self, orient_along):
        """
        Create the crystal.

        Parameters
        ----------
        orient_along : (3,) array_like
            Direction to orient the crystal along.

        Returns
        -------
        None

        """
        self.ribbon.create(orient_along)
        self.ribbon.translate_to_center()
        pmax = self.ribbon.grid.max(axis=(0,1,2), keepdims=False) 
        pmin = self.ribbon.grid.min(axis=(0,1,2), keepdims=False) 
        self._atoms.cell = np.diag(pmax-pmin)
        self.origo = pmin.copy()


    def get_mesh(self, tag=None):
        """
        Returns the grid with respect to :attr:`.origin`.

        Parameters
        ----------
        tag : None | `'ml'` | `'ms'`
            Name of the file to write out the atom positions.

        Returns
        -------
        ndarray

        """
        if tag is None:
            return self.ribbon.grid
        elif tag == 'ml':
            return self.ribbon.mline
        elif tag == 'ms':
            return self.ribbon.msurf


    def get_atom_positions(self):
        """
        Returns the atom positions with respect to :attr:`.origin`.

        """
        return self._atoms.positions


    def to_xyz(self, fn, title=''):
        """
        Write to a XYZ file.

        Parameters
        ----------
        fn : pathlike
            Name of the file to write out the atom positions.

        """
        self._atoms.positions -= self.origo
        #assert np.all(self._atoms.positions.min(axis=0) >= 0.0)
        #assert np.all(self._atoms.positions.max(axis=0) 
        #              <= self._atoms.cell.lengths())
        ase.io.write(fn, self._atoms, format='xyz', comment='', fmt='%22.15f')
        self._atoms.positions += self.origo
        

    def as_hoomd_convex_polyhedron_union(self):
        raise NotImplementedError


    def calc_scattering_pattern(self, fn=None, calculator='AESDebye',
                                pattern_type='x', beg=0.01, end=1.0,
                                step=0.01, nthreads=-1, ncells=15,
                                debyer_cmd='debyer'):
        """
        Calculates a scattering pattern.

        Parameters
        ----------
        fn : pathlike | None
            Name of the file to write out the scattering pattern. If ``None``,
            no output will be written.
        calculator : `'AESDebye'` | `'Debyer'`
            Name of the calculator to use.
        pattern_type : `'x'` | `'S'`
            `'x'` for X-ray scattering pattern, `'S'` for structure pattern.
        beg : float
            Starting value of Q
        end : float
            Ending value of Q
        step : float
            Step size in Q
        nthreads : int
            Number of OpenMP threads to be used by _AESDebye_
        ncells : int
            Number of cells in th lattice for _AESDebye_
        debyer_cmd : pathlike
            Path of the _debyer_ command.

        """
        if calculator == 'AESDebye':
            if pattern_type == 'S':
                symbols = ['None' for x in self._atoms.get_chemical_symbols()]
            elif pattern_type == 'x':
                symbols = self._atoms.get_chemical_symbols()
            positions = aesdebye.Positions(chemicalSymbols=symbols,
                                           coordinates=self._atoms.positions)
            dc = aesdebye.DebyeCalculator(nThreads=nthreads, nCells=ncells,
                                          useMPI=False, useGPU=False,
                                          verbose=False, binsResolution=0.001)
            results = dc.calculateProfile(positions, start=beg, end=end+step,
                                          steps=1+int((end-beg)/step))
            pdf, profile = results['Cd-Cd']
            Q = np.array(profile.q, copy=None)
            pattern_CdCd = np.array(profile.intensity, copy=None)

            pdf, profile = results['Se-Se']
            pattern_SeSe = np.array(profile.intensity, copy=None)

            pdf, profile = results['Cd-Se']
            pattern_CdSe = np.array(profile.intensity, copy=None)

            pattern = pattern_CdCd + pattern_SeSe + 2*pattern_CdSe

            #pdf, profile = results['total']
            #Q = np.array(profile.q, copy=None)
            #pattern = np.array(profile.intensity, copy=None)

            if fn is not None:
                profile.toCSV(fn)

        elif calculator == 'Debyer':
            with tempfile.TemporaryDirectory() as tmpdir:
                fn_xyz = Path(tmpdir, 'tmp.xyz')
                self.to_xyz(fn_xyz)
                fn_out = Path(tmpdir, 'out.txt') if fn is None else fn

                args = [debyer_cmd, f"-{pattern_type}", "-a0", "-b0", "-c0",
                        f"-f{beg}", f"-t{end}", f"-s{step}", f"-o{fn_out}",
                        f"{fn_xyz}"]
                cp = subprocess.run(args, capture_output=True, text=True,
                                shell=False, check=True)
                print(cp.stdout)
                if 'ERROR' in cp.stdout:
                    raise SystemExit()
                
                with open(f"{fn_out}", 'r') as fh:
                    lines = fh.readlines()
                Q = []
                pattern = []
                for line in lines:
                    stripped = line.strip(' \n')
                    if stripped.startswith('#'):
                        continue
                    else:
                        words = stripped.split()
                        Q.append(float(words[0]))
                        pattern.append(float(words[1]))
                Q = np.asarray(Q)
                pattern = np.asarray(pattern)
        else:
            raise ValueError(f"calculator (={calculator}) must be 'AESDebye'"
                             " or 'Debyer'.")
        return Q, pattern
