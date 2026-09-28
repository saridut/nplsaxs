#!/usr/bin/env python

import io
from pathlib import Path

def write(system, file, extended_cfg=True):
    if is_instance(file, io.TextIOBase):
        fh = file
        opened = False
    else:
        fh = open(Path(file), 'w')
        opened = True
    
    fh.write(f"Number of particles = {system.num_atoms}\n")
    fh.write('\n')
    fh.write('A = 1.0 Angstrom\n')
    fh.write('\n')
    for i in [1,2,3]:
        for j in [1,2,3]:
            fh.write(f"H(i,j) {} A\n")
    fh.write('R = 1.0 [ns^-1]\n')
    fh.write('\n')
    if extended_cfg:
        pass
        fh.write('.NO_VELOCITY.\n')
        fh.write('entry_count = 3\n')
    else:
        for iatm in range(cfg.num_atoms):
            fh.write(f"{mass} {name} {x} {y} {z} {vx} {vy} {vz}\n")
    if opened:
        fh.close()


def read(system, fn):
    """
    Extract box vectors from an AtomEye cfg file.

    """
    length_scale = 1.0 #in Angstrom
    rate_scale = 1.0   #in ns^-1
    H0 = np.zeros((3,3)); transform = np.identity(3); eta = np.zeros((3,3))
    atom_section = False
    with open(fn, 'r') as fh:
        while (linein := fh.readline()):
            #Remove blank and comment lines
            line = linein.strip('\n ')
            if len(line)==0 or line.startswith('#'):
                continue
            if line.startswith('Number of particles') :
                N = int(line.split()[4])
                coords = np.empty((N,3)); atom_names = []; iatm = 0
            elif line.startswith('A') and len(words := line.split()) > 2:
                length_scale = float(words[2])
            elif line.startswith('R') and len(words := line.split()) > 2:
                rate_scale = float(words[2])
            elif line.startswith('H0('):
                words = line.split()
                i = int(words[0][3]) - 1; j = int(words[0][5]) - 1
                H0[i,j] = float(words[2])*length_scale
            elif line.startswith('Transform('):
                words = line.split()
                i = int(words[0][10]) - 1; j = int(words[0][12]) - 1
                transform[i,j] = float(words[2])
            elif line.startswith('eta('):
                words = line.split()
                i = int(words[0][4]) - 1; j = int(words[0][6]) - 1
                eta[i,j] = float(words[2])
            elif line.startswith('.NO_VELOCITY.'):
                velocity_data = False
            elif line.startswith('entry_count'):
                entry_count = int(line.split()[2])
                if velocity_data:
                    num_aux = entry_count - 6
                else:
                    num_aux = entry_count - 3
            elif line.startswith('auxiliary'):
                continue
            else:
                #Atom data records
                if not atom_section:
                    atom_section = True
                    H = H0 @ transform
                words = line.split(); nwords = len(words)
                if nwords == 1 :
                    if words[0].isnumeric():
                        mass = float(words[0])
                    else:
                        nam = words[0]
                else:
                    atom_names.append(nam) 
                    s = [float(x) for x in words[0:3]]
                    coords[iatm,0] = s[0]*H[0,0] + s[1]*H[1,0] + s[2]*H[2,0]
                    coords[iatm,1] = s[0]*H[0,1] + s[1]*H[1,1] + s[2]*H[2,1]
                    coords[iatm,2] = s[0]*H[0,2] + s[1]*H[1,2] + s[2]*H[2,2]
                    iatm += 1

    return H0, atom_names, coords
