#!/usr/bin/env python

def write(system, fn, title='', unit='nm'):
    """
    Writes atom coordinates to an XYZ file.  
    
    Parameters
    ---------
    file_unit : 'nm' | 'angstrom'

    """
    if file_unit == 'nm':
        factor = 1
    elif file_unit == 'angstrom': 
        factor = 10
    else:
        raise ValueError()

    na = system.num_atoms
    with open(fn,'w') as fh:
        fh.write(str(na) + '\n')
        fh.write(title + '\n')
    
        for i in range(1, na+1):
            atm_nam = system.atoms[i].name
            coords = factor*system.atoms[i].coords
            fh.write( '%s  '%atm_nam 
                + '  '.join(['% .15g'%x for x in coords]) + '\n')



def read(system, fn, offset=0, aids=None, file_unit='nm'):
        """
        Reads atom positions from an XYZ file.
    
        """
        if file_unit == 'nm':
            factor = 1
        elif file_unit == 'angstrom': 
            factor = 0.1
        else:
            raise ValueError()        

        with open(fn, 'r') as fh:
            lines = fh.readlines()
    
        na = int(lines[0].strip('\n'))
        assert na <= system.num_atoms
        if aids:
            aids_ = aids
        else:
            aids_ = [offset+i for i in range(1, na+1)]
        for i in range(2, 2+na):
            words = lines[i].strip('\n').split()
            coords = np.array([float(x) for x in words[1:4]])
            iatm = aids_[i-2]
            system.atoms[iatm].coords[:] = factor*coords
