"""
Python implementation of a ring current correction to backbone amide
NMR shieldings based on work by Christensen et el (10.1021/ct2002607)

Copyright (c) 2016, Casper Steinmann <casper.steinmann@gmail.com>
"""

import sys

import numpy
import openbabel

# incomplete implementation of the OBAminoAcidProperty enum (as a dictionary)
# http://openbabel.org/dev-api/namespaceOpenBabel_1_1OBAminoAcidProperty.shtml
AminoAcidProperties = {'AROMATIC': 3}

def RingCenterNormal(residue, ATOMS):
    """ Calculates vector for center of ring for aromatic sidechains
        and normal to the ring plane for a residue.

        Uses the variable ATOMS to select ring atoms

        Arguments:
        ----------
        residue : OBResidue in a molecule
        ATOMS : list of atom names to filter through

        Returns:
        --------
        tuple with vector to center of ring and normal vector
    """
    R = numpy.zeros((len(ATOMS),3))

    i_atom = 0
    for atom in openbabel.OBResidueAtomIter(residue):
        name = residue.GetAtomID(atom).strip()
        if name in ATOMS:
            R[i_atom][0] = atom.GetX()
            R[i_atom][1] = atom.GetY()
            R[i_atom][2] = atom.GetZ()
            i_atom += 1

    # calculate center
    center = R.mean(axis=0)
    
    # generate normal
    v1 = R[0]
    v2 = R[3]
    n = numpy.cross(v1, v2)
    n /= numpy.linalg.norm(n)

    return center, n

def HistidoneCenterNormal(residue):
    """ Neutral histidone """
    raise NotImplementedError("HIS has not been implemented")

def HistidoneCenterNormalp(residue):
    """ Protonated histidone """
    raise NotImplementedError("HIS+ has not been implemented")

def PhenylAlanineCenterNormal(residue):
    """ Phenylalanine """
    PHE_ATOMS = ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]
    return RingCenterNormal(residue, PHE_ATOMS)

def Tryptophan5CenterNormal(residue):
    """ 5 member ring of Tyrosine """
    TRP_ATOMS = ["CG", "CD1", "CD2", "NE1", "CE2"]
    return RingCenterNormal(residue, TRP_ATOMS)

def Tryptophan6CenterNormal(residue):
    """ 6 member ring of Tyrosine """
    TRP_ATOMS = ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"]
    return RingCenterNormal(residue, TRP_ATOMS)

def TyrosineCenterNormal(residue):
    """ Tyrosine

        INFO: Tyrosine shares atom names with Phenylalanine
    """
    return PhenylAlanineCenterNormal(residue)


# Data from the paper (10.1021/ct2002607) for the point-dipole model of HN
RCPD_FN = {'PHE':  PhenylAlanineCenterNormal,
           'TYR':  TyrosineCenterNormal,
           'HIS+': HistidoneCenterNormalp,
           'HIS':  HistidoneCenterNormal,
           'TRP5': Tryptophan5CenterNormal,
           'TRP6': Tryptophan6CenterNormal}


# intensity factors from the paper (10.1021/ct2002607)
RCPD_HN = {'PHE': 1.0, 'TYR': 0.81, 'HIS+': 0.69, 'HIS': 0.68, 'TRP5': 0.57, 'TRP6': 1.02}

# overall scaling factor from the paper (10.1021/ct2002607)
RCPD_B = 30.42 # ppm Angstrom^3

class Ring(object):
    """ Simple representation of a ring.

        This structure has all the data needed to later evaluate the
        ring current effect from a given ring onto another HN atom.
    """

    def __init__(self, r, n, HN, B, name):
        self.r = numpy.array(r)
        self.n = numpy.array(n)
        self.HN = HN
        self.B = B
        self.name = name

    def __repr__(self):
        r = "[{0:9.4f}, {1:9.4f}, {2:9.4f}]".format(self.r[0], self.r[1], self.r[2])
        n = "[{0:9.4f}, {1:9.4f}, {2:9.4f}]".format(self.n[0], self.n[1], self.n[2])
        return "Ring({}, {}, {}, {}, '{}')".format(r, n, self.HN, self.B, self.name)

    def __str__(self):
        """ Print-friendly represenation """
        r = "({0:9.4f}, {1:9.4f}, {2:9.4f})".format(self.r[0], self.r[1], self.r[2])
        s = "{0:s} at R = {1:s}".format(self.name, r)
        return s

def load(filename):
    """ Loads a .pdb file into an openbabel molecule """
    obc = openbabel.OBConversion()
    obc.SetInFormat("pdb")

    mol = openbabel.OBMol()
    obc.ReadFile(mol, filename)
    return mol

def Rings(mol):
    """ Returns a list of ring data from mol.
 
        All data is based on the point-dipole model by Christensen et al (DOI: 10.1021/ct2002607)

        Arguments:
        ----------
        mol : OBMol molecule with molecule information (loaded from .pdb file)

        Returns:
        --------
        list of rings (instatiated from the Ring class)
    """

    rings = []
    B = RCPD_B
    for i_residue, residue in enumerate(openbabel.OBResidueIter(mol)):
        name = residue.GetName()
        if not name in ["SOL", "NA"]:
            aromatic_key = AminoAcidProperties['AROMATIC']
            is_aromatic = residue.GetAminoAcidProperty(aromatic_key)
            if is_aromatic:
                if name == 'PHE':
                    r, n = RCPD_FN[name](residue)
                    B = RCPD_B
                    HN = RCPD_HN[name]
                    rings.append( Ring(r, n, HN, B, name) )

                if name == 'TYR':
                    r, n = RCPD_FN[name](residue)
                    B = RCPD_B
                    HN = RCPD_HN[name]
                    rings.append( Ring(r, n, HN, B, name) )

                # Tryptophan has both a 5- and a 6-member ring
                if name == 'TRP':
                    r, n = RCPD_FN['TRP5'](residue)
                    B = RCPD_B
                    HN = RCPD_HN['TRP5']
                    rings.append( Ring(r, n, HN, B, name) )

                    r, n = RCPD_FN['TRP6'](residue)
                    HN = RCPD_HN['TRP6']
                    rings.append( Ring(r, n, HN, B, name) )
    return rings

def AtomLists(mol, atom_names):
    """ Returns indices and vectors of all HN atoms in an openbabel molecule

        Arguments:
        ----------
        mol : openbabel molecule

        Returns:
        --------
        tuple of indices and vectors to each HN atoms
    """
    atoms = []
    idx = []
    for i_residue, residue in enumerate(openbabel.OBResidueIter(mol)):
        name = residue.GetName()

        if not name in ["SOL", "NA"]:
            for atom in openbabel.OBResidueAtomIter(residue):
                atom_name = residue.GetAtomID(atom).strip()
                if atom_name in atom_names:
                    r = [atom.GetX(), atom.GetY(), atom.GetZ()]
                    atoms.append(numpy.array(r))
                    idx.append(atom.GetIdx())

    return idx, atoms

def RingCurrentCorrection(rings, h_r):
    """ Calculates the ring current correction based on rings and atom coordinates

        NOTE: This does not distinguish between atom types, i.e. HN and HA so be careful
              as this work is only done for HN.

              Results for HA has been attemped using the same data and errors of around
              5% - 10% is to be expected (private communication with first author of paper)
    """
    data = numpy.zeros((len(h_r), len(rings)))
    for i_ring, ring in enumerate(rings):
       for ii, r in enumerate(h_r):
           dr = r - ring.r
           tang = dr.dot(ring.n) / numpy.linalg.norm(dr)

           R2 = dr.dot(dr)
           R = numpy.sqrt(R2)
           G = (1.0 - 3.0*tang*tang) / (R2*R)  # geometric factor for point-dipole model (eq 1 in SI of 10.1021/ct2002607)
           rc = ring.HN * ring.B * G           # ring current correction (eq 8 in 10.1021/ct2002607)
           data[ii][i_ring] = rc

    m = data.mean(axis=1)
    return m

if __name__ == '__main__':

    if len(sys.argv) == 1:
        exit('usage: rc.py pdbfile1 [pdbfile2 ...]')

    l = []
    for file in sys.argv[1:]:
        mol = load(file)
        rings = Rings(mol)
        if len(rings) == 0:
            continue

        print "Found {0:3d} rings.".format(len(rings))
        for iring, ring in enumerate(rings):
            print "{0:5d} - {1:s}".format(iring+1, ring)

        hn_idx, hn_r = AtomLists(mol, ["HN", "HA"]) # search through HN and HA atoms
        data = RingCurrentCorrection(rings, hn_r)
        l.append(data)

    # output averaged ring current shift for each atom that was searched
    dd = numpy.array(l)
    mm = numpy.mean(dd, axis=0)
    
    print ""
    print "Final ring-current corrections:".format()
    for id, m in zip(hn_idx, mm):
        print "{0:4d}{1:9.3f}".format(id, m)
