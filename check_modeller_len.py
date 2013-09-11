#!/usr/bin/env python
# encoding: utf-8

import glob
import itertools
import subprocess
import os
import prody
from zam import protein
import numpy
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *


def main():
    pdb_ids = [pdb.strip() for pdb in open('proteins.txt').readlines()]

    for pdb_id in pdb_ids:
        print 'processing {}'.format(pdb_id)

        # load the system as setup
        p = protein.Protein(os.path.join('FixedModeller', pdb_id + '.pdb'))

        # get one of the rosetta models
        r = protein.Protein(glob.glob(os.path.join('final_dataset', pdb_id, 'S_*.pdb'))[0])

        print '    system:  {}'.format(len(p))
        print '    rosetta: {}'.format(len(r))

        if not len(p) == len(r):
            print '    WARNING: Lengths do not match!'
        print


if __name__ == '__main__':
    main()
