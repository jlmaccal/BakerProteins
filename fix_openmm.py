#!/usr/bin/env python
# encoding: utf-8


import os
import itertools
import prody
import tempfile
import shutil
import contextlib
import glob
import subprocess
from zam import protein
import numpy as np


def fix_openmm():
    # get the whole crystal structure
    # get only the ATOM records
    # and HETAM records for MSE
    # convert MSE to MET
    with open('no_smet.pdb', 'w') as outfile:
        with open('experimental.pdb') as infile:
            for line in infile:
                if line.startswith('ATOM'):
                    outfile.write(line)
                if line.startswith('HETATM'):
                    if line[17:20] == 'MSE':
                        atom_name = line[12:17]
                        if atom_name == 'SE   ':
                            atom_name = ' SD  '
                        line_fixed = 'ATOM  ' + line[6:12] + atom_name + 'MET' + line[20:67] + '\n'
                        outfile.write(line_fixed)

    # load the file into prody
    p = prody.parsePDB('no_smet.pdb')
    p = p.select('not hydrogen')

    # get one of the rosetta models
    r = prody.parsePDB('rosetta.pdb')

    # perform an alignment to find out what part of the crystal structure
    # corresponds to the rosetta file
    match = prody.matchChains(r, p, subset='all', overlap=25, pwalign=True)[0][1]
    print len(match)
    prody.writePDB('chain.pdb', match)

    # now clean it up with pdb fixer
    subprocess.check_call('python ~/Source/PdbFixer/pdbfixer.py chain.pdb', shell=True)

    # now load it with zam
    p = protein.Protein('output.pdb')
    p.Dehydrogen()
    disulfide_pairs = find_disulfide(p)
    for r1, r2 in disulfide_pairs:
        print '    added disulfide between {} and {}'.format(r1, r2)
        p.Res[r1].FullName = 'CYX'
        p.Res[r2].FullName = 'CYX'
    p.WritePdb('start.pdb')

    # now run tleap
    print '    running tleap'
    run_tleap(disulfide_pairs)


tleap_string = '''
set default PBradii mbondi3
source leaprc.ff12SB
sys = loadPdb start.pdb
{disulfide_string}
check sys
saveAmberParm sys system.top system.mdcrd
savePdb sys system.pdb
quit
'''


def run_tleap(disulfide_pairs):
    disulfide_string = gen_disulfide_string(disulfide_pairs)
    open('tleap.in', 'w').write(tleap_string.format(disulfide_string=disulfide_string))
    subprocess.check_call('tleap -f tleap.in > tleap.out', shell=True)


def gen_disulfide_string(disulfide_pairs):
    strings = []
    for r1, r2 in disulfide_pairs:
        strings.append('bond sys.{}.SG sys.{}.SG'.format(r1 + 1, r2 + 1))  # +1 for 0 vs 1-based indexing
    return ''.join(strings)


def find_disulfide(p):
    seq = p.Seq
    disulfides = []
    cysteines = [i for i, res in enumerate(seq) if res == 'CYS']
    for res1, res2 in itertools.combinations(cysteines, 2):
            index_1 = p.AtomInd(ResNum=res1, AtomName='SG')[0]
            index_2 = p.AtomInd(ResNum=res2, AtomName='SG')[0]
            dist = np.linalg.norm(p.Pos[index_1, :] - p.Pos[index_2, :])
            if dist < 2.10:
                # we have a disulfide bond
                disulfides.append((res1, res2))
    return disulfides


def stage_files(pdb_id, temp_dir):
    pdb_path = os.path.join('final_dataset', pdb_id, pdb_id + '.pdb')
    rosetta_path = glob.glob(os.path.join('final_dataset', pdb_id, 'S_*'))[0]

    shutil.copy(pdb_path, os.path.join(temp_dir, 'experimental.pdb'))
    shutil.copy(rosetta_path, os.path.join(temp_dir, 'rosetta.pdb'))


def recover_model(pdb_id, temp_dir):
    model_path = os.path.join(temp_dir, 'system.pdb')
    shutil.copy(model_path, os.path.join('FixedOpenMM', pdb_id + '.pdb'))


def get_pdb_ids():
    with open('proteins.txt') as infile:
        return [line.strip() for line in infile]


def main():
    os.mkdir('FixedOpenMM')

    pdb_ids = get_pdb_ids()
    for pdb_id in pdb_ids:
        with make_temp_dir() as temp_dir:
            stage_files(pdb_id, temp_dir)
            with in_dir(temp_dir):
                fix_openmm()
            recover_model(pdb_id, temp_dir)


@contextlib.contextmanager
def make_temp_dir():
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    shutil.rmtree(temp_dir)


@contextlib.contextmanager
def in_dir(directory):
    pwd = os.getcwd()
    os.chdir(directory)
    yield
    os.chdir(pwd)


if __name__ == '__main__':
    main()
