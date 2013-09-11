#!/usr/bin/env python
# encoding: utf-8


import modeller
import modeller.automodel
import os
import tempfile
import shutil
import contextlib
import glob


def perform_sequence_alignment():
    e = modeller.environ()
    m1 = modeller.model(e, file='experimental.pdb')
    m2 = modeller.model(e, file='rosetta.pdb')
    aln = modeller.alignment(e)
    aln.append_model(m1, align_codes='experimental', atom_files='experimental.pdb')
    aln.append_model(m2, align_codes='rosetta')
    aln.align2d()
    aln.write(file='align.ali', alignment_format='PIR')


def build_model():
    e = modeller.environ()
    a = modeller.automodel.automodel(e, alnfile='align.ali', knowns='experimental', sequence='rosetta')
    a.starting_model = 1
    a.ending_model = 1
    a.make()


def stage_files(pdb_id, temp_dir):
    pdb_path = os.path.join('final_dataset', pdb_id, pdb_id + '.pdb')
    rosetta_path = glob.glob(os.path.join('final_dataset', pdb_id, 'S_*'))[0]

    shutil.copy(pdb_path, os.path.join(temp_dir, 'experimental.pdb'))
    shutil.copy(rosetta_path, os.path.join(temp_dir, 'rosetta.pdb'))


def recover_model(pdb_id, temp_dir):
    model_path = glob.glob(os.path.join(temp_dir, 'rosetta.B*.pdb'))[0]
    shutil.copy(model_path, os.path.join('FixedModeller', pdb_id + '.pdb'))


def get_pdb_ids():
    with open('proteins.txt') as infile:
        return [line.strip() for line in infile]


def main():
    os.mkdir('FixedModeller')

    pdb_ids = get_pdb_ids()
    for pdb_id in pdb_ids:
        with make_temp_dir() as temp_dir:
            stage_files(pdb_id, temp_dir)
            with in_dir(temp_dir):
                perform_sequence_alignment()
                build_model()
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
