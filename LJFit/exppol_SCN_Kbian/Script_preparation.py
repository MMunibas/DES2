# Basics
import os
import json
import ctypes
import pandas
import numpy as np

# PyCHARMM
import pycharmm
import pycharmm.generate as gen
import pycharmm.ic as ic
import pycharmm.coor as coor
import pycharmm.energy as energy
import pycharmm.dynamics as dyn
import pycharmm.nbonds as nbonds
import pycharmm.minimize as minimize
import pycharmm.crystal as crystal
import pycharmm.image as image
import pycharmm.psf as psf
import pycharmm.read as read
import pycharmm.write as write
import pycharmm.settings as settings
import pycharmm.lingo as stream
import pycharmm.select as select
import pycharmm.shake as shake
import pycharmm.cons_fix as cons_fix
import pycharmm.cons_harm as cons_harm
from pycharmm.lib import charmm as libcharmm


# Sample directory
sampledir = "sample_files"
if not os.path.exists(sampledir):
    os.makedirs(sampledir)

# Reference data file
data_dir = "data"
data_file = "data_results_SCN_{:d}.json"
Ncomp = 4

# Load topology files
stream.charmm_script("stream source/toppar.str")

for ii in range(Ncomp):

    # Load data
    with open(os.path.join(data_dir, data_file.format(ii)), 'r') as f:
        ref_data = json.load(f)
    
    # Residue list
    residues = [item["residues"] for key, item in ref_data.items()]

    # Positions list
    positions = [item["positions"] for key, item in ref_data.items()]

    # Generate system
    residue_string = ("{:s} "*len(residues[0])).format(*residues[0])
    read.sequence_string(residue_string)
    gen.new_segment(seg_name='SYS')

    # Write psf file
    write.psf_card(
        "{:s}/sample_{:d}.psf".format(sampledir, ii), 
        title="SCN water cluster sample")

    for jj in np.array(list(ref_data.keys()), dtype=int):
        
        pos = np.array(positions[jj])
        pandas_pos = pandas.DataFrame({
            'x': pos[:, 0], 'y': pos[:, 1], 'z': pos[:, 2]})
        coor.set_positions(pandas_pos)
        
        # Write pdb and psf files
        write.coor_card(
            "{:s}/sample_{:d}_{:d}.crd".format(sampledir, ii, jj), 
            title="SCN water cluster sample")
        
    stream.charmm_script(
        "delete atom sele all end")
    
