# Test Script to import PhysNet as energy function in CHARMM via PyCHARMM

# Basics
import os
import sys
import json
import ctypes
import pandas
import subprocess
import numpy as np

# ASE basics (v 3.20.1 modified)
from ase import Atoms
from ase import io

# Optimization algorithms
from scipy.optimize import minimize, curve_fit

# Miscellaneous
from ase.visualize import view
import ase.units as units
import time

# Matplotlib
import matplotlib.pyplot as plt

# Step 0: Parameter definition
#-----------------------------------------------------------

Ha2kcalmol = units.mol*units.Hartree/units.kcal

# Step 1: Prepare data
#-----------------------------------------------------------

# Charmm run command
#charmm_command = 'charmm'

# Reference data file
data_dir = 'data'
data_file = 'data_results_SCN_{:d}.json'

# Data source files
source_dir = 'source'

# Topology data files
toppar_dir = 'source/toppar'

# Sample files directory
sample_dir = "sample_files"

# Number of sample systems
Nsystems = 4

# Initialize reference data lists
Nref = []
potential = []
residues = []
positions = []
symbols = []

# Iterate over reference samples
for ii in range(Nsystems):
    
    # Load data
    with open(os.path.join(data_dir, data_file.format(ii)), 'r') as f:
        ref_data = json.load(f)
    
    # Number of reference points
    Nref.append(np.sum([bool(ref_data[key]) for key in ref_data.keys()]))
    
    # Interaction potential list in Hartree
    potential.append(
        [item["E_int"]*Ha2kcalmol for key, item in ref_data.items()])
    
    # Residue list
    residues.append([item["residues"] for key, item in ref_data.items()])

    # Positions list
    positions.append([item["positions"] for key, item in ref_data.items()])

    # Atomic symbols list
    symbols.append([item["symbols"] for key, item in ref_data.items()])


# Step 2: Prepare van-der-Waals fit
#-----------------------------------------------------------

# Calculate potential energy
def clc_pot(*pars):
    
    # Read parameter file template
    with open(os.path.join(toppar_dir, "template_kscn.str"), 'r') as f:
        ltmp = f.read()
        
    # Write parameter
    ltmp = ltmp.replace("%epslS%", "{:.5f}".format(pars[0]))
    ltmp = ltmp.replace("%epslC%", "{:.5f}".format(pars[1]))
    ltmp = ltmp.replace("%epslN%", "{:.5f}".format(pars[2]))
    ltmp = ltmp.replace("%rmnhS%", "{:.4f}".format(pars[3]))
    ltmp = ltmp.replace("%rmnhC%", "{:.4f}".format(pars[4]))
    ltmp = ltmp.replace("%rmnhN%", "{:.4f}".format(pars[5]))
    
    # Save parameter file
    with open(os.path.join(toppar_dir, "toppar_kscn.str"), 'w') as f:
        f.write(ltmp)

    # Save parameter file
    with open("last_toppar_kscn.str", 'w') as f:
        f.write(ltmp)
    
    # Read input template
    with open(os.path.join(source_dir, "template_charmm_sample.inp"), 'r') as f:
        inplines = f.read()
    
    # Read run script template
    with open(os.path.join(source_dir, "template_charmm_sample.sh"), 'r') as f:
        runlines = f.read()
    
    # Prepare sample CHARMM inputs
    for ii in range(Nsystems):
        
        inpfile = os.path.join(sample_dir, "charmm_sample_{:d}.inp".format(ii))
        outfile = os.path.join(sample_dir, "charmm_sample_{:d}.out".format(ii))
        
        # Write sample charmm input file
        smpllines = inplines
        smpllines = smpllines.replace("%SSS%", "{:d}".format(ii))
        smpllines = smpllines.replace("%NNN%", "{:d}".format(Nref[ii]))
        with open(inpfile, 'w') as f:
            f.write(smpllines)
        
        # Add sample to run script
        runlines += "srun $my_charmm -i {:s} -o {:s}\n".format(inpfile, outfile)
        
    with open("charmm_sample.sh", 'w') as f:
        f.write(runlines)
            
    return run_pot()

# Run potential energy calculation
def run_pot():
    
    # Execute run file
    subprocess.run(["chmod", "764", "charmm_sample.sh"])
    subprocess.run(["bash",  "charmm_sample.sh"])
    
    # Interaction potential list
    V = []
    
    # Read energies
    for ii in range(Nsystems):
        
        outfile = os.path.join(sample_dir, "charmm_sample_{:d}.out".format(ii))
    
        # Read interaction energy
        Vnbond = []
        with open(outfile, 'r') as f:
            outlines = f.readlines()
        for line in outlines:
            if "ENER EXTERN>" in line:
                Vnbond.append(np.sum(np.array(line.split()[2:4], dtype=float)))
            elif "ENER DCMPOL>" in line:
                Vnbond[-1] += float(line.split()[2])
        # Add bonded and non-bonded term
        V.append(Vnbond)
        
    return V

def clc_rmse(pars, *args, **kwargs):
    
    # Get system potential
    V = clc_pot(*pars)
    
    #if len(args) == 1:
        
        #sigma_pot = args[1]
        
        ## Get potential error
        #Vrmse = []
        #for ii, poti in enumerate(potential):
            
            #Vdiff = poti - V[ii]
            #Vrmsei = np.sqrt(np.mean(sigma_pot[ii]*Vdiff**2))
            #Vrmse.append(sigma_sys[ii]*Vrmsei)
            
        #Vrmse = np.sum(Vrmse)
        
    #else:
    
    # Get potential error
    Vdiff = []
    Nkeep = 10
    Vbest = np.full((Nkeep, 5), 1000.0, dtype=float)
    Vworst = np.zeros((Nkeep, 5), dtype=float)
    for ii in range(Nsystems):
        for jj in range(Nref[ii]):
    
            #Vdiff.append(
                #potential[ii][jj] - V[ii][jj] + pars[9 + ii])
            
            #io.write(
                #os.path.join(sample_dir, "cluster_{:d}.xyz".format(ii)),
                #Atoms(symbols[ii], positions=positions[ii]),
                #comment="Fit: {:.2f}; Ref: {:.2f}; Diff: {:.2f}".format(
                    #potential[ii], V[ii], Vdiff[-1]))
                    
            Vdiff_i = potential[ii][jj] - V[ii][jj]
            Vdiff.append(Vdiff_i)
            
            # Check for best results
            for kk in range(Nkeep):
                if np.abs(Vdiff_i) < np.abs(Vbest[kk, 2]):
                    Vbest[kk] = (
                        potential[ii][jj], V[ii][jj], Vdiff_i, 
                        ii, jj)
                    break
            # Check for worst results
            for kk in range(Nkeep):
                if np.abs(Vdiff_i) > np.abs(Vworst[kk, 2]):
                    Vworst[kk] = (
                        potential[ii][jj], V[ii][jj], Vdiff_i, 
                        ii, jj)
                    break
            
    # Returm RMSE
    Vrmse = np.sqrt(np.mean(np.array(Vdiff)**2))
    
    with open('step_overview.txt', 'a') as f:
        f.write("Minimize: Current pars: " + str(pars) + "\n")
        f.write("Minimize: Current RMSE: " + str(Vrmse) + "\n")
        status = (
            "Best: Vref    Vfit    Vdiff    System/Index    " + 
            "Worst: Vref    Vfit    Vdiff    System/Index\n")
        for kk in range(Nkeep):
            status += (
                "  {: 7.2f}  {: 7.2f}  {: 7.2f}  {:2d}/{:2d}    ".format(
                    Vbest[kk, 0], Vbest[kk, 1], Vbest[kk, 2], 
                    int(Vbest[kk, 3]), int(Vbest[kk, 4]))
                + "  {: 7.2f}  {: 7.2f}  {: 7.2f}  {:2d}/{:2d}\n".format(
                    Vworst[kk, 0], Vworst[kk, 1], Vworst[kk, 2],  
                    int(Vworst[kk, 3]), int(Vworst[kk, 4]))
                )
        f.write(status)
    return Vrmse

# V0
epsilon = np.array([-0.3639, -0.1016, -0.0741], dtype=float)
rminhlf = np.array([1.9755, 1.8801, 1.8577], dtype=float)

#V1
epsilon = np.array([-1.83593909e-01, -1.0e-04, -2.22628666e-02], dtype=float)
rminhlf = np.array([ 2.42794766, 1.89344831, 2.28604331], dtype=float)

# Pack to parameter array
pars = np.concatenate((epsilon, rminhlf))

# Step 4: Start van-der-Waals fit
#-----------------------------------------------------------

#clc_rmse(pars)
#exit()

if True:
    
    # Start optimization
    bounds = np.array([
        (-1.0, -0.0001),
        (-1.0, -0.0001),
        (-1.0, -0.0001),
        (1.5, 3.0),
        (1.5, 3.0),
        (1.5, 3.0),
        ])

    res = minimize(
        clc_rmse, pars, 
        method='TNC',
        options={'eps': 1e-4},
        bounds=bounds)
    res = res.x

    pars = res
    epsilon = np.array([pars[0], pars[1], pars[2]])
    rminhlf = np.array([pars[3], pars[4], pars[5]])

    # Get interaction potential
    all_Vnbond = clc_pot(*pars)
    
all_Vnbond = run_pot()
print(len(all_Vnbond), [len(Vi) for Vi in all_Vnbond])
print(len(potential), [len(Vi) for Vi in potential])

# Save results
np.save(
    "ref_potential", 
    np.array(potential))
np.save(
    "fit_potential", 
    np.array(all_Vnbond))


# Step 5: Plot nonbonding potential
#-----------------------------------------------------------

# Plot properties

# Fontsize
SMALL_SIZE = 12
MEDIUM_SIZE = 14
BIGGER_SIZE = 16

#plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rc('font', size=SMALL_SIZE, weight='bold')  # controls default text sizes
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# Graphical output format type
gout_type = 'png'
dpi = 200

# Figure size
figsize = (6, 6)
sfig = float(figsize[0])/float(figsize[1])

# Figure
fig = plt.figure(figsize=figsize)

# Axes arrangement
left = 0.15
bottom = 0.15
column = [0.75, 0.10]
row = column[0]*sfig

# Axes
axs1 = fig.add_axes(
    [left, bottom, column[0], row])
    

# Potential range
flat_Vnbond = [Vi for Vlist in all_Vnbond for Vi in Vlist]
flat_potential = [Vi for Vlist in potential for Vi in Vlist]
emin = np.min([np.min(flat_Vnbond), np.min(flat_potential)])
emax = np.max([np.max(flat_Vnbond), np.max(flat_potential)])
de = emax - emin

all_rmse = np.sqrt(
    np.mean((np.array(flat_Vnbond) - np.array(flat_potential))**2))
color = ['blue', 'red', 'green', 'magenta']

# Iterate over systems
for isys in range(Nsystems):

    # Data
    sys_Vnbond = all_Vnbond[isys]
    sys_potential = potential[isys]
    
    # Range
    axs1.set_xlim(emin - 0.1*de, emax + 0.1*de)
    axs1.set_ylim(emin - 0.1*de, emax + 0.1*de)
    if isys == 0:
        axs1.plot(
            [emin - 0.1*de, emax + 0.1*de], 
            [emin - 0.1*de, emax + 0.1*de], '-k')
    
    
    # Plot
    sys_rmse = np.sqrt(
        np.mean((np.array(sys_Vnbond) - np.array(sys_potential))**2))
    label = "Sys{:d}, RMSE = {:.2f} kcal/mol".format(isys, sys_rmse)
    axs1.plot(
        sys_potential, sys_Vnbond, 'o', color=color[isys], 
        mfc='None', label=label)
    
    axs1.set_xlabel(
        r'$\Delta$E$_\mathrm{M06-2X}$ (kcal/mol)', fontweight='bold')
    axs1.get_xaxis().set_label_coords(0.5, -0.12)
    axs1.set_ylabel(r'$\Delta$E$_\mathrm{FF}$ (kcal/mol)', fontweight='bold')
    axs1.get_yaxis().set_label_coords(-0.12, 0.5)
    
    axs1.set_xticks(np.arange(int(np.floor(emin)), int(np.ceil(emax)), 10))
    axs1.set_yticks(np.arange(int(np.floor(emin)), int(np.ceil(emax)), 10))
    
    axs1.legend(loc='lower right')
    
plt.savefig(
    "fit_Ecorr.{:s}".format(gout_type),
    format=gout_type, dpi=dpi)
plt.close()










