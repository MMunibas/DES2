# Test Script to import PhysNet as energy function in CHARMM via PyCHARMM

# Basics
import os
import numpy as np

# Matplotlib
import matplotlib.pyplot as plt

# Miscellaneous
import ase.units as units

# Number of sample systems
systems_num = 4

# System label
systems_label = [
    r"16 H$_2$O",
    r"14 H$_2$O + 1 SCN$^-$",
    r"14 H$_2$O + 1 K$^+$",
    r"12 H$_2$O + 1 SCN$^-$ + 1 K$^+$"
    ]


# Source files
systems_ref = "ref_potential.npy"
systems_fit = "fit_potential.npy"

# Plot properties

# Fontsize
SMALL_SIZE = 10
MEDIUM_SIZE = 12
BIGGER_SIZE = 14

#plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rc('font', size=MEDIUM_SIZE, weight='bold')  # controls default text sizes
plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# Graphical output format type
gout_type = 'png'
dpi = 200

# Load data
data_ref = np.load(systems_ref)
data_fit = np.load(systems_fit)


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
emin = np.min([np.min(data_ref), np.min(data_fit)])
emax = np.max([np.max(data_ref), np.max(data_fit)])
de = emax - emin

# Compute RMSE
systems_rmse = [
    np.sqrt(np.mean((np.array(data_i) - np.array(data_j))**2))
    for (data_i, data_j) in zip(data_ref, data_fit)]
systems_rmse_all = np.sqrt(
    np.mean(
        (np.array(data_ref.reshape(-1)) - np.array(data_fit.reshape(-1)))**2
        )
    )

# Color and Marker code
systems_color = ['blue', 'red', 'green', 'magenta']
systems_marker = ['o', 's', 'd', '^']

# Plot diagonal
axs1.plot(
    [emin - 0.1*de, emax + 0.1*de], 
    [emin - 0.1*de, emax + 0.1*de], '-k')

# Iterate over systems
for isys in range(systems_num):
    
    # Plot data
    label = (
        systems_label[isys] 
        + "\nRMSE = {:.2f} kcal/mol".format(systems_rmse[isys]))
    axs1.plot(
        data_ref[isys], data_fit[isys],
        color=systems_color[isys], 
        marker=systems_marker[isys],
        mfc='None', 
        ls='None',
        label=label)

# Range
axs1.set_xlim(emin - 0.1*de, emax + 0.1*de)
axs1.set_ylim(emin - 0.1*de, emax + 0.1*de)

axs1.set_xlabel(
    r'$\Delta E_\mathrm{M06-2X}$ (kcal/mol)', fontweight='bold')
axs1.get_xaxis().set_label_coords(0.5, -0.12)
axs1.set_ylabel(r'$\Delta E_\mathrm{FF}$ (kcal/mol)', fontweight='bold')
axs1.get_yaxis().set_label_coords(-0.12, 0.5)

axs1.set_xticks(np.arange(
    int(np.floor(emin/10.)*10.), int(np.ceil(emax/10.)*10.), 25))
axs1.set_yticks(np.arange(
    int(np.floor(emin/10.)*10.), int(np.ceil(emax/10.)*10.), 25))

axs1.legend(loc='upper left')

plt.savefig(
    "fit_Ecorr.{:s}".format(gout_type),
    format=gout_type, dpi=dpi)
plt.close()
