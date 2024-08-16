
#---------
# Modules
#---------

module load charmm/c45a1-gcc9.2.0-ompi4.0.2

my_charmm="/data/toepfer/Project_Eutectic/FDCM_calc/dev-release-dcm/build/cmake/charmm"
ulimit -s 10420



srun $my_charmm -i sample_files/charmm_sample_0.inp -o sample_files/charmm_sample_0.out
srun $my_charmm -i sample_files/charmm_sample_1.inp -o sample_files/charmm_sample_1.out
srun $my_charmm -i sample_files/charmm_sample_2.inp -o sample_files/charmm_sample_2.out
srun $my_charmm -i sample_files/charmm_sample_3.inp -o sample_files/charmm_sample_3.out
