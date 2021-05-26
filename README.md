This repo iss an attempt to recreate the IllustrisTNG evolving baryonic mass fraction plots from Figure 2 in Martizzi+19 (https://arxiv.org/abs/1810.01883) and Figure A.1 from Artale+21 (https://arxiv.org/abs/2102.01092)

mass_fraction_check_1.ipynb loads the baryonic masses in halos, filaments and voids in snap_099.0.hdf5
mass_fraction_check_2.ipynb loads the baryonic masses in halos, filaments and voids in snap_099.3.hdf5
mass_fraction_check_3.ipynb loads the baryonic masses in halos, filaments and voids in snap_067.0.hdf5
mass_fraction_check_4.ipynb loads the baryonic masses in halos, filaments and voids in snap_067.3.hdf5

When comparing the results of these, the baryonic mass distribution is identical to 8 significant figures for snap_099.0.hdf5 and snap_099.3.hdf5, and identical to 8 significant figures for snap_067.0.hdf5 and snap_067.3.hdf5

plot_mass_fractions.ipynb loads the baryonic masses for all full snapshots and chunks, and calculates the corresponding baryonic mass fractions as a function of redshift.

Howeverr, upon inspection, it appears that in any given snapshot, each chunk yields identical results for the mass in its halos, filaments, and voids.
