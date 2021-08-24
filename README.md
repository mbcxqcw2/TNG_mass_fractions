

This repo attempts to recreate the IllustrisTNG evolving baryonic mass fraction plots from Figure 2 in Martizzi+19 (https://arxiv.org/abs/1810.01883) and Figure A.1 from Artale+21 (https://arxiv.org/abs/2102.01092)


plot_mass_fractions_artale.ipynb loads the data in the author’s original way.


The resulting output of this notebook is evolving_mass_fraction_TNG100-3_artale.png.


plot_mass_fractions_walker.ipynb loads the data using yt.


The resulting output of this notebook is evolving_mass_fraction_TNG100-3_walker_onechunk.png.


plot_mass_fractions.ipynb also loads the data using yt. 


The resulting output of this notebook is evolving_mass_fraction_TNG100-3_walker_onechunk.png.


The difference between plot_mass_fractions.ipynb and plot_mass_fractions_walker.ipynb is that TNG100-3 is stored as 7 chunks on our computers, and the former function loads them all individually, while the latter only loads chunk 0. However, upon analysis, it appears to us that even when specifying loading a single chunk, yt actually still loads the entire simulation snapshot. When you load any individual chunk, the output looks the same as any other chunk, and you can see that the results of the two notebooks themselves are identical.


#Older codes are described below#


mass_fraction_check_1.ipynb loads the baryonic masses in halos, filaments and voids in snap_099.0.hdf5 mass_fraction_check_2.ipynb loads the baryonic masses in halos, filaments and voids in snap_099.3.hdf5 mass_fraction_check_3.ipynb loads the baryonic masses in halos, filaments and voids in snap_067.0.hdf5 mass_fraction_check_4.ipynb loads the baryonic masses in halos, filaments and voids in snap_067.3.hdf5

When comparing the results of these, the baryonic mass distribution is identical to 8 significant figures for snap_099.0.hdf5 and snap_099.3.hdf5, and identical to 8 significant figures for snap_067.0.hdf5 and snap_067.3.hdf5

plot_mass_fractions.ipynb loads the baryonic masses for all full snapshots and chunks, and calculates the corresponding baryonic mass fractions as a function of redshift.

Howeverr, upon inspection, it appears that in any given snapshot, each chunk yields identical results for the mass in its halos, filaments, and voids.

