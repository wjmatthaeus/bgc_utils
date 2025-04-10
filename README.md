# bgc_utils
bgc_utils is a collection of utilities for the Biome-BGC family of models. This software is in development and comes with no guarantees.

**Input file driven experiments**: bgcLauncher takes directory names containing MET and EPC files creates directories for all model runs, calls iniWriter to create ini files for every combination of MET and EPC files, runs individual models.

**Direct parameter variation driven experiments**: in development

**MET File creation from GCM climate data**: nc_to_mtc43.sh is a set of shell commands to convert GCM climate data in .nc files to .mtc43 for *BGC input 

**Paleo-BGC+ analysis R Script**: companion to data repository 10.5281/zenodo.6588761 and manuscript "Stems Matter: Xylem Physiological Limits Are an Accessible and Critical Improvement to Models of Plant Gas Exchange in Deep Time"

**binary_grid_..R** processes a directory of binary outputs produced by running x_BGC using x_Launcher.sh

**EEDT_cesm..R** reproduces the analysis for the manuscript Matthaeus et al. "Vegetation turnovers reduced water availability during the last icehouse" using the data hosted at osf.io
