# GrainRecon
GrainRecon is a code package for reconstructing 3D grain map, namely grain orientations, positions and shapes, from diffraction patterns of laboratory X-ray diffraction contrast tomography (LabDCT). The code was developed by Haixing Fang while he work at CNRS SIMaP laboratory and ESRF in Grenoble, France, together with Dr. Pierre Lhuissier and Prof. dr. Wolfgang Ludwig. This work is part of the project Advanced Laboratory X-ray Microtomography funded by the Agence Nationale de la Recherche (ANR-18-CE42-0005). The code may be continuously updated as work progresses.

# Brief introduction of LabDCT
LabDCT is an emerging technique for 3D grain mapping using lab X-ray source, inspired by the principle of synchrotron DCT. In LabDCT, a conical white beam illuminates a sample and diffracted beams from grains fulfilling the Bragg's condition are recorded by a 2D detector behind the sample. Each diffraction spot contains a small variation of X-ray wavelenths as incident angles are slightly different given the same {hkl} of the same grain because of the cone beam geometry. In addition, spots projected onto the detector have different magnifications in directions perpendicular and parallel to the diffraction vector. These make LabDCT rather different from synchrotron DCT.

# Four different methods have been developed
1) 'cpu' - run on cpu cores only and does not require GPU, slow. It may take ~7 days to complete a volume with 6 million voxels;
2) 'gpu_cuda_Gt' - run on NVIDIA gpu and require CUDA driver, fast and robust. More suitable for a roughly know geometry;
3) 'gpu_cuda_comp' - run on NVIDIA gpu and require CUDA driver, fast and robust. Suitable for projections with serious spot overlapping;
4) 'gpu_cuda_index_compete' - run on NVIDIA gpu and require CUDA driver, based on brute force indexing and competing completeness for assign orientation for un-indexed voxels. However, this method has not be optimized to be robust yet.

# How to use the code?
1) Prepare for recontruction:
  - run get_spots.m for spot segmentation from LabDCT projections;
  - run get_tomo_slices_create_h5.m to obtain an hdf5 file for tomo volume mask;
  - run get_geometry.m to get geometry information and/or set geometry parameters manually in setup_exp.m;
  - set all parameters, including sample name, output folder, recontruction parameter etc. in setup_para.m.
2) Run grain mapping:
  - for a fresh reconstruction, execute the function run_GrainMapping_fullvol_fun.m by assinging a correct SampleFlag.
  - for an interrupted reconstruction, execute the function run_GrainMapping_fullvol_fun_continue.m to continue where it has been interrupted.
  - in the output folder, files in the format of e.g. 1DS.mat, pos1.mat, SeedingIndexing1.mat will be stored and the number corresponds to the iteration number.
3) Merge regions and identify grains:
  - run GrainMapping_assemble_vol.m, it will generate an .h5 file, a .dream3d file and a .xdmf file, which you can use ParaView for 3D visualization
  - Optionally, you can force this function to use the strategy of "completeness competing" to correct suspicious voxels and compute for empty voxels. By default, grains with fewer than 5 voxels or with seeding completeness lower than a threshold will be checked by this strategy.
4) Geometry fitting:
  - If geometry is needed to be fitting (this is usually the case with the first time grain reconstruction), run GrainMapping_geo_fit.m to fit 7 parameters, including P0y, P0z, dety00, detz00, tilt_x, tilt_y and tilt_z.
  - After geometry fitting, you may repeat steps of 2-3 to redo/refine grain mapping results.

All the codes have been tested executable with a Matlab version 2020a or later.

# Example
An example of a virtual Fe sample containing 6 grains have been included. The reconstruction time is about 20 min with gpu methods.


# License
This package is free to use, ditribute and adapt for non-commercial use only. See LICENSE for license rights and limitations (CC BY-NC 4.0).

# Reference
A manuscript has been prepared and now is under review:
H. Fang, W. Ludwig, P. Lhuissier, Reconstruction algorithms for grain mapping by laboratory X-ray diffraction contrast tomography (in review).
Please cite this article if you use or get inspired by the code presented here.

## Contact emails: haixing.fang@grenoble-inp.fr, haixing0a@esrf.fr or haixingfang868@gmail.com
