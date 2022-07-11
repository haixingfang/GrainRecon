# GrainRecon
GrainRecon [(initial site at Univ. Grenoble Alpes)](https://gricad-gitlab.univ-grenoble-alpes.fr/TomoX_SIMaP/GrainRecon) is a code package for reconstructing 3D grain map, namely grain orientations, positions and shapes, from diffraction patterns of laboratory X-ray diffraction contrast tomography (LabDCT). The code was developed by [Dr. Haixing Fang](https://orcid.org/0000-0001-8114-5276) while he works at [CNRS SIMaP laboratory](https://simap.grenoble-inp.fr/) and [ESRF](https://www.esrf.fr/UsersAndScience/Experiments/StructMaterials/ID11) in Grenoble, France, together with [Dr. Pierre Lhuissier](https://simap.grenoble-inp.fr/fr/equipes/m-lhuissier-pierre) at SIMaP and [Prof. dr. Wolfgang Ludwig](https://scholar.google.fr/citations?user=f8-PwEMAAAAJ&hl=fr) at [MatéIS](https://mateis.insa-lyon.fr/fr/content/ludwig-wolfgang) and [ESRF](https://www.esrf.fr/UsersAndScience/Experiments/StructMaterials/ID11). This work is part of the project Advanced Laboratory X-ray Microtomography funded by the Agence Nationale de la Recherche (ANR-18-CE42-0005). The code may be continuously updated as work progresses.

# Brief introduction of LabDCT
LabDCT is a novel technique for 3D grain mapping using lab X-ray source, inspired by the principle of synchrotron DCT. In LabDCT, a conical white beam illuminates a sample and diffracted beams from grains fulfilling the Bragg's condition are recorded by a 2D detector behind the sample. Each diffraction spot contains a small variation of X-ray wavelenths as incident angles are slightly different given the same {hkl} of the same grain because of the cone beam geometry. In addition, spots projected onto the detector have different magnifications in directions perpendicular and parallel to the diffraction vector. These make LabDCT rather different from synchrotron DCT.

# 3(4) different methods have been developed
1) 'cpu' - run on cpu cores only and does not require GPU, slow. It may take ~7 days to complete a volume with 6 million voxels;
2) 'gpu_cuda_Gt' - run on NVIDIA gpu and require CUDA driver, fast and robust. More suitable for a roughly know geometry;
3) 'gpu_cuda_comp' - run on NVIDIA gpu and require CUDA driver, fast and robust. Suitable for projections with serious spot overlapping;
4) 'gpu_cuda_index_compete' - run on NVIDIA gpu and require CUDA driver, based on brute force indexing and competing completeness for assign orientation for un-indexed voxels. However, this method has not be optimized to be robust enough yet and thus so far not recommended to use.

# How to use the code? Generally, 4 steps:
## 1) Prepare for recontruction
  - run [get_spots.m](https://github.com/haixingfang/GrainRecon/blob/main/get_spots.m) for spot segmentation from LabDCT projections;
  - run [get_tomo_slices_create_h5.m](https://github.com/haixingfang/GrainRecon/blob/main/get_tomo_slices_create_h5.m) to obtain an hdf5 file for tomo volume mask;
  - run [get_geometry.m](https://github.com/haixingfang/GrainRecon/blob/main/get_geometry.m) to get geometry information and/or set geometry parameters manually in setup_exp.m;
  - set all parameters, including sample name, output folder, recontruction parameter etc. in [setup_para.m](https://github.com/haixingfang/GrainRecon/blob/main/setup_para.m).
## 2) Run grain mapping
  - for a fresh reconstruction, execute the function [run_GrainMapping_fullvol_fun.m](https://github.com/haixingfang/GrainRecon/blob/main/run_GrainMapping_fullvol_fun.m) by assinging a correct SampleFlag.
  - for an interrupted reconstruction, execute the function [run_GrainMapping_fullvol_fun_continue.m](https://github.com/haixingfang/GrainRecon/blob/main/run_GrainMapping_fullvol_fun_continue.m) to continue where it has been interrupted.
  - in the output folder, files in the format of e.g. 1DS.mat, pos1.mat, SeedingIndexing1.mat will be stored and the number corresponds to the iteration number. A file of para.mat records all the parameters.
## 3) Merge regions and identify grains
  - run [GrainMapping_assemble_vol.m](https://github.com/haixingfang/GrainRecon/blob/main/GrainMapping_assemble_vol.m), it will generate a .h5 file, a .dream3d file and a .xdmf file, which you can use ParaView for 3D visualization
  - Optionally, you can force this function to use the strategy of "completeness competing" to correct suspicious voxels and compute for empty voxels. By default, grains with fewer than 5 voxels or with seeding completeness lower than a threshold will be checked by this strategy.
## 4) Geometry fitting
  - If geometry is needed to be fitting (this is usually the case with the first time grain reconstruction), run [GrainMapping_geo_fit.m](https://github.com/haixingfang/GrainRecon/blob/main/GrainMapping_fit_geo.m) to fit 7 parameters, including P0y, P0z, dety00, detz00, tilt_x, tilt_y and tilt_z.
  - After geometry fitting, you may repeat steps of 2-3 to redo/refine grain mapping results.<br>
### All the codes have been tested executable with a Matlab version 2020a or later.

## Run on linux cluster
Examples of bash file for submitting jobs to ESRF slurm cluster:
1) [go_SLURM_job_fullvol.slurm](https://github.com/haixingfang/GrainRecon/blob/main/go_SLURM_job_fullvol.slurm) shows an example to submit a job to a cluster equipped with only CPU for reconstructing a full volume.<br>
2) [go_SLURM_job_subvol_array.slurm](https://github.com/haixingfang/GrainRecon/blob/main/go_SLURM_job_subvol_array.slurm) shows an example to submit a job to a cluster equipped with only CPU for reconstructing all subvolumes partitioned from a full volume.<br>
3) [go_SLURM_gpu_single_job.slurm](https://github.com/haixingfang/GrainRecon/blob/main/go_SLURM_gpu_single_job.slurm) shows an example to submit a job to a cluster equipped with GPU for reconstructing a full volume. <br>
### command example for reconstructing the sample "virtual_Fe_100um_6grains", which is available in the folder [/Examples/virtual_Fe_100um_6grains/](https://github.com/haixingfang/GrainRecon/blob/main/Examples/virtual_Fe_100um_6grains)
```
sbatch go_SLURM_gpu_single_job.slurm
```

# Example dataset
Three examples are available in [./Examples](https://github.com/haixingfang/GrainRecon/blob/main/Examples):
1) a virtual Fe sample containing 6 grains in a magnified geometry; <br>
2) a virtual Fe sample containing 144 grains in both Laue-focusing and magnified geometries (exists only on Gitlab site); <br>
3) an experimental AlCu alloy sample in a magnified geometry (exists only on Gitlab site).<br>

# Note
By default the code is designed for reconstructing cubic crystals. For other crystal symmetries, You may need to modify some code lines related to the crystal symmetry name.

# License
This package is free to use, ditribute and adapt for non-commercial use only. See [LICENSE](https://github.com/haixingfang/GrainRecon/blob/main/LICENSE) for license rights and limitations (CC BY-NC 4.0).
A mirror site is accessible [here](https://github.com/haixingfang/GrainRecon) on Github.

# Reference
A manuscript has been prepared and now is under review:<br>
H. Fang, W. Ludwig, P. Lhuissier, Reconstruction algorithms for grain mapping by laboratory X-ray diffraction contrast tomography (in review).<br>
Please cite this article if you use or get inspired by the code presented here. <br>

# Know more about LabDCT or synchrotron DCT?
## Read the following books and articles.
[1] [Poulsen, H.F. (2004). Three-dimensional X-ray diffraction microscopy: mapping polycrystals and their dynamics, Springer, Berlin.](https://books.google.fr/books?hl=zh-CN&lr=&id=_jzrH20Qu6cC&oi=fnd&pg=PA1&dq=Three-dimensional+X-ray+diffraction+microscopy:+mapping+polycrystals+and+their+dynamics&ots=fuKB6aOUDR&sig=X1FLzGThZC5dBig_TmHRcPR34Jk&redir_esc=y#v=onepage&q=Three-dimensional%20X-ray%20diffraction%20microscopy%3A%20mapping%20polycrystals%20and%20their%20dynamics&f=false)<br>
[2] [Ludwig, W., Schmidt, S., Lauridsen, E.M. & Poulsen, H.F. (2008). J. Appl. Cryst. 41, 302-309.](https://onlinelibrary.wiley.com/doi/pdf/10.1107/S0021889808001684?casa_token=R34uKE0yZ-kAAAAA:nAWCkh8VEcvYkcdsX7gUqB3C05qQDH-5WrJ-OtSuBEiqf_iT1I3s2nCKz4sVOUSEvPYmzXJiOWmrBbH0)<br>
[3] [Ludwig, W., Reischig, P., King, A. Herbig, M., Lauridsen, E.M., Johnson, G., Marrow, T.J. & Buffière, J.Y. (2009). Rev. Sci. Instrum. 80, 033905.](https://aip.scitation.org/doi/full/10.1063/1.3100200?casa_token=P5TD352wKKgAAAAA:JQJrFf2zposYugxPD1u7j_TInetWxNG8cojaDD_Xd8VfJi4IyYkLGf5gXEv-m1YwWH49zBCS9WRO)<br>
[4] [King, A., Reischig, P., Adrien, J. & Ludwig, W. (2013). J. Appl. Cryst. 46, 1734-1740.](https://onlinelibrary.wiley.com/doi/pdf/10.1107/S0021889813022553?casa_token=qNuPs8Cl0HYAAAAA:cdd2pUDdX4zQnAXdeM47NNfu_A2KUeFLcCvSQL37allmTNCuks3_Uqq7idWahDsFgfliuTYttIfvfFPT)<br>
[5] [Van Aarle, W., Ludwig, W., King, A. & Penumadu, D. (2015). J. Appl. Cryst. 48, 334-343.](https://onlinelibrary.wiley.com/doi/full/10.1107/S1600576715000928?casa_token=2NJbHkPcSqAAAAAA:E8Y8bRglog_x8aa2csR4KwR4ElfHcs3AiV6fdhVwerqJ2jptIwxXW1p7Rfrq0HPf5OfFFHNalBfPoiq2)<br>
[6] [Bachmann, F., Bale, H., Gueninchault, N., Holzner, C. & Lauridsen, E.M. (2019). J. Appl. Cryst. 52, 643-651.](https://journals.iucr.org/j/issues/2019/03/00/nb5238/nb5238.pdf)<br>
[7] [Fang, H., Juul Jensen, D. & Zhang, Y. (2020). Acta Cryst. A 76, 652-663.](https://journals.iucr.org/a/issues/2020/06/00/iv5008/iv5008.pdf)<br>
[8] [Fang, H., Juul Jensen, D. & Zhang, Y. (2021). IUCrJ 8, 559-573.](https://journals.iucr.org/m/issues/2021/04/00/fc5052/index.html)<br>

## Contact emails: haixing.fang@grenoble-inp.fr, haixing0a@esrf.fr or haixingfang868@gmail.com


