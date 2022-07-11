% call mex_forward_comp.c
% return Nr_simu, Nr_match, dis_median
% mex, April 8, 2022

function [Nr_simu,Nr_intersect,dis_median] = forward_comp_mex(pos,U,S,B,Ahkl,nrhkl, ...
                        proj_bin_bw,rot_angles,RotDet,Lsam2sou,Lsam2det, ...
                        P0y,P0z,dety00,detz00,pixelysize,pixelzsize,detysize,detzsize, ...
                        BeamStopY,BeamStopZ,thetamax,lambda_min,lambda_max)
% % for testing
% pos=pos_indexing;
hkl = [Ahkl(1:nrhkl,1) Ahkl(1:nrhkl,2) Ahkl(1:nrhkl,3)]';

nr_hkl = size(hkl,2);
nr_rot = length(rot_angles);
U_mex = reshape(U',1,9);
S_mex = reshape(S',1,9);
B_mex = reshape(B',1,9);
hkl_mex = reshape(hkl,[1 size(hkl,1)*size(hkl,2)]);
RotDet_mex = reshape(RotDet',1,9);
proj_bin_bw_mex = double(reshape(proj_bin_bw,1,[])); % new version using 1D input
param_mex = double([nr_hkl nr_rot Lsam2sou Lsam2det P0y P0z ...
                    dety00 detz00 pixelysize pixelzsize detysize detzsize BeamStopY ...
                    BeamStopZ thetamax lambda_min lambda_max]);

Output_mex = mex_forward_comp(pos, U_mex, S_mex, B_mex, hkl_mex, proj_bin_bw_mex, ...
                        rot_angles, RotDet_mex, param_mex);
Nr_simu = Output_mex(1);
Nr_intersect = Output_mex(2);
dis_median = Output_mex(4);

    
