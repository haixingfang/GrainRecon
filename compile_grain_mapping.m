% compile mex cuda functions
% April 14, 2022

function compile_grain_mapping()

mexcuda cuda_forward_Gt.cu;
mexcuda cuda_forward_comp.cu;
mexcuda cuda_grow_indexed_region.cu;
mexcuda cuda_forward_comp_compete.cu;

fprintf('Successfully compiled cuda functions.\n');
fprintf(['Note: you may need to modify the pre-assigned size for "float Gt_match" in "cuda_forward_Gt.cu".\n' ...
    'It is recomended to set it as "float Gt_match[nr_spots*7]" where nr_spots = nr_projs*maximum spots per projection. \n']);
