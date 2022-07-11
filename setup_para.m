% set up all parameters: geometry, detector, sample, reconstruction, filefolders
% Nov 23, 2021

%%%%%%%%%%%%%% different samples
switch(SampleName)
    %%%%%%%%%%%%%%%%%%%%% AlCu8wt_middle_thinned 20210930 LabDCT
    case 'AlCu8wt_middle_thinned_0930'
    %     input_al;
    [cell,sgno,atomparam,space_group_IT_number]=input_al_fun();
    simap_data_flag=1; % using simap data = 1; otherwise = 0
    sample_no=2; % sample number
    fitted_geo_already=1;  % fitted geometry already? 1-yes; 0-no
    not_remove_spot=0; % flag for not to remove assigned spots from the spot list

    % load experimental parameters: source, detector and distances
    setup_exp;
    L=Lsam2sou+Lsam2det;
    RotX=[1 0 0; 0 cosd(tilt_x) -sind(tilt_x); 0 sind(tilt_x) cosd(tilt_x)];
    RotY=[cosd(tilt_y) 0 sind(tilt_y);0 1 0;-sind(tilt_y) 0 cosd(tilt_y)];
    RotZ=[cosd(tilt_z) -sind(tilt_z) 0;sind(tilt_z) cosd(tilt_z) 0;0 0 1];
    RotDet=RotX*RotY*RotZ;

    % tomo file
    FileFolder='./Examples/AlCu_8wt_middle_thinned_0930';
    fname_prefix=FileFolder(3:end);
    if ~exist(OutputFolder,'dir')
        mkdir(OutputFolder); % to store output files
    end
    h5Folder_tomo=FileFolder;
    h5FileName_tomo='tomo_2021_09_30_AlCu_8wt_middle_thinned.h5';

    % Spots file
    SpotsFile='Spots_AlCu_8wt_middle_thinned_0930.mat';

    % indexing criteria
    TrustComp=0.85; % trust completeness
    minComp=0.3;   % minimum completeness
    drop_off=0.05;  % drop-off value for growing indexed region
    minEucDis=0.5*mean([pixelysize pixelzsize])*1; % minimum tolerate euclidien distance [mm]
    maxD=3; % maximum acceptable completeness weighted center difference, it needs update if larger than this value [pixel]
    maxDmedian=19.9; % maximum acceptable median distance betwee forward and exp spots [pixel]
    if simap_data_flag==1
        S=[1 0 0;0 1 0;0 0 1];
    else
        S=[1 0 0;0 -1 0;0 0 1];
    end
    
    %%%%% Fe_100um_11_11_simu
    case 'simu_Fe'
%     input_fe;
    [cell,sgno,atomparam,space_group_IT_number]=input_fe_fun();
    simap_data_flag=1; % using simap data = 1; otherwise = 0
    if simap_data_flag==1
        sample_no=3; % sample number
    else
        sample_no=2; % sample number
    end
    fitted_geo_already=1;  % fitted geometry already? 1-yes; 0-no
    not_remove_spot=0; % flag for not to remove assigned spots from the spot list

    % load experimental parameters: source, detector and distances
    setup_exp;
    L=Lsam2sou+Lsam2det;
    RotX=[1 0 0; 0 cosd(tilt_x) -sind(tilt_x); 0 sind(tilt_x) cosd(tilt_x)];
    RotY=[cosd(tilt_y) 0 sind(tilt_y);0 1 0;-sind(tilt_y) 0 cosd(tilt_y)];
    RotZ=[cosd(tilt_z) -sind(tilt_z) 0;sind(tilt_z) cosd(tilt_z) 0;0 0 1];
    RotDet=RotX*RotY*RotZ;

    % tomo file
    FileFolder='./Examples/Fe_100um_11_11_simu';
    fname_prefix=FileFolder(3:end);
    if ~exist(OutputFolder,'dir')
        mkdir(OutputFolder); % to store output files
    end
    h5Folder_tomo=FileFolder;
    if simap_data_flag~=1
        h5FileName_tomo='tomo_Fe_100um_11_11_simu.h5';
    else
        h5FileName_tomo='tomo_Fe_100um_11_11_simu_simap.h5';
    end

    % Spots file
    if simap_data_flag==1
         SpotsFile='Spots_Fe_100um_simu_geo_simap_s2_no_plus_4hkl.mat'; % May 16, 2022
%         SpotsFile='Spots_Fe_100um_simu_geo_simap_s2_181projs.mat'; % June 23, 2022
    else
         SpotsFile='Spots_Fe_100um_11_11_simu_LoG.mat';
    end

    % indexing criteria
    if simap_data_flag==1
        TrustComp=0.85; % trust completeness
        minComp=0.45;   % minimum completeness
        maxDmedian=15;
    else
        TrustComp=0.85; % trust completeness
        minComp=0.55;   % minimum completeness
        maxDmedian=10;
    end
    drop_off=0.02;  % drop-off value for growing indexed region
    minEucDis=0.5*mean([pixelysize pixelzsize])*1; % minimum tolerate euclidien distance [mm]
    maxD=3; % maximum acceptable completeness weighted center difference, it needs update if larger than this value [pixel]
    if simap_data_flag==1
        S=[1 0 0;0 1 0;0 0 1];
    else
        S=[1 0 0;0 -1 0;0 0 1];
    end
    
case 'virtual_Fe_100um_6grains'
    [cell,sgno,atomparam,space_group_IT_number]=input_fe_fun();
    simap_data_flag=1; % using simap data = 1; otherwise = 0
    sample_no=4; % sample number
    fitted_geo_already=1;  % fitted geometry already? 1-yes; 0-no
    not_remove_spot=1; % flag for not to remove assigned spots from the spot list

    % load experimental parameters: source, detector and distances
    setup_exp;
    L=Lsam2sou+Lsam2det;
    RotX=[1 0 0; 0 cosd(tilt_x) -sind(tilt_x); 0 sind(tilt_x) cosd(tilt_x)];
    RotY=[cosd(tilt_y) 0 sind(tilt_y);0 1 0;-sind(tilt_y) 0 cosd(tilt_y)];
    RotZ=[cosd(tilt_z) -sind(tilt_z) 0;sind(tilt_z) cosd(tilt_z) 0;0 0 1];
    RotDet=RotX*RotY*RotZ;

    % tomo file
    FileFolder='./Examples/virtual_Fe_100um_6grains';
    fname_prefix=FileFolder(3:end);
    if ~exist(OutputFolder,'dir')
        mkdir(OutputFolder); % to store output files
    end
    h5Folder_tomo=FileFolder;
    if simap_data_flag~=1
         h5FileName_tomo='tomo_virtual_Fe_100um_6grains.h5';
    else
         h5FileName_tomo='tomo_virtual_Fe_100um_6grains_simap.h5'; % rotate by 180 deg counter-clockwise
    end

    % Spots file
    if simap_data_flag==1
        SpotsFile='Spots_virtual_Fe_100um_6grains_simu_geo_simap_5hkl.mat';
    else
        error('Spots file is not available yet');
    end

    % indexing criteria
    TrustComp=0.85; % trust completeness
    minComp=0.55;   % minimum completeness
    drop_off=0.02;  % drop-off value for growing indexed region
    minEucDis=0.5*mean([pixelysize pixelzsize])*1; % minimum tolerate euclidien distance [mm]
    maxD=3; % maximum acceptable completeness weighted center difference, it needs update if larger than this value [pixel]
    maxDmedian=5; % maximum acceptable median distance betwee forward and exp spots [pixel]
    if simap_data_flag==1
        S=[1 0 0;0 1 0;0 0 1];
    else
        S=[1 0 0;0 -1 0;0 0 1];
    end
end
hklnumber=3; % maximum is 10, recommended be at least >= 3, 4
% iter_end=10; % maximum iterations which determines the finess of gridding for indexing seeds
sprintf('hklnumber = %.0f, C_trust = %.2f, C_min = %.2f, drop_off = %.2f, minEucDis = %.3f mm, maxD = %.1f pixel and maxDmedian = %.1f pixel', ...
    hklnumber,TrustComp,minComp,drop_off,minEucDis,maxD,maxDmedian)
B=FormB(cell);
V = cellvolume(cell); % [Angs^3]


