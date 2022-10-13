% please use get_geometry.m to acquire the gemotry information

% experimental parameters
if simap_data_flag~=1
    if sample_no==1
        % Al, Lss = 14 mm, Lsd = 14 mm, 2017-02-24
        Lsam2sou = 13.96; % sample-to-source distance [mm] [11 14]
        Lsam2det = 14.01; % sample-to-detector distance [mm][11 14 18 30]
        dety0 = 1016;	% beamcenter, y in pixel coordinatees
        detz0 = 1016;	% beamcenter, z in pixel coordinatees
        detysize = 2032;      % detector y size [pixels]
        detzsize = 2032;      % detector z size [pixels]
        pixelysize=3.3623e-3; % pixel size of the projection image [mm]
        pixelzsize=3.3623e-3; % pixel size of the projection image [mm]
        BeamStopY=[585 1460];
        BeamStopZ=[545 1480];
    elseif sample_no==2
        % Simu_Fe_100um_11_11
        Lsam2sou = 11.0; % sample-to-source distance [mm] [11 14]
        Lsam2det = 11.0; % sample-to-detector distance [mm][11 14 18 30]
        dety0 = 1016;	% beamcenter, y in pixel coordinatees
        detz0 = 1016;	% beamcenter, z in pixel coordinatees
        detysize = 2032;      % detector y size [pixels]
        detzsize = 2032;      % detector z size [pixels]
        pixelysize=3.3623e-3; % pixel size of the projection image [mm]
        pixelzsize=3.3623e-3; % pixel size of the projection image [mm]
        BeamStopY=[545 1506];
        BeamStopZ=[536 1488];
    elseif sample_no==3
        % Iron_CR50_dtu_13_13
        Lsam2sou = 13.022; % sample-to-source distance [mm] [11 14]
        Lsam2det = 12.958; % sample-to-detector distance [mm][11 14 18 30]
        dety0 = 1016;	% beamcenter, y in pixel coordinatees
        detz0 = 1016;	% beamcenter, z in pixel coordinatees
        detysize = 2032;      % detector y size [pixels]
        detzsize = 2032;      % detector z size [pixels]
        pixelysize=3.3623e-3; % pixel size of the projection image [mm]
        pixelzsize=3.3623e-3; % pixel size of the projection image [mm]
        BeamStopY=[560 1455];
        BeamStopZ=[545 1460];
    end
    aquisition_center=[0 0 0];
    % source offset
    P0=[-Lsam2sou 0 0]; % [mm]
    P0=P0-aquisition_center;
    P0y=P0(2); 
    P0z=P0(3);
	RotAxisOffset=0;
    % detector offset and tilt
    det_center=[Lsam2det 0 0];
    det_normal=[1 0 0];
	detdiru=[0 -1 0]; % pointing to the right
	detdirv=[0 0 -1]; % pointing downwards
	[RotDet,tilt_x,tilt_y,tilt_z]=get_det_tilt(det_normal,detdiru);
    dety00=det_center(2)-aquisition_center(2); % [mm]
    detz00=det_center(3)-aquisition_center(3); % [mm]
    ExpTime=500; % exposure time [s]
else
    if sample_no==1
    % CCD_AlCu8wt%_DCT_2021_07_22_121projs
    % D:\ExpData_process\CCD_2021_07_22_AlCu8wt\DCT_1_6_all
    Lsam2sou = 7.1810; % sample-to-source distance [mm]
    Lsam2det = 57.3852; % sample-to-detector distance [mm]
    dety0 = 1020;	% beamcenter, y in pixel coordinatees
    detz0 = 1020;	% beamcenter, z in pixel coordinatees
    detysize = 2040;      % detector y size [pixels]
    detzsize = 2040;      % detector z size [pixels]
    pixelysize=24e-3; % pixel size of the projection image [mm]
    pixelzsize=24e-3; % pixel size of the projection image [mm]
    BeamStopY=[531 1670];
    BeamStopZ=[450 1580];
    BeamStopCenter=[1107 996];
    BeamStopRadius=575;
    aquisition_center=[0 0 -110.1430]; % from tomo file: Y=>Z, X=>-Y, Z=>-X (left: simap setup; right: conventional system)
    % source offset
    P0=[-7.1779 0.0294 -110.1430]; % from get_geometry.m [mm]
    P0=P0-aquisition_center;
    P0y=P0(2);
    P0z=P0(3);
    % detector offset and tilt
    det_center=[57.3633 -0.2345 -108.3650]; % from get_geometry.m
    det_normal=[0.9999 -0.0041 -0.0096];   % from get_geometry.m
	detdiru = [-0.004080 -0.999992 0.000041];  % pointing to the right, from get_geometry.m
	detdirv = [-0.009640 -0.000008 -0.999953]; % pointing downwards, from get_geometry.m
	[RotDet,tilt_x,tilt_y,tilt_z]=get_det_tilt(det_normal,detdiru);
    dety00=det_center(2)-aquisition_center(2); % [mm]
    detz00=det_center(3)-aquisition_center(3); % [mm]
    ExpTime=60*6; % exposure time [s]

    % when first time indexing, check this
    fitted_geo_already=1;
    if fitted_geo_already==1
        ParaFit=[7.1253   57.3497   -0.2219    1.7374    0.5475    0.0131    0.3556];
%         ParaFit=[7.1501   57.3800   -0.1424    1.7135    0.4323   -0.0210   -0.2582];
        Lsam2sou = ParaFit(1);
        Lsam2det = ParaFit(2);
        dety00 = ParaFit(3);
        detz00 = ParaFit(4);
        tilt_x = ParaFit(5);
        tilt_y = ParaFit(6);
        tilt_z = ParaFit(7);
		if length(ParaFit)==10
			P0y=ParaFit(8);
			P0z=ParaFit(9);
			RotAxisOffset=ParaFit(10);
		end
    end
    
    elseif sample_no==2
    % CCD_AlCu8wt%_DCT_middle_thinned_2021_09_30_121projs
    % D:\ExpData_process\CCD_2021_09_30_AlCu_8wt_middle_thinned\3_DCT
    Lsam2sou = 6.60779; % sample-to-source distance [mm]
    Lsam2det = 53.0147; % sample-to-detector distance [mm]
    dety0 = 1020;	% beamcenter, y in pixel coordinatees
    detz0 = 1020;	% beamcenter, z in pixel coordinatees
    detysize = 2040;      % detector y size [pixels]
    detzsize = 2040;      % detector z size [pixels]
    pixelysize=24e-3; % pixel size of the projection image [mm]
    pixelzsize=24e-3; % pixel size of the projection image [mm]
    BeamStopY=[534 1670];
    BeamStopZ=[430 1590];
    BeamStopCenter=[1095 999];
    BeamStopRadius=574;
    aquisition_center=[0 0 -121.433]; % from tomo file: Y=>Z, X=>-Y, Z=>-X (left: simap setup; right: conventional system)
    % source offset
    P0=[-6.6034 0.0276 -121.4330]; % from get_geometry.m [mm]
    P0=P0-aquisition_center;
    P0y=P0(2);
    P0z=P0(3);
    % detector offset and tilt
    det_center=[52.9931 -0.2163 -119.684]; % from get_geometry.m
    det_normal=[0.999945 -0.004076 -0.009639];   % from get_geometry.m
	detdiru = [-0.004077 -0.999992 -0.000041];
	detdirv = [-0.009639 0.000079 -0.999954]
	[RotDet,tilt_x,tilt_y,tilt_z]=get_det_tilt(det_normal,detdiru);
    dety00=det_center(2)-aquisition_center(2); % [mm]
    detz00=det_center(3)-aquisition_center(3); % [mm]
    ExpTime=60*6; % exposure time [s]
    % when first time indexing, check this
    if fitted_geo_already==1
        ParaFit=[6.2018   52.8844   -0.2641    1.6125    0.0255    0.6953    0.3482]; % 4.47 pixel, Nov 25, fitted by 8 grains
        Lsam2sou = ParaFit(1);
        Lsam2det = ParaFit(2);
        dety00 = ParaFit(3);
        detz00 = ParaFit(4);
        tilt_x = ParaFit(5);
        tilt_y = ParaFit(6);
        tilt_z = ParaFit(7);
		if length(ParaFit)==10
			P0y=ParaFit(8);
			P0z=ParaFit(9);
			RotAxisOffset=ParaFit(10);
		end
    end
    elseif sample_no==3
        % simu_Fe using the geometry from CCD_AlCu8wt%_DCT_middle_thinned_2021_09_30_121projs
        % D:\ExpData_process\CCD_2021_09_30_AlCu_8wt_middle_thinned\3_DCT
        Lsam2sou = 6.60779; % sample-to-source distance [mm]
        Lsam2det = 53.0147; % sample-to-detector distance [mm]
        dety0 = 1020;	% beamcenter, y in pixel coordinatees
        detz0 = 1020;	% beamcenter, z in pixel coordinatees
        detysize = 2040;      % detector y size [pixels]
        detzsize = 2040;      % detector z size [pixels]
        pixelysize=24e-3; % pixel size of the projection image [mm]
        pixelzsize=24e-3; % pixel size of the projection image [mm]
        BeamStopY=[534 1670];
        BeamStopZ=[430 1590];
        BeamStopCenter=[1095 999];
        BeamStopRadius=574;
        aquisition_center=[0 0 -121.433]; % from tomo file: Y=>Z, X=>-Y, Z=>-X (left: simap setup; right: conventional system)
        % source offset
        P0=[-6.6034 0.0276 -121.4330]; % from get_geometry.m [mm]
        P0=P0-aquisition_center;
        P0y=P0(2);
        P0z=P0(3);
        % detector offset and tilt
        det_center=[52.9931 -0.2163 -119.684]; % from get_geometry.m
		det_normal=[0.999945 -0.004076 -0.009639];   % from get_geometry.m
		detdiru = [-0.004077 -0.999992 -0.000041];
		detdirv = [-0.009639 0.000079 -0.999954]
		[RotDet,tilt_x,tilt_y,tilt_z]=get_det_tilt(det_normal,detdiru);
        dety00=det_center(2)-aquisition_center(2); % [mm]
        detz00=det_center(3)-aquisition_center(3); % [mm]
        ExpTime=60*6; % exposure time [s]
        % when first time indexing, check this
        if fitted_geo_already==1
            ParaFit=[6.1361   52.8971   -0.2426    1.5924    0.0036    0.6366    0.3528]; % true geometry
            Lsam2sou = ParaFit(1);
            Lsam2det = ParaFit(2);
            dety00 = ParaFit(3);
            detz00 = ParaFit(4);
            tilt_x = ParaFit(5);
            tilt_y = ParaFit(6);
            tilt_z = ParaFit(7);
			if length(ParaFit)==10
                P0y=ParaFit(8);
                P0z=ParaFit(9);
                RotAxisOffset=ParaFit(10);
            end
        end
    elseif sample_no==4
        % simu_Fe using the geometry from CCD_AlCu8wt%_DCT_middle_thinned_2021_09_30_121projs
        % D:\Documents\Matlab\GrainRecon\virtual_Fe_100um_6grains
        Lsam2sou = 6.60779; % sample-to-source distance [mm]
        Lsam2det = 53.0147; % sample-to-detector distance [mm]
        dety0 = 1020;	% beamcenter, y in pixel coordinatees
        detz0 = 1020;	% beamcenter, z in pixel coordinatees
        detysize = 2040;      % detector y size [pixels]
        detzsize = 2040;      % detector z size [pixels]
        pixelysize=24e-3; % pixel size of the projection image [mm]
        pixelzsize=24e-3; % pixel size of the projection image [mm]
        BeamStopY=[534 1670];
        BeamStopZ=[430 1590];
        BeamStopCenter=[1095 999];
        BeamStopRadius=574;
        aquisition_center=[0 0 -121.433]; % from tomo file: Y=>Z, X=>-Y, Z=>-X (left: simap setup; right: conventional system)
        % source offset
        P0=[-6.6034 0.0276 -121.4330]; % from get_geometry.m [mm]
        P0=P0-aquisition_center;
        P0y=P0(2);
        P0z=P0(3);
        % detector offset and tilt
        det_center=[52.9931 -0.2163 -119.684]; % from get_geometry.m
		det_normal=[0.999945 -0.004076 -0.009639];   % from get_geometry.m
		detdiru = [-0.004077 -0.999992 -0.000041];
		detdirv = [-0.009639 0.000079 -0.999954]
		[RotDet,tilt_x,tilt_y,tilt_z]=get_det_tilt(det_normal,detdiru);
        dety00=det_center(2)-aquisition_center(2); % [mm]
        detz00=det_center(3)-aquisition_center(3); % [mm]
        ExpTime=60*6; % exposure time [s]
        % when first time indexing, check this
        if fitted_geo_already==1
            ParaFit=[6.1361   52.8971   -0.2426    1.5924    0.0036    0.6366    0.3528]; % true geometry
            Lsam2sou = ParaFit(1);
            Lsam2det = ParaFit(2);
            dety00 = ParaFit(3);
            detz00 = ParaFit(4);
            tilt_x = ParaFit(5);
            tilt_y = ParaFit(6);
            tilt_z = ParaFit(7);
			if length(ParaFit)==10
                P0y=ParaFit(8);
                P0z=ParaFit(9);
                RotAxisOffset=ParaFit(10);
            end
        end
    end
end



    
    

