% % grain verification and pairing spots
% % shapes of forward simulated spots are also computed
% % September 16, 2021
function [SpotNr_simu,SpotNr_obs,SpotsPair_all,GrainIndex_all,SubGrain,rot_angles,SmallGrID, ...
    SpotNr_simu_hkl,SpotNr_obs_hkl]=Forward_simu_spots_exp(DS,Rsample,RecVolumePixel, ...
    tomo_scale,ExpTime,atomparam,proj,Spots,rot_start,rot_step,S,B,Ahkl,nrhkl,hkl_square, ...
    Energy,lambda,V,K1,I0E,RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,dety00,detz00,P0y,P0z, ...
    pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ,simap_data_flag, ...
    OutputFolder,rot_angles,simu_grainno)
% for testing
% DS=DS_fit;

if nargin<=41
    simu_grainno=1:length(DS.SeedID);
end
L = Lsam2sou + Lsam2det;

% add on Nov 25, 2021, to calculate the grain centroid based on weighted completeness
for k=simu_grainno%1:length(DS.SeedID)
    grainID=DS.SeedID(k);
    ind=[];
    [ind(:,1),ind(:,2),ind(:,3)] = ind2sub(size(DS.GIDvol),find(DS.GIDvol == grainID));
    ind1=find(DS.GIDvol == grainID);
    id(k,:)=sum((DS.CompVol(ind1).*ind))/sum(DS.CompVol(ind1));
    
    grain_centroid(k,1)=((id(k,1)+RecVolumePixel(1,1)-1)-tomo_scale.Dimension(1)/2).*DS.VoxSize(1)+DS.Center(1);
    grain_centroid(k,2)=((id(k,2)+RecVolumePixel(2,1)-1)-tomo_scale.Dimension(2)/2).*DS.VoxSize(2)+DS.Center(2); % centroid coordinate y
    grain_centroid(k,3)=((id(k,3)+RecVolumePixel(3,1)-1)-tomo_scale.Dimension(3)/2).*DS.VoxSize(3)+DS.Center(3);
    if simap_data_flag==1
        grain_centroid(k,1)=-grain_centroid(k,1); % [mm]
        grain_centroid(k,2)=-grain_centroid(k,2); % [mm]
    end
end

% voxelize individual grains => "subgrain"
SmallGrID=[];
for jj=simu_grainno%1:length(DS.SeedID)
    yy=DS.GIDvol==DS.SeedID(jj);
    yy=uint8(yy);
    if DS.nVox(jj)<4000
        Scaling=1;
    elseif DS.nVox(jj)<8000
        Scaling=0.5;
    else
        Scaling=0.25;
    end
    zz=imresize3(yy,Scaling,'nearest'); % function valid after 2017a
    Binning=size(yy)./size(zz);
    clear id;
    [id(:,1), id(:,2), id(:,3)]=ind2sub(size(zz),find(zz==1));
    if ~isempty(id)
        for k=1:length(id(:,1))
            SubGrain{jj}(k,1)=k;
            SubGrain{jj}(k,2)=((id(k,1)+RecVolumePixel(1,1)/Binning(1)-1)-tomo_scale.Dimension(1)/Binning(1)/2).*DS.VoxSize(1)*Binning(1)+DS.Center(1); % centroid coordinate x
            SubGrain{jj}(k,3)=((id(k,2)+RecVolumePixel(2,1)/Binning(2)-1)-tomo_scale.Dimension(2)/Binning(2)/2).*DS.VoxSize(2)*Binning(2)+DS.Center(2); % centroid coordinate y
            SubGrain{jj}(k,4)=((id(k,3)+RecVolumePixel(3,1)/Binning(3)-1)-tomo_scale.Dimension(3)/Binning(3)/2).*DS.VoxSize(3)*Binning(3)+DS.Center(3); % centroid coordinate z
            SubGrain{jj}(k,5)=DS.VoxSize(1)*DS.VoxSize(2)*DS.VoxSize(3)*Binning(1)*Binning(2)*Binning(3); % volume
            SubGrain{jj}(k,6)=2*(3*SubGrain{jj}(k,5)/(4*pi))^(1/3); % EqDiameter
            SubGrain{jj}(k,7)=2*(3*DS.nVox(jj)*DS.VoxSize(1)*DS.VoxSize(2)*DS.VoxSize(3)/(4*pi))^(1/3); % EqDiameter of the parent grain
        end
    else
        clear id;
        [id(:,1), id(:,2), id(:,3)]=ind2sub(size(yy),find(yy==1));
        for k=1:length(id(:,1))
            SubGrain{jj}(k,1)=k;
            SubGrain{jj}(k,2)=((id(k,1)+RecVolumePixel(1,1)-1)-tomo_scale.Dimension(1)/2).*DS.VoxSize(1)+DS.Center(1); % centroid coordinate x
            SubGrain{jj}(k,3)=((id(k,2)+RecVolumePixel(2,1)-1)-tomo_scale.Dimension(2)/2).*DS.VoxSize(2)+DS.Center(2); % centroid coordinate y
            SubGrain{jj}(k,4)=((id(k,3)+RecVolumePixel(3,1)-1)-tomo_scale.Dimension(3)/2).*DS.VoxSize(3)+DS.Center(3); % centroid coordinate z
            SubGrain{jj}(k,5)=DS.VoxSize(1)*DS.VoxSize(2)*DS.VoxSize(3); % volume
            SubGrain{jj}(k,6)=2*(3*SubGrain{jj}(k,5)/(4*pi))^(1/3); % EqDiameter
            SubGrain{jj}(k,7)=2*(3*DS.nVox(jj)*DS.VoxSize(1)*DS.VoxSize(2)*DS.VoxSize(3)/(4*pi))^(1/3); % EqDiameter of the parent grain
        end
        SmallGrID=[SmallGrID;jj];        
    end
    if ~isempty(id)
        if ~isempty(SubGrain{jj})
            if simap_data_flag==1
                SubGrain{jj}(:,2)=-SubGrain{jj}(:,2); % [mm]
                SubGrain{jj}(:,3)=-SubGrain{jj}(:,3); % [mm]
            end
        end
    end
    jj;
end

% estimate sample diameter by assuming the sample has approximately a cylinder shape
if ~exist('Rsample','var') 
    Rsample=0.4; %[mm]
end
if ~exist('ExpTime','var')
    ExpTime=60*6;
end
% check DA_cmp and TFT_cmp folders exist
if ispc
    if ~exist(strcat(OutputFolder,'\TFT_cmp'), 'dir')
       mkdir(strcat(OutputFolder,'\TFT_cmp')); % TFT_cmp folder is to store each output projection
       direc = strcat(OutputFolder,'\TFT_cmp');  % save frames in this directory
    else
       direc = strcat(OutputFolder,'\TFT_cmp');
    end
    if ~exist(strcat(OutputFolder,'\DA_cmp'), 'dir')
       mkdir(strcat(OutputFolder,'\DA_cmp')); % DA_cmp folder is to store data record for each projection
       direc2 = strcat(OutputFolder,'\DA_cmp');
    else
       direc2 = strcat(OutputFolder,'\DA_cmp');
    end
elseif isunix || ismac
    if ~exist(strcat(OutputFolder,'/TFT_cmp'), 'dir')
       mkdir(strcat(OutputFolder,'/TFT_cmp')); % TFT_cmp folder is to store each output projection
       direc = strcat(OutputFolder,'/TFT_cmp');  % save frames in this directory
    else
       direc = strcat(OutputFolder,'/TFT_cmp');
    end
    if ~exist(strcat(OutputFolder,'/DA_cmp'), 'dir')
       mkdir(strcat(OutputFolder,'/DA_cmp')); % DA_cmp folder is to store data record for each projection
       direc2 = strcat(OutputFolder,'/DA_cmp');
    else
       direc2 = strcat(OutputFolder,'/DA_cmp');
    end
end

prefix = 'proj'; % prefix of frame names, default is 'frame'
% basic grain information
grains=length(DS.SeedID); % number of grains
grainvolume=DS.nVox*DS.VoxSize(1)*DS.VoxSize(2)*DS.VoxSize(3)*1e9; % [um^3]
grainsize=2*(3*grainvolume/(4*pi)).^(1/3); % equivalent diameter of grain [um]
euler_grains=DS.EulerZXZ;% Euler angles [deg]: ([0 2*pi], [0 pi], [0 2*pi])

% read transmission data and CsI scintillator data, add on June 22, 2020
atomparam_atomno=atomparam.atomno;
[Transmission, rou]=ReadTransData(atomparam_atomno); % [(-), g/cm^3]
[CsI, Swank]=ReadCsI();
hkl_color=[255 0 0;0 255 0;0 0 255;255 255 0;0 255 255;255 0 255;0 128 0;128 0 128;0 0 128;255 140 0; ...
    128 0 0;0 128 128;255 215 0;240 230 140;255 69 0;255 140 0];
% hkl color: red, green, blue, yellow, cyan, magenta, olive,
% purple, navy, dark orange
% maroon, teal, gold, khaki,orange red, dark orange
hkl = [Ahkl(1:nrhkl,1) Ahkl(1:nrhkl,2) Ahkl(1:nrhkl,3)]';

counter=0;
SpotsPair_all=[];
GrainIndex_all=[];
% for rot = rot_start:rot_step:rot_end
for rot = rot_angles
    counter=counter+1;
    rot_number=(rot-rot_start)/rot_step+1; % recording number of rotations
    % experimental LabDCT projections
    imdata=uint16(proj{rot_number});

    omega=rot*pi/180; % [rad]
    Omega=[cos(omega) -sin(omega) 0;sin(omega) cos(omega) 0;0 0 1]; 
    graininfo = zeros(abs(grains),16);
    A=[];
    % Generate orientations of the grains and loop over all grains
    for grainno = simu_grainno%1:abs(grains)
        phi1 = euler_grains(grainno,1)*pi/180;
        Phi = euler_grains(grainno,2)*pi/180;
        phi2 = euler_grains(grainno,3)*pi/180;
        U = euler2u(phi1,Phi,phi2);
        graininfo(grainno,1:6) = [grainno grainsize(grainno) grainvolume(grainno) phi1*180/pi Phi*180/pi phi2*180/pi];
        graininfo(grainno,7:15) = [U(1,1) U(1,2) U(1,3) U(2,1) U(2,2) U(2,3) U(3,1) U(3,2) U(3,3)];
        graininfo(grainno,16)=length(SubGrain{grainno}(:,1)); % number of 3D cells for calculation
        reshape(graininfo(grainno,7:15),3,3)';
        
        % Calculate matrix A with (1:totalnr, 2:grain, 3:refno, 4-6:h,k,l, 7:F^2, 8:phi1, 9:PHI, 10:phi2,
        % 11-13:Gw(1),Gw(2),Gw(3), 14:omega, 15:2theta, 16:eta, 17:dety, 18:detz, 19:Lorentz, 20:Polarization, 21:Int)
        % Gw is the G-vector in the omega-system (w=0)
        % Gt is the G-vector in the tilted system (identical to the lab-system except for the tilt of sample stage)
        % All angles in A are in degrees
        pos=SubGrain{grainno}(:,2:4); % [mm]
        SamposW_all=Omega*S*pos';
        center_all = [L*ones(1,size(SamposW_all,2)); (SamposW_all(2,:)-P0y)*L./(Lsam2sou+SamposW_all(1,:)); ...
            (SamposW_all(3,:)-P0z)*L./(Lsam2sou+SamposW_all(1,:))]'; % sample center projected to the position of the detector
        alpha_all = atan(sqrt((SamposW_all(2,:)-P0y).^2+(SamposW_all(3,:)-P0z).^2)./(Lsam2sou+SamposW_all(1,:)));
        grainpos_all = [Lsam2sou+SamposW_all(1,:); SamposW_all(2,:)-P0y; SamposW_all(3,:)-P0z]';

        Gw = S*U*B*hkl;
        Gt=Omega*Gw;
        v1 = [zeros(1,size(hkl,2));Gt(2,:);Gt(3,:)];
        Glen = (Gt(1,:).^2 + Gt(2,:).^2 + Gt(3,:).^2).^0.5;
        nr=1;
        nrefl = 1;
        SubA{grainno}=[];
        for subgrainno=1:length(SubGrain{grainno}(:,1))
            SamposW=SamposW_all(:,subgrainno);
            center=center_all(subgrainno,:);
            alpha=alpha_all(subgrainno);
            grainpos=grainpos_all(subgrainno,:);
            beta = acos(dot((ones(size(hkl))'.*grainpos/norm(grainpos))',Gt./Glen)); % [rad]

            theta = beta-pi/2;
            sintth = sin(2*theta);
            costth = cos(2*theta);
            d = 1./Glen*2*pi;
            lambdahkl = 2 * d .*sin(theta);
            Energy_hkl=12.398./lambdahkl; % [keV]

            phix = acos(dot(v1./sqrt(v1(1,:).^2+v1(2,:).^2+v1(3,:).^2),(ones(size(hkl))'.*center/norm(center))'));
            phiy = phix-2*theta;
            L2 = (Lsam2det-SamposW(1))/cos(alpha);
            diffvec = L2*sintth./sin(phiy); % [mm]
            konst=sqrt(Gt(2,:).^2+Gt(3,:).^2);

            dety22 = (center(2)+ (diffvec.*Gt(2,:)./konst)); % dety [mm]
            detz22 = (center(3)+ (diffvec.*Gt(3,:)./konst)); % detz [mm]

            K_out_unit = ([ones(1,size(hkl,2))*Lsam2det;dety22;detz22]-ones(1,size(hkl,2)).*SamposW) ...
                ./((Lsam2det-SamposW(1)).^2+(dety22-SamposW(2)).^2+(detz22-SamposW(3)).^2).^(1/2);
            t = (RotDet(1,1)*(Lsam2det-SamposW(1))-RotDet(2,1)*(dety00-SamposW(2))-RotDet(3,1)*(detz00-SamposW(3)))./ ...
            (RotDet(1,1)*K_out_unit(1,:)+RotDet(2,1)*K_out_unit(2,:)+RotDet(3,1)*K_out_unit(3,:));

            dety22 = [RotDet(1,2) RotDet(2,2) RotDet(3,2)]*(t.*K_out_unit+[SamposW(1)-Lsam2det SamposW(2)-dety00 SamposW(3)-detz00]');
            detz22 = [RotDet(1,3) RotDet(2,3) RotDet(3,3)]*(t.*K_out_unit+[SamposW(1)-Lsam2det SamposW(2)-dety00 SamposW(3)-detz00]');
%             dety = -round(dety22/pixelysize-0.5)+dety0; % [pixel]
%             detz = -round(detz22/pixelzsize-0.5)+detz0; % [pixel]
            dety = -round(dety22/pixelysize)+dety0; % [pixel]
            detz = -round(detz22/pixelzsize)+detz0; % [pixel]

            select3=find(beta > pi/2 & beta < (90+thetamax*4)/180*pi & ...
                lambdahkl > lambda_min & lambdahkl < lambda_max & ...
                dety>=1 & dety<=detysize & detz>=1 & detz<=detzsize & ...
                ~(dety>=BeamStopY(1) & dety<=BeamStopY(2) & ...
                    detz>=BeamStopZ(1) & detz<=BeamStopZ(2)));
             if ~isempty(select3)
                for j=1:length(select3)          
                    SubA{grainno}(nr,1) = nr;
                    SubA{grainno}(nr,2) = grainno;
                    SubA{grainno}(nr,3) = nrefl;
                    SubA{grainno}(nr,4:6) = hkl(:,select3(j))';
                    SubA{grainno}(nr,7) = Ahkl(select3(j),5);
                    SubA{grainno}(nr,8) = phi1*180/pi;
                    SubA{grainno}(nr,9) = Phi*180/pi;
                    SubA{grainno}(nr,10) = phi2*180/pi;
                    SubA{grainno}(nr,11:13) = Gt(:,select3(j))';
                    SubA{grainno}(nr,14) = rot; % omega [deg]
                    SubA{grainno}(nr,15)= 2*theta(select3(j))*180/pi;
                    eta=acos(dot([0 Gt(2,select3(j)) Gt(3,select3(j))]/norm([0 Gt(2,select3(j)) Gt(3,select3(j))]), ...
                        [0 0 1]/norm([0 0 1])));
                    SubA{grainno}(nr,16) = eta*180/pi;% [0 360] eta [deg]

                    SubA{grainno}(nr,17) = dety(select3(j));
                    SubA{grainno}(nr,18) = detz(select3(j));
               
                    Lorentz=1./(sin(2*theta(select3(j))));
                    SubA{grainno}(nr,19)=Lorentz;
                    P=(1+costth(select3(j))^2)/2;
                    SubA{grainno}(nr,20)=P;

                    %Diffracted intensity
                    SubA{grainno}(nr,21)=0;
                    ee=min(find(Energy>(Energy_hkl(select3(j))-(Energy(2)-Energy(1))) & Energy<(Energy_hkl(select3(j))+(Energy(2)-Energy(1)))));
                    [A_Ehkl, L_total]=beam_attenuation(SamposW,Lsam2sou,Lsam2det,dety(select3(j)),detz(select3(j)), ...
                    atomparam,Transmission,rou,Energy_hkl(select3(j)),Rsample);       % attenuation intensity factor, June 22, 2020
                    [DQE_Ehkl]=Detector_efficiency(CsI,Swank,Energy_hkl(select3(j))); % DQE, June 22, 2020                        
                    if SubGrain{grainno}(subgrainno,6)==Inf % few cases the grain volume is Inf due to meshing
                        SubGrain{grainno}(subgrainno,6)=mean(setdiff(SubGrain{grainno}(:,6),Inf,'rows'));
                    end
                    if SubGrain{grainno}(subgrainno,6)>1 && SubGrain{grainno}(subgrainno,6)<400 % identify unit as um
                        K2(ee) = lambda(ee)^3*SubGrain{grainno}(subgrainno,5)*10^12/V^2; % [dimensionless]
                    else
                        K2(ee) = lambda(ee)^3*SubGrain{grainno}(subgrainno,5)*10^21/V^2; % [dimensionless] % identify unit as mm
                    end
                    K2(ee) = A_Ehkl*DQE_Ehkl*K2(ee); % consider attenuation and detector efficiency
                    SubA{grainno}(nr,21) = SubA{grainno}(nr,21)+K1*K2(ee)*abs(I0E(ee))*Lorentz*P*Ahkl(select3(j),5)*ExpTime; % intensity [photons]
                    SubA{grainno}(nr,22) = Energy_hkl(select3(j));
                    SubA{grainno}(nr,23) = subgrainno;
                    nr=nr+1;
                    nrefl=nrefl+1;
                end
            end
        end % Loop over subgrains

        SubA_eff{grainno}=[];
        if ~isempty(SubA{grainno})
            for kk=1:length(SubA{grainno}(:,1))
                if (~(all(SubA{grainno}(kk,:))==0) || SubA{grainno}(kk,21)>0) ...
                        && (SubA{grainno}(kk,17)>=1 && SubA{grainno}(kk,17)<detysize ...
                        && SubA{grainno}(kk,18)>=1 && SubA{grainno}(kk,18)<detzsize)
                    SubA_eff{grainno}=[SubA_eff{grainno};SubA{grainno}(kk,:)]; % select the data contributing to the intensity on the detector
                end
            end
        end
        A=[A;SubA_eff{grainno}];
        grainno;
    end % loop over grains
    
    if ~isempty(A)
        %Make diffraction images
        [SpotsPair,GrainIndex,GrainIndex_unique]=make_projection(A,imdata,Spots,SubGrain,hkl_square,hkl_color, ...
                    Lsam2sou,Lsam2det,pixelysize,pixelzsize,detysize,detzsize,direc,prefix,rot_number,rot,grain_centroid);
        SpotsPair_all=[SpotsPair_all;SpotsPair];
        GrainIndex_all=[GrainIndex_all;GrainIndex];

        A_rot{rot_number}=A;
        header_A_rot = ['ReflectionNo.' ' ' 'GrainNo.' ' ' 'NumberOfReflection' ' ' 'h' ' ' 'k' ' ' 'l' ' ' 'F^2' ' ' ...
            'phi1' ' ' 'Phi' ' ' 'phi2' ' ' 'Gw(1)' ' ' 'Gw(2)' ' ' 'Gw(3)' ' ' 'Omega' ' ' '2-theta' ' ' ...
            'eta' ' ' 'det_y' ' ' 'det_z' ' ' 'LorentzFactor' ' ' 'PolarizationFactor' ' ' 'IntegratedIntensity' ...
            ' ' 'EnergyHKL' ' ' 'Subgrainno'];
        if ~isempty(A_rot{rot_number})
            fid=fopen(fullfile(direc2,strcat(num2str(rot_number-1),'A_rot.txt')),'wt');
            fprintf(fid, [header_A_rot '\n']);
            for pp=1:length(A_rot{rot_number}(:,1))
                fprintf(fid, '%d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d\n', A_rot{rot_number}(pp,:));
            end
            fclose(fid);
        end

        header_SpotsPair = ['GrainNo' ' ' 'posX' ' ' 'posY' ' ' 'posZ' ' ' 'h' ' ' 'k' ' ' 'l' ' ' ...
            'phi1' ' ' 'Phi' ' ' 'phi2' ' ' 'Energy_hkl' ' ' 'rot' ' ' 'dety' ' ' 'detz' ' ' ...
            'spotpair_exp_ID' ' ' 'dety_exp' ' ' 'detz_exp' ' ' 'disY' ' ' 'disZ' ' ' 'dis_euclidian' ...
            ' ' 'spotsize' ' ' 'spotsize_exp'];
        if ~isempty(SpotsPair)
            fid=fopen(fullfile(direc2,strcat(num2str(rot_number-1),'SpotsPair.txt')),'wt');
            fprintf(fid, [header_SpotsPair '\n']);
            for pp=1:length(SpotsPair(:,1))
                fprintf(fid, '%d %f %f %f %d %d %d %f %f %f %f %d %f %f %d %f %f %f %f %f %f %f\n', SpotsPair(pp,:));
            end
            fclose(fid);
        end

        if exist('GrainIndex','var')
            A_GrainIndex{rot_number}=GrainIndex;
            header_A_GrainIndex = ['SpotID' ' ' 'SpotSize' ' ' 'GrainID' ' ' 'h' ' ' 'k' ' ' 'l' ' ' ...
                'AverageIntensiy' ' ' 'ReflectionNo' ' ' 'OverlapFraction' ' ' 'IntInt' ' ' 'rot'];
            if ~isempty(A_GrainIndex{rot_number})
                fid=fopen(fullfile(direc2,strcat(num2str(rot_number-1),'GrainIndex.txt')),'wt');
                fprintf(fid, [header_A_GrainIndex '\n']);
                for pp=1:length(A_GrainIndex{rot_number}(:,1))
                    fprintf(fid, '%d %d %d %d %d %d %f %d %d %f %d\n', A_GrainIndex{rot_number}(pp,:));
                end
                fclose(fid);
            end
            A_GrainIndex_unique{rot_number}=GrainIndex_unique;
            if ~isempty(A_GrainIndex_unique{rot_number})
                fid=fopen(fullfile(direc2,strcat(num2str(rot_number-1),'GrainIndex_unique.txt')),'wt');
                fprintf(fid, [header_A_GrainIndex '\n']);
                for pp=1:length(A_GrainIndex_unique{rot_number}(:,1))
                    fprintf(fid, '%d %d %d %d %d %d %f %d %d %f %d\n', A_GrainIndex_unique{rot_number}(pp,:));
                end
                fclose(fid);
            end
        end
    end
    rot
end % loop over rotations

% save graininfo
header_graininfo = ['GrainNo.' ' ' 'GrainDiameter' ' ' 'GrainVolume' ' ' 'EulerAngle(phi1)' ' ' 'EulerAngle(Phi)' ...
    ' ' 'EulerAngle(phi2)' ' ' 'U11' ' ' 'U12' ' ' 'U13' ' ' 'U21' ' ' 'U22' ' ' 'U23' ' ' ...
    'U31' ' ' 'U32' ' ' 'U33' ' ' 'SubGrainNo'];
if ~isempty(graininfo)
    fid=fopen(fullfile(direc2,'graininfo.txt'),'wt');
    fprintf(fid, [header_graininfo '\n']);
    for i=1:length(graininfo(:,1))
        fprintf(fid, '%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d\n', graininfo(i,:));
    end
    fclose(fid);
end
header_SpotsPair = ['GrainNo' ' ' 'posX' ' ' 'posY' ' ' 'posZ' ' ' 'h' ' ' 'k' ' ' 'l' ' ' ...
    'phi1' ' ' 'Phi' ' ' 'phi2' ' ' 'Energy_hkl' ' ' 'rot' ' ' 'dety' ' ' 'detz' ' ' ...
    'spotpair_exp_ID' ' ' 'dety_exp' ' ' 'detz_exp' ' ' 'disY' ' ' 'disZ' ' ' 'dis_euclidian' ...
    ' ' 'spotsize' ' ' 'spotsize_exp'];
if ~isempty(SpotsPair_all)
    fid=fopen(fullfile(direc2,'SpotsPair_all.txt'),'wt');
    fprintf(fid, [header_SpotsPair '\n']);
    for pp=1:length(SpotsPair_all(:,1))
        fprintf(fid, '%d %f %f %f %d %d %d %f %f %f %f %d %f %f %d %f %f %f %f %f %f %f\n', SpotsPair_all(pp,:));
    end
    fclose(fid);
end
header_A_GrainIndex = ['SpotID' ' ' 'SpotSize' ' ' 'GrainID' ' ' 'h' ' ' 'k' ' ' 'l' ' ' 'AverageIntensiy' ...
    ' ' 'ReflectionNo' ' ' 'OverlapFraction' ' ' 'IntInt' ' ' 'rot'];
if ~isempty(GrainIndex_all)
    fid=fopen(fullfile(direc2,'GrainIndex_all.txt'),'wt');
    fprintf(fid, [header_A_GrainIndex '\n']);
    for pp=1:length(GrainIndex_all(:,1))
        fprintf(fid, '%d %d %d %d %d %d %f %d %d %f %d\n', GrainIndex_all(pp,:));
    end
    fclose(fid);
end

% number of spots
SpotNr_simu=zeros(1,abs(grains));
SpotNr_obs=zeros(1,abs(grains));
for grainno=simu_grainno%1:abs(grains)
    if ~isempty(GrainIndex_all)
        SpotNr_simu(grainno)=length(find(GrainIndex_all(:,3)==grainno));
    end
    if ~isempty(SpotsPair_all)
        SpotNr_obs(grainno)=length(find(SpotsPair_all(:,1)==grainno));
    end
end

% count the ideal spot number for each hkl
Spots_all=[];
for i=1:length(A_rot)
    if ~isempty(A_rot{i})
        [A_rot_unique,ia,ic] = unique(A_rot{i}(:,[2 4 5 6]),'rows');
        Spots_all=[Spots_all;A_rot{i}(ia,:)];
    end
    i;
end
hklReflection_all=Spots_all(:,4:6);
hklReflection_all=abs(hklReflection_all);
hklReflection_all=sort(hklReflection_all,2);
[hklReflection,ia,ic]=unique(hklReflection_all,'rows');
hklReflection(:,4)=hklReflection(:,1).^2+hklReflection(:,2).^2+hklReflection(:,3).^2;
hklReflection = sortrows(hklReflection,4);
hklReflection=hklReflection(:,1:3);
for grainno = simu_grainno
    i=grainno;
    SpotsID{i}=Spots_all(find(Spots_all(:,2)==i),:);
    for j=1:length(hklReflection(:,1))
        SpotNr_simu_hkl(i,j)=0;
        for k=1:length(SpotsID{i}(:,1))
            if (sum(sort(abs(SpotsID{i}(k,4:6)),2)==hklReflection(j,:))==3 && ...
                ~(SpotsID{i}(k,17)>BeamStopY(1) && SpotsID{i}(k,17)<BeamStopY(2) && ...
                SpotsID{i}(k,18)>BeamStopZ(1) && SpotsID{i}(k,18)<BeamStopZ(2)))
               SpotNr_simu_hkl(i,j)=SpotNr_simu_hkl(i,j)+1;
            end
        end
        SpotNr_obs_hkl(i,j)=0;
        if ~isempty(SpotsPair_all)
            for k=1:length(SpotsPair_all(:,1))
                if sum(sort(abs(SpotsPair_all(k,5:7)),2)==hklReflection(j,:))==3 && i==SpotsPair_all(k,1)
                    SpotNr_obs_hkl(i,j)=SpotNr_obs_hkl(i,j)+1;
                end
            end
        end
    end
end

