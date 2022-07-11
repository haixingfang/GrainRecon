% compute the distribution of completeness values and intensity correlation
% within each voxel in a specific grain
% input: ground truth simulated data, simu_Fe
% Nov 30, 2021

% read the 16-bit images
ImageFolder='D:\Documents\Matlab\LabDCT_simap\virtual_Fe_100um_11_11_simu\TFT_cmp_simap_s2_geo1_16bit';
for j=1:121
    proj16{j}=imread(fullfile(ImageFolder,['proj' num2str(j-1,'%.4d') '_grey.tif']));
end

% % prepare for calculating intensities
% estimate sample diameter by assuming the sample has approximately a cylinder shape
Rsample=sqrt(SampleVolumeDim(1)^2+SampleVolumeDim(2)^2)/2; % [mm]
if ~exist('Rsample','var') 
    Rsample=0.4; %[mm]
end
if ~exist('ExpTime','var')
    ExpTime=60*6;
end
% read transmission data and CsI scintillator data, add on June 22, 2020
atomparam_atomno=atomparam.atomno;
[Transmission, rou]=ReadTransData(atomparam_atomno); % [(-), g/cm^3]
[CsI, Swank]=ReadCsI();

%%% get input structure
FileFolder='D:\Documents\Matlab\GrainRecon\Fe_100um_11_11_simu\input_structure';
h5FileName='Grain100um_400_400_600_input.h5';
Prefix_Name=h5FileName(1:end-3);
DS=readLabDCT(fullfile(FileFolder,h5FileName));

% add on Nov 25, 2021, to calculate the grain centroid based on weighted completeness
for k=1:length(DS.SeedID)
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
for jj=1:length(DS.SeedID)
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
            SubGrain{jj}(k,5)=DS.VoxSize(1)*DS.VoxSize(2)*DS.VoxSize(3); % volume [m^3]
            SubGrain{jj}(k,6)=2*(3*SubGrain{jj}(k,5)/(4*pi))^(1/3); % EqDiameter
            SubGrain{jj}(k,7)=2*(3*DS.nVox(jj)*DS.VoxSize(1)*DS.VoxSize(2)*DS.VoxSize(3)/(4*pi))^(1/3); % EqDiameter of the parent grain
        end
        SmallGrID=[SmallGrID;jj];        
    end
    if ~isempty(SubGrain{jj})
        if simap_data_flag==1
            SubGrain{jj}(:,2)=-SubGrain{jj}(:,2); % [mm]
            SubGrain{jj}(:,3)=-SubGrain{jj}(:,3); % [mm]
        end
    end
    jj;
end
grains=length(DS.SeedID); % number of grains
grainvolume=DS.nVox*DS.VoxSize(1)*DS.VoxSize(2)*DS.VoxSize(3)*1e9; % [um^3]
grainsize=2*(3*grainvolume/(4*pi)).^(1/3); % equivalent diameter of grain [um]
euler_grains=DS.EulerZXZ;

for grainno = 1%:5%abs(grains)
    yy=DS.GIDvol==DS.SeedID(grainno);
    yy=uint8(yy);
    if DS.nVox(grainno)<4000
        Scaling=1;
    elseif DS.nVox(grainno)<8000
        Scaling=0.5;
    else
        Scaling=0.25;
    end
    zz=imresize3(yy,Scaling,'nearest'); % function valid after 2017a
    Binning=size(yy)./size(zz);
    clear id;
    [id(:,1), id(:,2), id(:,3)]=ind2sub(size(zz),find(zz==1));
    crop_dim{grainno}=[min(id); max(id)];
    gs{grainno}.mask=zeros(size(zz));
    gs{grainno}.Completeness=zeros(size(zz));
    gs{grainno}.dis_median=zeros(size(zz));
    gs{grainno}.Icorr=zeros(size(zz));
    gs{grainno}.Iexp_sum=zeros(size(zz));
    
    gs_idn{grainno}.mask=zeros(size(zz));
    gs_idn{grainno}.Completeness=zeros(size(zz));
    gs_idn{grainno}.dis_median=zeros(size(zz));
    gs_idn{grainno}.Icorr=zeros(size(zz));
    gs_idn{grainno}.Iexp_sum=zeros(size(zz));
    
    phi1 = euler_grains(grainno,1)*pi/180;
    Phi = euler_grains(grainno,2)*pi/180;
    phi2 = euler_grains(grainno,3)*pi/180;
    U = euler2u(phi1,Phi,phi2);
    C_Icorr{grainno}=zeros(length(SubGrain{grainno}(:,1)),9);
        
    id_vol=DS.GIDvol==grainno; % must be one single connected region
    id_dismap=bwdist(id_vol);
    id_neigb_ind=find(id_dismap>0 & id_dismap<3); % distances follow: 1, sqrt(2), sqrt(3), 2, sqrt(5), ...
    id_neigb=unique(DS.GIDvol(id_neigb_ind));
    id_neigb=double(id_neigb(id_neigb>0)); % neighbors ID
    for k=1:length(SubGrain{grainno}(:,1))
        pos=SubGrain{grainno}(k,2:4); % [mm]
        [Nr_simu,Nr_intersect,Icorr,Iexp_sum,dis_median,SimuSpots,HittedSpots]=index_verify_grey(U,proj16,proj_bin_bw,pos,rot_angles,S,B,Ahkl,nrhkl, ...
            RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
            pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
            K1,V,I0E,Energy,ExpTime,atomparam,Transmission,rou,Rsample,CsI,Swank,lambda,VoxSize*mean(Binning));
        C_Icorr{grainno}(k,:)=[Nr_simu Nr_intersect Nr_intersect/Nr_simu dis_median Icorr Iexp_sum euler_grains(grainno,:)];
        gs{grainno}.mask(id(k,1),id(k,2),id(k,3))=1;
        gs{grainno}.Completeness(id(k,1),id(k,2),id(k,3))=Nr_intersect/Nr_simu;
        gs{grainno}.dis_median(id(k,1),id(k,2),id(k,3))=dis_median;
        gs{grainno}.Icorr(id(k,1),id(k,2),id(k,3))=Icorr;
        gs{grainno}.Iexp_sum(id(k,1),id(k,2),id(k,3))=Iexp_sum;
    end
    % suppose the grain voxel extends beyong the true grain boundary
    clear idn idn_pos;
    id_dismap=imresize3(id_dismap,Scaling,'nearest');
    id_neigb_ind=find(id_dismap>0 & id_dismap<2); % distances follow: 1, sqrt(2), sqrt(3), 2, sqrt(5), ...   
    [idn(:,1),idn(:,2),idn(:,3)]=ind2sub(size(zz),id_neigb_ind);
    C_Icorr_idn{grainno}=zeros(length(id_neigb_ind),9);
    for k=1:length(id_neigb_ind)
        idn_pos(k,1)=((idn(k,1)+RecVolumePixel(1,1)/Binning(1)-1)-tomo_scale.Dimension(1)/Binning(1)/2).*DS.VoxSize(1)*Binning(1)+DS.Center(1); % centroid coordinate x
        idn_pos(k,2)=((idn(k,2)+RecVolumePixel(2,1)/Binning(2)-1)-tomo_scale.Dimension(2)/Binning(2)/2).*DS.VoxSize(2)*Binning(2)+DS.Center(2); % centroid coordinate y
        idn_pos(k,3)=((idn(k,3)+RecVolumePixel(3,1)/Binning(3)-1)-tomo_scale.Dimension(3)/Binning(3)/2).*DS.VoxSize(3)*Binning(3)+DS.Center(3); % centroid coordinate z
        if simap_data_flag==1
            idn_pos(k,1)=-idn_pos(k,1); % [mm]
            idn_pos(k,2)=-idn_pos(k,2); % [mm]
        end
        [Nr_simu,Nr_intersect,Icorr,Iexp_sum,dis_median,SimuSpots,HittedSpots]=index_verify_grey(U,proj,proj_bin_bw,idn_pos(k,:),rot_angles,S,B,Ahkl,nrhkl, ...
            RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
            pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
            K1,V,I0E,Energy,ExpTime,atomparam,Transmission,rou,Rsample,CsI,Swank,lambda,VoxSize*mean(Binning));
        C_Icorr_idn{grainno}(k,:)=[Nr_simu Nr_intersect Nr_intersect/Nr_simu dis_median Icorr Iexp_sum euler_grains(grainno,:)];
        gs_idn{grainno}.mask(idn(k,1),idn(k,2),idn(k,3))=1;
        gs_idn{grainno}.Completeness(idn(k,1),idn(k,2),idn(k,3))=Nr_intersect/Nr_simu;
        gs_idn{grainno}.dis_median(idn(k,1),idn(k,2),idn(k,3))=dis_median;
        gs_idn{grainno}.Icorr(idn(k,1),idn(k,2),idn(k,3))=Icorr;
        gs_idn{grainno}.Iexp_sum(idn(k,1),idn(k,2),idn(k,3))=Iexp_sum;
    end
    grainno
end

% plot property as a function of distance to the grain centroid
% prop=gs{grainno}.Completeness;
prop=gs_idn{grainno}.dis_median;
% prop=gs{grainno}.Icorr;
% prop=gs{grainno}.Iexp_sum;
profile=prop_profile(prop,DS,grainno,id);
figure;
errorbar(profile(:,1),profile(:,2),profile(:,3),'ro-')
xlabel('Distance to grain centroid (pixel)');
ylabel('Property value');


% 3D visualization
figure('Name','Grain shape');
plot3D(gs{grainno}.mask,1,[1 0 0],'pasive');
% xlim([crop_dim{grainno}(1,1)-1 crop_dim{grainno}(2,1)+1]);
% ylim([crop_dim{grainno}(1,2)-1 crop_dim{grainno}(2,2)+1]);
% zlim([crop_dim{grainno}(1,3)-1 crop_dim{grainno}(2,3)+1]);
view(3);
set(gca,'visible','off');
print('grain_3D','-dtiff','-r300');

prop=gs_idn{grainno}.Completeness;
prop=gs_idn{grainno}.dis_median;
prop=gs_idn{grainno}.Icorr;
% prop=gs{grainno}.Iexp_sum;
figure('Name','Visualization in 3D');
prop_plot=permute(prop,[2 1 3]);
prop_plot(prop_plot==0)=NaN;
h1 = vol3d('cdata',prop_plot,'texture','3D');
xlim([crop_dim{grainno}(1,1)-1 crop_dim{grainno}(2,1)+1]);
ylim([crop_dim{grainno}(1,2)-1 crop_dim{grainno}(2,2)+1]);
zlim([crop_dim{grainno}(1,3)-1 crop_dim{grainno}(2,3)+1]);
view(3);
% alphamap('vup');
colormap jet;
colorbar('FontSize',16);
set(gca,'visible','on');
xlabel('X (pixel)')
ylabel('Y (pixel)')
zlabel('Z (pixel)')
print('prop_3D','-dtiff','-r300');




[value,ind]=min(C_Icorr{grainno}(:,3))
[value,ind]=max(C_Icorr{grainno}(:,3))
[min(C_Icorr{grainno}(:,3)) max(C_Icorr{grainno}(:,3)) min(C_Icorr{grainno}(:,3))/max(C_Icorr{grainno}(:,3))]
[min(C_Icorr{grainno}(:,4)) max(C_Icorr{grainno}(:,4)) min(C_Icorr{grainno}(:,4))/max(C_Icorr{grainno}(:,4))]
[min(C_Icorr{grainno}(:,5)) max(C_Icorr{grainno}(:,5)) min(C_Icorr{grainno}(:,5))/max(C_Icorr{grainno}(:,5))]

% suppose the voxel has been indexed by an orientation from one of the
% neighboring grain, calculate the completeness
ind=2;
count=0;
for ind=1:length(C_Icorr{grainno}(:,1))
    C_Icorr_false{grainno}=zeros(length(id_neigb),14);
    pos=SubGrain{grainno}(ind,2:4); % [mm]
    phi1 = euler_grains(grainno,1)*pi/180;
    Phi = euler_grains(grainno,2)*pi/180;
    phi2 = euler_grains(grainno,3)*pi/180;
    U = euler2u(phi1,Phi,phi2);
    for j=1:length(id_neigb)
        phi1 = euler_grains(id_neigb(j),1)*pi/180;
        Phi = euler_grains(id_neigb(j),2)*pi/180;
        phi2 = euler_grains(id_neigb(j),3)*pi/180;
        U1 = euler2u(phi1,Phi,phi2);
        [ang(j),~]=misori2(U,U1);
        [Nr_simu,Nr_intersect,Icorr,Iexp_sum,dis_median,SimuSpots,HittedSpots]=index_verify_grey(U1,proj,proj_bin_bw,pos,rot_angles,S,B,Ahkl,nrhkl, ...
            RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
            pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
            K1,V,I0E,Energy,ExpTime,atomparam,Transmission,rou,Rsample,CsI,Swank,lambda,VoxSize);
        C_Icorr_false{grainno}(j,:)=[Nr_simu Nr_intersect Nr_intersect/Nr_simu dis_median Icorr ...
                    Iexp_sum grainno euler_grains(grainno,:) id_neigb(j) euler_grains(id_neigb(j),:)];
    end
    if ~isempty(find(C_Icorr_false{grainno}(:,3)>C_Icorr{grainno}(ind,3))) && ...
           ~isempty(find(C_Icorr_false{grainno}(:,4)<C_Icorr{grainno}(ind,4))) %&& ...
          %  ~isempty(find(C_Icorr_false{grainno}(:,5)>C_Icorr{grainno}(ind,5)))
%     if ~isempty(find(C_Icorr_false{grainno}(:,5)>C_Icorr{grainno}(ind,5)))
        count=count+1;
    end
    ind
end
figure;
subplot(2,2,1);
hold all;
plot(ang,C_Icorr_false{grainno}(:,3),'rx','MarkerSize',14);
plot(0,C_Icorr{grainno}(ind,3),'bo','LineWidth',1.5);
xlabel('misOR (^{o})');
ylabel('Completeness');
set(gca,'LineWidth',1.5);
set(gca,'FontSize',18);
box on;grid on;axis tight;
subplot(2,2,2);
hold all;
plot(ang,C_Icorr_false{grainno}(:,4),'rx','MarkerSize',14);
plot(0,C_Icorr{grainno}(ind,4),'bo','LineWidth',1.5);
xlabel('misOR (^{o})');
ylabel('dis median (pixel)');
set(gca,'LineWidth',1.5);
set(gca,'FontSize',18);
box on;grid on;axis tight;
subplot(2,2,3);
hold all;
plot(ang,C_Icorr_false{grainno}(:,5),'rx','MarkerSize',10);
plot(0,C_Icorr{grainno}(ind,5),'bo','LineWidth',1.5);
xlabel('misOR (^{o})');
ylabel('I_{corr}');
set(gca,'LineWidth',1.5);
set(gca,'FontSize',18);
box on;grid on;axis tight;


[SpotNr_simu_all,SpotNr_obs_all,SpotsPair_all,GrainIndex_all,SubGrain_all,~,SmallGrID]=Forward_simu_spots_exp(DS, ...
    Rsample,RecVolumePixel,tomo_scale,ExpTime,atomparam,proj,Spots,rot_start,rot_step,rot_end, ...
    S,B,Ahkl,nrhkl,hkl_square,Energy,lambda,V,K1,I0E,RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det, ...
    dety00,detz00,P0y,P0z,pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
    simap_data_flag,strcat(OutputFolder,'_true_geo'),[rot_start:2*rot_step:rot_end-180]);

