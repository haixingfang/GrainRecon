% get geometry information
% % input: geometry of DCT images from Xact;
%          xml file from Xact tomo reconstruction
% August 16, 2021

clear all;
close all;
% geoFolder='D:\ExpData_process\CCD_2021_07_22_AlCu8wt\DCT_1_6_all';
% geoFolder='D:\ExpData\2022_01_12_DCT\Pure_Fe_D120_micro_source_DCT\dct_tomo\3_DCT';
% geoFolder='D:\ExpData_process\CCD_2022_01_19_AlCu_8wt_middle_thinned_micro_source\3_DCT';
% geoFolder='D:\ExpData_process\CCD_2022_01_20_AlCu_8wt_pin_580C_1h\3_DCT';
% geoFolder='D:\ExpData_process\CCD_2022_03_10_AlCu_8wt_pin_580C_1h_nanoS\3_DCT';
% geoFolder='D:\ExpData_process\Varian_2022_07_20_AlCu_8wt_580C_1h_nanoS\3_DCT';
% geoFolder='E:\ExpData_process\Varian_2022_07_20_AlCu_8wt_580C_1h_nanoS\3_DCT';
% geoFolder='E:\ExpData_process\Pixirad_2022_07_20_AlCu_8wt_580C_1h_nanoS\dct_scan_10_50keV';
% geoFolder='E:\2022_07_20_DCT\AlCu_8wt_580C_1h_pixirad\tomo_recon\Proj';
% geoFolder='E:\2022_07_20_DCT\Al2024_T26_350C_15min_flat_panel\auto_scan\dct_tomo\3_DCT';
% geoFolder='E:\ExpData_process\CCD_2022_09_14_AlCu_8wt_pin_580C_1h_nanoS\3_DCT';
geoFolder='E:\ExpData_process\CCD_2022_09_23_AlCu_8wt_pin_580C_1h_nanoS_after_calib\3_DCT';
% geoFolder='E:\2022_09_23_test\AlCu_8wt_580C_1h_varian\dct_tomo\dct_tomo\3_DCT';
geoName='geometry_all.csv';
geoTable=readtable(fullfile(geoFolder,geoName),'NumHeaderLines',2);
P0=[];P1=[];P2=[];P3=[];
nframe=6;
omega_step=3;
j=0;
for i=1:nframe:length(geoTable.Var1)
    j=j+1;
    omega=(i-1)/nframe*omega_step+0;
    Omega=[cosd(omega) -sind(omega) 0;sind(omega) cosd(omega) 0;0 0 1];
    pos=Omega*[-geoTable.Var4(i) -geoTable.Var2(i) geoTable.Var3(i)]';
    P0=[P0;pos'];    
    pos=Omega*[-geoTable.Var7(i) -geoTable.Var5(i) geoTable.Var6(i)]';
    P1=[P1;pos'];
    pos=Omega*[-geoTable.Var10(i) -geoTable.Var8(i) geoTable.Var9(i)]';
    P2=[P2;pos'];
    pos=Omega*[-geoTable.Var13(i) -geoTable.Var11(i) geoTable.Var12(i)]';
    P3=[P3;pos'];
    
    P21=P2(end,:)-P1(end,:);
    P32=P3(end,:)-P2(end,:);
    P31=P3(end,:)-P1(end,:);
    ang(j)=acosd(dot(P21/norm(P21),P31/norm(P31)));
    det_normal(j,:)=cross(P21,P31)./norm(cross(P21,P31));
    
%     tilt(j,1)=acosd(dot(det_normal(j,:)/norm(det_normal(j,:)),[1 0 0])); % [deg]
%     tilt(j,2)=acosd(dot(det_normal(j,:)/norm(det_normal(j,:)),[0 1 0]))-90; % [deg]
%     tilt(j,3)=acosd(dot(det_normal(j,:)/norm(det_normal(j,:)),[0 0 1]))-90; % [deg]
    
    % using the method as GrainRecon
    detdiru(j,:)=P31/norm(P31);
    detdirv(j,:)=-P21/norm(-P21);
    tilt(j,1) = acosd(dot(detdiru(j,:),[0 0 1]))-90;
    tilt(j,2) = acosd(dot(det_normal(j,:)/norm(det_normal(j,:)).*[1 0 1],[0 0 1]))-90; % equivalent to acosd(dot(det_normal.*[1 0 1],[1 0 0]))-90
    tilt(j,3) = 90-acosd(dot(det_normal(j,:)/norm(det_normal(j,:)).*[1 1 0],[0 1 0])); % equivalent to acosd(dot(det_normal.*[1 1 0],[1 0 0]))
end

figure;
hold all;
plot([1:length(tilt(:,1))],tilt(:,1),'ro');
plot([1:length(tilt(:,1))],tilt(:,2),'bx');
plot([1:length(tilt(:,1))],tilt(:,3),'m^');
box on;
legend('\psi_x','\psi_y','\psi_z');
xlabel('Projection No.');
ylabel('Tilt angle (^{o})');
set(gca,'FontSize',18);
set(gca,'LineWidth',1.5);

figure;
subplot(2,3,1);
plot([1:length(P0(:,1))],P0(:,1),'ro');
xlabel('Projection No.');
ylabel('Position (mm)');
set(gca,'FontSize',18);
set(gca,'LineWidth',1.5);
title('(a) P0_{x}');
subplot(2,3,2);
plot([1:length(P0(:,1))],P0(:,2),'bx');
xlabel('Projection No.');
ylabel('Position (mm)');
set(gca,'FontSize',18);
set(gca,'LineWidth',1.5);
title('(b) P0_{y}');
subplot(2,3,3);
plot([1:length(P0(:,1))],P0(:,3),'m^');
xlabel('Projection No.');
ylabel('Position (mm)');
set(gca,'FontSize',18);
set(gca,'LineWidth',1.5);
title('(c) P0_{z}');

subplot(2,3,4);
plot([1:length(P1(:,1))],P1(:,1),'ro');
xlabel('Projection No.');
ylabel('Position (mm)');
set(gca,'FontSize',18);
set(gca,'LineWidth',1.5);
title('(d) P1_{x}');
subplot(2,3,5);
plot([1:length(P1(:,1))],P1(:,2),'bx');
xlabel('Projection No.');
ylabel('Position (mm)');
set(gca,'FontSize',18);
set(gca,'LineWidth',1.5);
title('(e) P1_{y}');
subplot(2,3,6);
plot([1:length(P1(:,1))],P1(:,3),'m^');
xlabel('Projection No.');
ylabel('Position (mm)');
set(gca,'FontSize',18);
set(gca,'LineWidth',1.5);
title('(f) P1_{z}');
% print('output','-dtiff','-r300')
% [std(P0(:,1)) std(P0(:,2)) std(P0(:,3))]
% [std(P1(:,1)) std(P1(:,2)) std(P1(:,3))]


%%%%% acquisition file
% xmlFolder='D:\ExpData\2021_07_22_AlCu8wt\tomo-0723-no-pinhole';
% xmlFolder='D:\ExpData\2022_01_19_DCT\Al8Cu_pin_580C_1h\tomo_no_pinhole';
% xmlFolder = uigetdir('D:\ExpData\');
xmlFolder='E:\2022_09_23_test\AlCu_8wt_580C_1h_CCD_after_calibration\tomo_no_pinhole';
% xmlFolder='E:\2022_09_23_test\AlCu_8wt_580C_1h_varian\tomo_no_pinhole_800projs';
xmlName='unireconstruction.xml';
xmlfile=xml2struct(fullfile(xmlFolder,xmlName));
% from tomo file: Y=>Z, X=>-Y, Z=>-X (left: simap setup; right: conventional system)
aquisition_center(1)=-str2num(xmlfile.unireconstruction.conebeam.volume_acquisition.offset.Attributes.Z);
aquisition_center(2)=-str2num(xmlfile.unireconstruction.conebeam.volume_acquisition.offset.Attributes.X);
aquisition_center(3)=str2num(xmlfile.unireconstruction.conebeam.volume_acquisition.offset.Attributes.Y);
Lsam2sou=str2num(xmlfile.unireconstruction.conebeam.acquisitioninfo.geometry.Attributes.sod); % [mm]
Lsam2det=str2num(xmlfile.unireconstruction.conebeam.acquisitioninfo.geometry.Attributes.sdd)-Lsam2sou; % [mm]
RotAxisOffset=-str2num(xmlfile.unireconstruction.conebeam.correction.Attributes.offsetX); % Rotation axis offset [detector pixel], Aug 17, 2022

%% output the geometry information for reconstruction
% in conventional lab coordinate system
fprintf('Lsam2sou = %0.4f mm\n',Lsam2sou)
fprintf('Lsam2det = %0.4f mm\n',Lsam2det)
fprintf('P0 = %0.4f, %0.4f, %0.4f (mm)\n',mean(P0))
fprintf('acquisition_center = %0.4f, %0.4f, %0.4f (mm)\n',aquisition_center)
fprintf('det_center = %0.4f, %0.4f, %0.4f (mm)\n',mean(P1))
fprintf('detdiru = [%0.6f %0.6f %.6f]; detdirv = [%0.6f %0.6f %.6f]\n', ...
                mean(detdiru)./norm(mean(detdiru)),mean(detdirv)./norm(mean(detdirv)));
fprintf('det_normal = %0.6f, %0.6f, %0.6f\n',mean(det_normal))
fprintf('RotAxisOffset = %0.2f pixel (detector)\n',RotAxisOffset)
fprintf('tilt_angles = %0.4f, %0.4f, %0.4f (deg)\n',mean(tilt))


