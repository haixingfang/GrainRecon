% get geometry information
% % input: geometry of DCT images from Xact;
%          xml file from Xact tomo reconstruction
% August 16, 2021

clear all;
close all;
% geoFolder='D:\ExpData_process\CCD_2021_07_22_AlCu8wt\DCT_1_6_all';
% geoFolder='D:\ExpData\2022_01_12_DCT\Pure_Fe_D120_micro_source_DCT\dct_tomo\3_DCT';
% geoFolder='D:\ExpData_process\CCD_2022_01_19_AlCu_8wt_middle_thinned_micro_source\3_DCT';
geoFolder='D:\ExpData_process\CCD_2022_01_20_AlCu_8wt_pin_580C_1h\3_DCT';
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
    tilt(j,1)=acosd(dot(det_normal(j,:)/norm(det_normal(j,:)),[1 0 0])); % [deg]
    tilt(j,2)=acosd(dot(det_normal(j,:)/norm(det_normal(j,:)),[0 1 0]))-90; % [deg]
    tilt(j,3)=acosd(dot(det_normal(j,:)/norm(det_normal(j,:)),[0 0 1]))-90; % [deg]
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
xmlFolder = uigetdir('D:\ExpData\');
xmlName='unireconstruction.xml';
xmlfile=xml2struct(fullfile(xmlFolder,xmlName));
% from tomo file: Y=>Z, X=>-Y, Z=>-X (left: simap setup; right: conventional system)
aquisition_center(1)=-str2num(xmlfile.unireconstruction.conebeam.volume_acquisition.offset.Attributes.Z);
aquisition_center(2)=-str2num(xmlfile.unireconstruction.conebeam.volume_acquisition.offset.Attributes.X);
aquisition_center(3)=str2num(xmlfile.unireconstruction.conebeam.volume_acquisition.offset.Attributes.Y);
Lsam2sou=str2num(xmlfile.unireconstruction.conebeam.acquisitioninfo.geometry.Attributes.sod); % [mm]
Lsam2det=str2num(xmlfile.unireconstruction.conebeam.acquisitioninfo.geometry.Attributes.sdd)-Lsam2sou; % [mm]

%% output the geometry information for reconstruction
sprintf('Lsam2sou = %0.4f mm',Lsam2sou)
sprintf('Lsam2det = %0.4f mm',Lsam2det)
sprintf('P0 = %0.4f, %0.4f, %0.4f (mm)',mean(P0))
sprintf('acquisition_center = %0.4f, %0.4f, %0.4f (mm)',aquisition_center)
sprintf('det_center = %0.4f, %0.4f, %0.4f (mm)',mean(P1))
sprintf('det_normal = %0.4f, %0.4f, %0.4f',mean(det_normal))
sprintf('tilt_angles = %0.4f, %0.4f, %0.4f (deg)',mean(tilt))


