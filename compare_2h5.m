clear all;
close all;

%%% change the default settings
setMTEXpref('xAxisDirection','east'); % default: 'north'
setMTEXpref('zAxisDirection','outOfPlane'); % same as default
setMTEXpref('bAxisDirection','north'); % default: 'east'
setMTEXpref('aAxisDirection',''); % same as default
setMTEXpref('FontSize',44); % default: 15
% define symmetries
cs = crystalSymmetry('cubic');
ss = specimenSymmetry('orthorhombic');
cK = ipfHSVKey(cs);
cK.inversePoleFigureDirection=zvector;
    
% FileFolder='D:\Documents\Matlab\LabDCT_simap\SR_DCT_AlCu6wt_middle_thinned_20210915';
% h5FileName='SR_DCT_AlCu6wt_middle_thinned_trans2.h5';
FileFolder='D:\Documents\Matlab\LabDCT_simap\SR_DCT_AlCu6wt_middle_thinned_20210915_wolfgang';
h5FileName='SR_DCT_AlCu6wt_middle_thinned.h5';
Prefix_Name=h5FileName(1:end-3);
DS_ref=readLabDCT(fullfile(FileFolder,h5FileName));
euler = Euler(rodrigues2quat(vector3d(double(DS_ref.RodVec'))));
euler_grains=euler.*180/pi;% Euler angles [deg]: ([0 2*pi], [0 pi], [0 2*pi])
grainvolume=DS_ref.nVox*DS_ref.VoxSize(1)*DS_ref.VoxSize(2)*DS_ref.VoxSize(3)*1e9; % [um^3]
grainsize=2*(3*grainvolume/(4*pi)).^(1/3); % equivalent diameter of grain [um]
BoxDim=DS_ref.Dimension.*DS_ref.VoxSize'*1000; % [um]
grains=length(nonzeros(grainsize));
rot_ref = rotation('Euler',euler_grains(:,1)*degree,euler_grains(:,2)*degree,euler_grains(:,3)*degree);
o_ref=orientation(rot_ref,cs,ss);
GrainColor=IPFplot(o_ref,cs,ss,0); % plot IPF001

% FileFolder_rec='D:\Documents\Matlab\GrainRecon\AlCu_8wt_middle_thinned_0930_rec\subvol_27_assemble_drop0p02_OR2deg_angle_1_minComp_0p4_2nd';
FileFolder_rec='D:\Documents\Matlab\GrainRecon\AlCu_8wt_middle_thinned_0930_rec\fullvol_v3';
h5FileName='AlCu_8wt_middle_thinned_0930.h5';
Prefix_Name=h5FileName(1:end-3);
DS_rec=readLabDCT(fullfile(FileFolder_rec,h5FileName));

euler = Euler(rodrigues2quat(vector3d(double(DS_rec.RodVec'))));
euler_grains_rec=euler.*180/pi;% Euler angles [deg]: ([0 2*pi], [0 pi], [0 2*pi])
grainvolume_rec=DS_rec.nVox*DS_rec.VoxSize(1)*DS_rec.VoxSize(2)*DS_rec.VoxSize(3)*1e9; % [um^3]
grainsize_rec=2*(3*grainvolume_rec/(4*pi)).^(1/3); % equivalent diameter of grain [um]
BoxDim_rec=DS_rec.Dimension.*DS_rec.VoxSize'*1000; % [um]
grains_rec=length(nonzeros(grainsize_rec));
rot_rec = rotation('Euler',euler_grains_rec(:,1)*degree,euler_grains_rec(:,2)*degree,euler_grains_rec(:,3)*degree);
o_rec = orientation(rot_rec,cs,ss);
GrainColor_rec=IPFplot(o_rec,cs,ss,0); % plot IPF001

% Orientations in IPF map
figure('Name','Overlay orientations');
plotIPDF(o_ref,GrainColor,vector3d(0,0,1),'points','all','MarkerSize',6,'resolution',0.5*degree);
hold on;
plotIPDF(o_rec,GrainColor_rec,vector3d(0,0,1),'Marker','o','all','MarkerFaceColor','none', ...
    'LineWidth',1.5,'MarkerSize',10,'resolution',0.5*degree);
hold off;

% find grain pairs
% Su_angle=[226.6219    0.0002  132.6778]; % Nov 10
Su_angle=[0 0 0];
Su=euler2u(Su_angle(1)*pi/180,Su_angle(2)*pi/180,Su_angle(3)*pi/180);
count=0;
PairIndex=[];
clear ORdiff ORdiff1;
for i=1:length(GrainColor(:,1))
    if grainsize(i)>0
        for j=1:length(GrainColor_rec(:,1))
            grainsize_diff=abs(grainsize(i)-grainsize_rec(j))/grainsize(i);
            color_diff=abs(GrainColor(i,:)-GrainColor_rec(j,:));
            E1=euler_grains(i,:);
            E2=euler_grains_rec(j,:);
            U1=euler2u(E1(1)*pi/180,E1(2)*pi/180,E1(3)*pi/180);
            if exist('Su','var')
                U1=Su*U1;
            end
            U2=euler2u(E2(1)*pi/180,E2(2)*pi/180,E2(3)*pi/180);
            [ang,ax]=misori2(U1,U2);
			ang1=angle(o_ref(i),o_rec(j))/degree;
            if ang1<4
                count=count+1;
                PairIndex(count,:)=[i j DS_ref.EulerZXZ(i,:) rod2euler(DS_rec.RodVec(j,:)) grainsize(i) grainsize_rec(j) ang1];
                ORdiff(count)=ang1;% consider crystal and sample symmetry
                ORdiff1(count)=ang;
                Gsizediff(count)=grainsize_diff;
            end
        end
    end
    i
end
[length(PairIndex(:,1)) length(unique(PairIndex(:,1))) length(unique(PairIndex(:,2)))]
% slice_ref=slice_show(DS_ref,round(DS_ref.Dimension(2)/2),0,1,0)
% slice_rec=slice_show(DS_rec,round(DS_rec.Dimension(2)/2),0,1,0)
problem_ID=setdiff([1:length(DS_rec.SeedID)],unique(PairIndex(:,2)))



% % optimize the rotation matrix by fitting
select_ind=setdiff(1:length(PairIndex(:,1)),[16    22    23    24    25    26    59]);
PairIndex_select=PairIndex(select_ind,:);



% % optimize the rotation matrix by fitting
PairIndex_select=PairIndex(PairIndex(:,10)>20,:);
[~,ind]=unique(PairIndex_select(:,1));
PairIndex_select=PairIndex_select(ind,:);

x=[0 0 0];
Su=euler2u(x(1)*pi/180,x(2)*pi/180,x(3)*pi/180);
Err=[];
for j=1:length(PairIndex_select(:,1))
    u0=euler2u(PairIndex_select(j,3)*pi/180,PairIndex_select(j,4)*pi/180,PairIndex_select(j,5)*pi/180);
    u1=Su*u0;
    u_ref=euler2u(PairIndex_select(j,6)*pi/180,PairIndex_select(j,7)*pi/180,PairIndex_select(j,8)*pi/180);
    [ang,~]=misori2(u1,u_ref);
%         Err(j)=sum(sum((u1-u_ref).^2));
%         Err(j)=sqrt(Err(j));
    Err(j)=ang;
end
%%%% further remove unqualified pairs
ind=find(Err<4);
PairIndex_select=PairIndex_select(ind,:);
% then loop initial values from the whole space and select the optimal one
    
% % fit U
cs = crystalSymmetry('cubic');
clear x;
Err2Sum=@(x)FitOR(x,PairIndex_select(:,3:5),PairIndex_select(:,6:8),cs);
% x0=[226.6219 0.0002 132.6778]; % Nov 10
% x0=[216.6219   -0.2161  142.6201];
% x0=[226.2761    -2  130.9372];
x0=[332.9059    0.4280   26.3354]; % initial value obtained by looping the whole OR space
% x0=[339.3184    0.1608   20.4896];
% optimized solution [337.6452    0.4312 21.5984] => < ang > = 0.059 deg
LB=[max([x0(1)-5 0]) max([x0(2)-5 -5]) max([x0(3)-5 0])];
UB=[x0(1)+5 x0(2)+5 x0(3)+5];
opts=optimset('Display','iter','Algorithm',{'levenberg-marquardt',0.005}, ...
    'MaxFunEvals',2000,'MaxIter',2000,'TolFun',1e-10, ...
    'PlotFcns',@optimplotfval);        
[x,fval,exitflag,output] = fminsearchbnd(Err2Sum,x0,LB,UB,opts);
sprintf('The resulted residual is: %0.5f',fval)




ang=[];
% Su_angle=[226.6219 0.0002 132.6778]; % Nov 10
% Su_angle=[337.6452 0.4312 21.5984]; % Nov 30, fitting to minimize the misOR between SR-DCT and LabDCT
Su_angle=[93.4422   -0.5489  263.0114];% Jan 10, new recon by wolfgang
Su=euler2u(Su_angle(1)*pi/180,Su_angle(2)*pi/180,Su_angle(3)*pi/180);
for j=1:length(PairIndex_select(:,1))
    U1=euler2u(PairIndex_select(j,3)*pi/180,PairIndex_select(j,4)*pi/180,PairIndex_select(j,5)*pi/180);
    U1=Su*U1;
    U2=euler2u(PairIndex_select(j,6)*pi/180,PairIndex_select(j,7)*pi/180,PairIndex_select(j,8)*pi/180);
    [ang(j),~]=misori2(U1,U2);
end

% assign initial values by looping the whole OR space
input_al;
OR_folder='./ori_set_hyperspherical';
OR=get_ori_set(OR_folder,sgno);
tic;
cs = crystalSymmetry('cubic');
clear x;
Err2Sum=@(x)FitOR(x,PairIndex_select(:,3:5),PairIndex_select(:,6:8),cs);
for j=1:OR.len
    x0=OR.euler_angles(j,:);
    LB=[max([x0(1)-2 0]) max([x0(2)-2 -3]) max([x0(3)-2 0])];
    UB=[x0(1)+2 x0(2)+2 x0(3)+2];
    opts=optimset('Display','off','Algorithm',{'levenberg-marquardt',0.005}, ...
        'MaxFunEvals',2000,'MaxIter',2000,'TolFun',1e-10);
    [x,fval(j),exitflag,output] = fminsearchbnd(Err2Sum,x0,LB,UB,opts);
    j
end
toc
[value,ind]=min(fval);
x0=OR.euler_angles(ind,:);

% % fit U for SR_DCT_AlCu6wt_middle_thinned_20210915_wolfgang
cs = crystalSymmetry('cubic');
clear x;
Err2Sum=@(x)FitOR(x,PairIndex_select(:,3:5),PairIndex_select(:,6:8),cs);
% initial value obtained by looping the whole OR space, Jan 10, 2021
x0=[96.9088    1.4224  261.9760];
% x0=[265.7139    1.1427   88.0202];
% optimized solution [93.4421   -0.5489  263.0115] => < ang > = 0.0509 deg
LB=[max([x0(1)-5 0]) max([x0(2)-5 -5]) max([x0(3)-5 0])];
UB=[x0(1)+5 x0(2)+5 x0(3)+5];
opts=optimset('Display','iter','Algorithm',{'levenberg-marquardt',0.005}, ...
    'MaxFunEvals',2000,'MaxIter',2000,'TolFun',1e-10, ...
    'PlotFcns',@optimplotfval);        
[x,fval,exitflag,output] = fminsearchbnd(Err2Sum,x0,LB,UB,opts);
sprintf('The resulted residual is: %0.5f',fval)


% verification
ang=[];
Su_angle=x;% Jan 10, new recon by wolfgang
Su=euler2u(Su_angle(1)*pi/180,Su_angle(2)*pi/180,Su_angle(3)*pi/180);
for j=1:length(PairIndex_select(:,1))
    U1=euler2u(PairIndex_select(j,3)*pi/180,PairIndex_select(j,4)*pi/180,PairIndex_select(j,5)*pi/180);
    U1=Su*U1;
    U2=euler2u(PairIndex_select(j,6)*pi/180,PairIndex_select(j,7)*pi/180,PairIndex_select(j,8)*pi/180);
    [ang(j),~]=misori2(U1,U2);
end
[length(find(ang<0.1)) mean(ang(ang<0.1))]
