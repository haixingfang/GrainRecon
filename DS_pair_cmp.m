% pair and compare the grains in two different datasets
% DS_ref: as reference (normally it is ground truth)
% DS_rec: to be compared
% Dec 10, 2021
function [PairIndex,Unpaired]=DS_pair_cmp(DS_ref,DS_rec)

%%% change the default settings
setMTEXpref('xAxisDirection','east'); % default: 'north'
setMTEXpref('zAxisDirection','outOfPlane'); % same as default
setMTEXpref('bAxisDirection','north'); % default: 'east'
setMTEXpref('aAxisDirection',''); % same as default
setMTEXpref('FontSize',44); % default: 15
% define symmetries
cs = crystalSymmetry('cubic');
cK = ipfHSVKey(cs);
cK.inversePoleFigureDirection=zvector;

% DS_ref
euler = Euler(rodrigues2quat(vector3d(double(DS_ref.RodVec'))));
euler_grains=euler.*180/pi;% Euler angles [deg]: ([0 2*pi], [0 pi], [0 2*pi])
grainvolume=DS_ref.nVox*DS_ref.VoxSize(1)*DS_ref.VoxSize(2)*DS_ref.VoxSize(3)*1e9; % [um^3]
grainsize=2*(3*grainvolume/(4*pi)).^(1/3); % equivalent diameter of grain [um]
BoxDim=DS_ref.Dimension.*DS_ref.VoxSize'*1000; % [um]
grains=length(nonzeros(grainsize));
rot_ref = rotation('Euler',euler_grains(:,1)*degree,euler_grains(:,2)*degree,euler_grains(:,3)*degree);
o_ref=orientation(rot_ref,cs);
% GrainColor=IPFplot(o_ref,cs,ss,0); % plot IPF001

% DS_rec
euler = Euler(rodrigues2quat(vector3d(double(DS_rec.RodVec'))));
euler_grains_rec=euler.*180/pi;% Euler angles [deg]: ([0 2*pi], [0 pi], [0 2*pi])
grainvolume_rec=DS_rec.nVox*DS_rec.VoxSize(1)*DS_rec.VoxSize(2)*DS_rec.VoxSize(3)*1e9; % [um^3]
grainsize_rec=2*(3*grainvolume_rec/(4*pi)).^(1/3); % equivalent diameter of grain [um]
BoxDim_rec=DS_rec.Dimension.*DS_rec.VoxSize'*1000; % [um]
grains_rec=length(nonzeros(grainsize_rec));
rot_rec = rotation('Euler',euler_grains_rec(:,1)*degree,euler_grains_rec(:,2)*degree,euler_grains_rec(:,3)*degree);
o_rec = orientation(rot_rec,cs);
% GrainColor_rec=IPFplot(o_rec,cs,ss,0); % plot IPF001

X1=euler_grains;
X2=euler_grains_rec;
rot1 = rotation('Euler',X1(:,1)*degree,X1(:,2)*degree,X1(:,3)*degree);
rot2 = rotation('Euler',X2(:,1)*degree,X2(:,2)*degree,X2(:,3)*degree);
o1=orientation(rot1,cs);
o2=orientation(rot2,cs);
count=0;
for i=1:length(o1)
%     ori=o1(i)*cs;
%     ori_ang=[ori.phi1;ori.Phi;ori.phi2]'*180/pi;
    for j=1:length(o2)
        ang1=angle(o1(i),o2(j))/degree;
%         ang_all(j)=angle(o1(i),o2(j))/degree;
        if ang1<1
            count=count+1;
%             if all(DS_ref.Dimension==DS_rec.Dimension)
%                 dis=sqrt(sum((DS_ref.Coord(i,:)-DS_rec.Coord(j,:)).^2)); % [pixel]
%             else
                COM_ref=(DS_ref.Coord(i,:)-DS_ref.Dimension/2)+DS_ref.Center'./DS_ref.VoxSize(1); % [pixel]
                COM_rec=(DS_rec.Coord(j,:)-DS_rec.Dimension/2)+DS_rec.Center'./DS_rec.VoxSize(1); % [pixel]
                dis=sqrt(sum((COM_ref-COM_rec).^2)); % [pixel]
%             end
            PairIndex(count,:)=[i j DS_ref.EulerZXZ(i,:) DS_rec.EulerZXZ(j,:) ...
                grainsize(i) grainsize_rec(j) DS_rec.SeedComp(j) ang1 dis];
        end
    end
%     i
    if mod(i,50)==0 || i==length(o1)
        sprintf('%d grains have been compared.', i)
    end
end
sprintf('Number of paired and uniquely paired grains in reference dataset and uniquely paired grains in rec dataset:')
[length(PairIndex(:,1)) length(unique(PairIndex(:,1))) length(unique(PairIndex(:,2)))]
sprintf('Unpaired grain IDs in reference dataset:')
setdiff(DS_ref.SeedID,PairIndex(:,1))
sprintf('Unpaired grain IDs in rec dataset:')
setdiff(DS_rec.SeedID,PairIndex(:,2))
problem_ID=setdiff([1:length(o2)],unique(PairIndex(:,2)));
Unpaired=[problem_ID' DS_rec.EulerZXZ(problem_ID,:) DS_rec.SeedComp(problem_ID,:) DS_rec.nVox(problem_ID)];



