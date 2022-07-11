% created by Haixing Fang on Jan 20, 2020
function [DS euler_grains grainvolume grainsize BoxDim grains rot o q GrainColor]=ReadGrainFrom_h5(FileFolder,h5FileName,Prefix_Name)

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

Dream3D_FileName=[Prefix_Name '.dream3d']; % create a new name for dream3D
DS = readLabDCT(fullfile(FileFolder,h5FileName)); % read h5 file that is exported by LabDCT [X*Y*Z]
dfile=fullfile(FileFolder,Dream3D_FileName); % file path for creating the dream3D file
if exist(dfile,'file')~=0
    delete(dfile);
end
Dream3DWriter(DS,dfile);
Xdmf_FileName=[Prefix_Name '.xdmf']; % create a new name for xdmf file
LabDCT_XdmfWriter(FileFolder,Dream3D_FileName,Xdmf_FileName);

euler = Euler(rodrigues2quat(vector3d(double(DS.RodVec'))));
euler_grains=euler.*180/pi;% Euler angles [deg]: ([0 2*pi], [0 pi], [0 2*pi])
grainvolume=DS.nVox*DS.VoxSize(1)*DS.VoxSize(2)*DS.VoxSize(3)*1e9; % [um^3]
grainsize=2*(3*grainvolume/(4*pi)).^(1/3); % equivalent diameter of grain [um]
BoxDim=DS.Dimension.*DS.VoxSize'*1000; % [um]
grains=length(nonzeros(grainsize));

rot = rotation('Euler',euler_grains(:,1)*degree,euler_grains(:,2)*degree,euler_grains(:,3)*degree);
o = orientation(rot,cs,ss);
q=euler2quat(euler_grains*degree);
GrainColor=IPFplot(o,cs,ss,1); % plot IPF001



