function  rgb = ipf_color(U,symmetry_name,axis)
% Compute the IPF (inverse pole figure) colour for this orientation.
% symmetry_name: cubic, hexagonal, orthorhombic, tetragonal, triclinic
% axis = [0 0 1] by default, otherwise it can be [1 0 0],[0 1 0]
% adapt from pymicro/crystal/microstructure.py
% Jan 5, 2021

% % for testing
% x1=[5.27301095747675,8.91251170013055,349.797916410157];
% U=euler2u(x1(1)*pi/180,x1(2)*pi/180,x1(3)*pi/180); % mtex rgb = [0.9842	0.3982	0.1685]
% symmetry_name='cubic';
% axis=[0 0 1];

% default setting
if nargin<2
    symmetry_name='cubic';
    axis=[0 0 1];
elseif nargin<3
    axis=[0 0 1]; 
end
saturate_flag=0; % a flag to saturate the RGB values.
axis=axis./norm(axis);
if size(axis,2)>1
    axis=reshape(axis,[3 1]);
end
Vc=U'*axis; % note that U derived from euler2U is the inverse of actual U derived by conventional method

% get the symmetry operators
syms = symmetry_operator(symmetry_name);
syms_nr = size(syms,3);
for i=syms_nr+1:2*syms_nr
    syms(:,:,i)=-syms(:,:,i-syms_nr);
end
for i=1:size(syms,3)
    Vc_syms(:,i) = syms(:,:,i)*Vc;
end

% phi: rotation around 001 axis, from 100 axis to Vc vector, projected on (100,010) plane
Vc_phi=atan2(Vc_syms(2,:),Vc_syms(1,:))*180/pi; % [deg]
% chi: rotation around 010 axis, from 001 axis to Vc vector, projected on (100,001) plane
Vc_chi=atan2(Vc_syms(1,:),Vc_syms(3,:))*180/pi; % [deg]
% psi : angle from 001 axis to Vc vector
Vc_psi=acos(Vc_syms(3,:))*180/pi; % [deg]
switch symmetry_name
    case 'cubic'
        angleR = 45 - Vc_chi;  % red color proportional to (45 - chi)
        minAngleR = 0;
        maxAngleR = 45;
        angleB = Vc_phi;  % blue color proportional to phi
        minAngleB = 0;
        maxAngleB = 45;
    case 'hexagonal'
        angleR = 90 - Vc_psi; % red color proportional to (90 - psi)
        minAngleR = 0;
        maxAngleR = 90;
        angleB = Vc_phi; % blue color proportional to phi
        minAngleB = 0;
        maxAngleB = 30;
    otherwise
        error('Error: symmetry not supported');
end

% find the axis lying in the fundamental zone
i_SST = find(angleR >= minAngleR & angleR < maxAngleR & ...
                angleB >= minAngleB & angleB < maxAngleB);
if isempty(i_SST)
    error('problem moving to the fundamental zone');
end
r = angleR(i_SST) / maxAngleR;
g = (maxAngleR - angleR(i_SST)) / maxAngleR * (maxAngleB - angleB(i_SST)) / maxAngleB;
b = (maxAngleR - angleR(i_SST)) / maxAngleR * angleB(i_SST) / maxAngleB;
rgb = [r, g, b];

if saturate_flag==1
    rgb = rgb / max(rgb);
end

