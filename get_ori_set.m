function OR=get_ori_set(OR_folder,sgno,angular_resolution)
% generate orientation via solution of the hyperspherical covering problem
% Peter Mahler Larsen and Søren Schmidt
% Improved orientation sampling for indexing diffraction patterns of polycrystalline materials
% Journal of Applied Crystallography. (2017). 50, 1571-1582

if nargin<3
    angular_resolution='2';
end
if sgno>=221 % 4/m -3 2/m, Schoenflies: Oh, cubic materials such as Fe, Al, Si, Ni etc.
    switch angular_resolution
        case '3'
            q=dlmread(fullfile(OR_folder,'Ori_Oh_3deg.txt'),' ',1,0); % 9218 for 3 deg
        case '2'
            q=dlmread(fullfile(OR_folder,'Ori_Oh_2deg.txt'),' ',1,0); % 32768 for 2 deg
        case '1'
            q=dlmread(fullfile(OR_folder,'Ori_Oh_1deg.txt'),' ',1,0); % 262144 for 1 deg
    end
elseif sgno>=191 && sgno<=194 % 6/m 2/m 2/m, Schoenflies: D6h, hexagonal such as Mg-194, Zn-194, Zr-194
    switch angular_resolution
        case '3'
            q=dlmread(fullfile(OR_folder,'Ori_D6h_3deg.txt'),' ',1,0); % 9218 for 3 deg
        case '2'
            q=dlmread(fullfile(OR_folder,'Ori_D6h_2deg.txt'),' ',1,0); % 32768 for 2 deg
        case '1'
            q=dlmread(fullfile(OR_folder,'Ori_D6h_1deg.txt'),' ',1,0); % 262144 for 1 deg
    end
elseif sgno>=200 && sgno<=206 % 2/m-3, Schoenflies: Th
    switch angular_resolution
        case '3'
            q=dlmread(fullfile(OR_folder,'Ori_Th_3deg.txt'),' ',1,0); % 9218 for 3 deg
        case '2'
            q=dlmread(fullfile(OR_folder,'Ori_Th_2deg.txt'),' ',1,0); % 32768 for 2 deg
        case '1'
            q=dlmread(fullfile(OR_folder,'Ori_Th_1deg.txt'),' ',1,0); % 262144 for 1 deg
    end
elseif sgno>=10 && sgno<=15     % 2/m, Schoenflies: C2h
    switch angular_resolution
        case '3'
            q=dlmread(fullfile(OR_folder,'Ori_C2h_3deg.txt'),' ',1,0); % 9218 for 3 deg
        case '2'
            q=dlmread(fullfile(OR_folder,'Ori_C2h_2deg.txt'),' ',1,0); % 32768 for 2 deg
        case '1'
            q=dlmread(fullfile(OR_folder,'Ori_C2h_1deg.txt'),' ',1,0); % 262144 for 1 deg
    end
else
    error('No orientation set is available currently');
end
[euler_angles(:,1), euler_angles(:,2), euler_angles(:,3)] = quat2angle(q(:,1:4),'ZXZ');

euler_angles(euler_angles(:,1)<0,1)=euler_angles(euler_angles(:,1)<0,1)+2*pi;
euler_angles(euler_angles(:,2)<0,2)=euler_angles(euler_angles(:,2)<0,2)+2*pi;
euler_angles(euler_angles(:,3)<0,3)=euler_angles(euler_angles(:,3)<0,3)+2*pi;
rod=angle2rod(euler_angles(:,1),euler_angles(:,2),euler_angles(:,3),'ZXZ');
euler_angles=euler_angles*180/pi;

OR.len=length(q(:,1));
OR.euler_angles=euler_angles;
OR.rod=rod;
OR.q=q(:,1:4);
OR.q_b=q(:,2);
OR.q_c=q(:,3);
OR.q_d=q(:,4);
