% minimize the difference in distances between simu and exp
% fit the grain center of mass
% input: x is the position of the voxel [mm], euler angle fitting is optional
% Dec 3, 2021
function [FitOutput,fval]=fit_grain_COM(hittedSpots_pair,S,B,Lsam2sou,Lsam2det,P0y,P0z, ...
            pixelysize,pixelzsize,dety0,detz0,dety00,detz00,RotDet)
    U=reshape(hittedSpots_pair(1,2:10),3,3);
    clear x;
    fit_grain_COM_only=0;
    if fit_grain_COM_only==1
        x0=hittedSpots_pair(1,11:13); % position
        ErrMean=@(x)fit_grain_COM_fun(x,hittedSpots_pair,S,B,Lsam2sou,Lsam2det,P0y,P0z, ...
            pixelysize,pixelzsize,dety0,detz0,dety00,detz00,RotDet,U);
        LB=x0-0.1;
        UB=x0+0.1;
    else
        x0=[hittedSpots_pair(1,11:13) u2euler_corr(U)]; % position + euler_angle
        ErrMean=@(x)fit_grain_COM_fun(x,hittedSpots_pair,S,B,Lsam2sou,Lsam2det,P0y,P0z, ...
            pixelysize,pixelzsize,dety0,detz0,dety00,detz00,RotDet);
        LB=x0-0.1;
        UB=x0+0.1;
        LB=[x0(1:3)-0.1 max([x0(4)-3 0]) max([x0(5)-3 0]) max([x0(6)-3 0])];
        UB=[x0(1:3)+0.1 min([x0(4)+3 360]) min([x0(5)+3 360]) min([x0(6)+3 360])];
    end
    opts=optimset('Display','off','Algorithm','trust-region-reflective', ...
        'MaxFunEvals',2000,'MaxIter',500,'TolFun',1e-10);
    [FitOutput,fval,~,~] = fminsearchbnd(ErrMean,x0,LB,UB,opts);
%     sprintf('The resulted residual is: %0.5f pixel',fval)
end


function ErrMean=fit_grain_COM_fun(x,hittedSpots_pair,S,B,Lsam2sou,Lsam2det,P0y,P0z, ...
    pixelysize,pixelzsize,dety0,detz0,dety00,detz00,RotDet,U)

pos=[x(1) x(2) x(3)]; % [mm]
if nargin<=15
    U=euler2u(x(4)*pi/180,x(5)*pi/180,x(6)*pi/180);
end
L=Lsam2sou+Lsam2det; % [mm]
Err=[];
for i=1:length(hittedSpots_pair(:,1)) 
    dety_exp=hittedSpots_pair(i,19); % [pixel]
    detz_exp=hittedSpots_pair(i,20); % [pixel]   
    
    omega=hittedSpots_pair(i,18)*pi/180; % [rad]
    Omega=[cos(omega) -sin(omega) 0;sin(omega) cos(omega) 0;0 0 1];
    hkl=[hittedSpots_pair(i,14) hittedSpots_pair(i,15) hittedSpots_pair(i,16)]';

    SamposW=Omega*S*pos';
    center = [L, (SamposW(2)-P0y)*L/(Lsam2sou+SamposW(1)), ...
        (SamposW(3)-P0z)*L/(Lsam2sou+SamposW(1))]; % sample center projected to the position of the detector
    alpha = atan(sqrt((SamposW(2)-P0y)^2+(SamposW(3)-P0z)^2)/(Lsam2sou+SamposW(1)));
    grainpos = [Lsam2sou+SamposW(1) SamposW(2)-P0y SamposW(3)-P0z];
    Gw = S*U*B*hkl;
    Gt=Omega*Gw;
    v1 = [0 Gt(2) Gt(3)];
    Glen = (Gt(1)^2 + Gt(2)^2 + Gt(3)^2)^0.5;
    beta = acos(dot(grainpos/norm(grainpos),Gt/Glen)); % [rad]
     
    theta = beta-pi/2;
    sintth = sin(2*theta);
    costth = cos(2*theta);
    d = 1/Glen*2*pi;
    lambdahkl = 2 * d *sin(theta);
    Energy_hkl=12.398/lambdahkl; % [keV]    
    phix = acos(dot(v1/norm(v1),center/norm(center)));
    phiy = phix-2*theta;
    L2 = (Lsam2det-SamposW(1))/cos(alpha);
    diffvec = L2*sintth/sin(phiy); % [mm]
    konst = norm([0 Gt(2) Gt(3)]);
    dety22 = (center(2)+ (diffvec*Gt(2)/konst)); % dety [mm]
    detz22 = (center(3)+ (diffvec*Gt(3)/konst)); % detz [mm]
            
    %%% include detector tilt and center offset
    K_out_unit = ([Lsam2det dety22 detz22]'-SamposW)./norm([Lsam2det dety22 detz22]'-SamposW);
    t = (RotDet(1,1)*(Lsam2det-SamposW(1))-RotDet(2,1)*(dety00-SamposW(2))-RotDet(3,1)*(detz00-SamposW(3)))./ ...
    (RotDet(1,1)*K_out_unit(1)+RotDet(2,1)*K_out_unit(2)+RotDet(3,1)*K_out_unit(3));
    dety22 = [RotDet(1,2) RotDet(2,2) RotDet(3,2)]*(t*K_out_unit+[SamposW(1)-Lsam2det SamposW(2)-dety00 SamposW(3)-detz00]');
    detz22 = [RotDet(1,3) RotDet(2,3) RotDet(3,3)]*(t*K_out_unit+[SamposW(1)-Lsam2det SamposW(2)-dety00 SamposW(3)-detz00]');
%     dety = -round(dety22/pixelysize-0.5)+dety0; % [pixel]
%     detz = -round(detz22/pixelzsize-0.5)+detz0; % [pixel]
    dety = -round(dety22/pixelysize)+dety0; % [pixel]
    detz = -round(detz22/pixelzsize)+detz0; % [pixel]
            
    dis=sqrt((dety-dety_exp).^2+(detz-detz_exp).^2); % [pixel]
    Err=[Err;dis]; % [pixel]
end
ErrMean=mean(Err); % minimize the Euclidian distance [mm]
end


