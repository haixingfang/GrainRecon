% fit RotAxisOffset
% updated on September 12, 2022
function ErrMean=geo_fitting_RotAxisOffset(x,hittedSpots_pair,S,B,P0y,P0z, ...
                        pixelysize,pixelzsize,dety0,detz0,Lsam2sou,Lsam2det,dety00, ...
                        detz00,tilt_x,tilt_y,tilt_z)

RotAxisOffset=x;
RotDet=get_det_R(tilt_x,tilt_y,tilt_z);

L=Lsam2sou+Lsam2det; % [mm]
Err=[];
Su_angle=[0 0 0];
Su=euler2u(Su_angle(1)*pi/180,Su_angle(2)*pi/180,Su_angle(3)*pi/180);
trans_pos=[0 0 0];
for i=1:length(hittedSpots_pair(:,1)) 
    dety_exp=hittedSpots_pair(i,12); % [pixel]
    detz_exp=hittedSpots_pair(i,13); % [pixel]   
    
    omega=hittedSpots_pair(i,5)*pi/180; % [rad]
    Omega=[cos(omega) -sin(omega) 0;sin(omega) cos(omega) 0;0 0 1];
    hkl=[hittedSpots_pair(i,6) hittedSpots_pair(i,7) hittedSpots_pair(i,8)]';

    pos=hittedSpots_pair(i,2:4);
    pos=pos+trans_pos;
    pos(:,2)=pos(:,2)-RotAxisOffset;
    U=euler2u(hittedSpots_pair(i,9)*pi/180,hittedSpots_pair(i,10)*pi/180,hittedSpots_pair(i,11)*pi/180);
    U=Su*U;

    SamposW=Omega*S*pos';
    center = [L, (SamposW(2)-P0y+RotAxisOffset)*L/(Lsam2sou+SamposW(1)), ...
        (SamposW(3)-P0z)*L/(Lsam2sou+SamposW(1))]; % sample center projected to the position of the detector
    alpha = atan(sqrt((SamposW(2)-P0y+RotAxisOffset)^2+(SamposW(3)-P0z)^2)/(Lsam2sou+SamposW(1)));
    grainpos = [Lsam2sou+SamposW(1) SamposW(2)-P0y+RotAxisOffset SamposW(3)-P0z];
    Gw = S*U*B*hkl;
    Gt=Omega*Gw;
    v1 = [0 Gt(2) Gt(3)];
    Glen = (Gt(1)^2 + Gt(2)^2 + Gt(3)^2)^0.5;
    beta = acos(dot(grainpos/norm(grainpos),Gt/Glen)); % [rad]
     
%     if beta > pi/2 && beta < (90+thetamax*4)/180*pi
        theta = beta-pi/2;
        sintth = sin(2*theta);
        costth = cos(2*theta);
        d = 1/Glen*2*pi;
        lambdahkl = 2 * d *sin(theta);
        Energy_hkl=12.398/lambdahkl; % [keV]    
%         if lambdahkl > lambda_min && lambdahkl < lambda_max
            phix = acos(dot(v1/norm(v1),center/norm(center)));
            phiy = phix-2*theta;
            L2 = (Lsam2det-SamposW(1))/cos(alpha);
            diffvec = L2*sintth/sin(phiy); % [mm]
            konst = norm([0 Gt(2) Gt(3)]);
            dety22 = (center(2)+ (diffvec*Gt(2)/konst)); % dety [mm]
            detz22 = (center(3)+ (diffvec*Gt(3)/konst)); % detz [mm]
            
            %%% include detector tilt and center offset
            K_out_unit = ([Lsam2det dety22 detz22]'-SamposW)./norm([Lsam2det dety22 detz22]'-SamposW);
            t = (RotDet(1,1)*(Lsam2det-SamposW(1))+RotDet(2,1)*(dety00-RotAxisOffset-SamposW(2))+RotDet(3,1)*(detz00-SamposW(3)))./ ...
                (RotDet(1,1)*K_out_unit(1)+RotDet(2,1)*K_out_unit(2)+RotDet(3,1)*K_out_unit(3));
            dety22 = [RotDet(1,2) RotDet(2,2) RotDet(3,2)]*(t*K_out_unit+[SamposW(1)-Lsam2det SamposW(2)-dety00+RotAxisOffset SamposW(3)-detz00]');
            detz22 = [RotDet(1,3) RotDet(2,3) RotDet(3,3)]*(t*K_out_unit+[SamposW(1)-Lsam2det SamposW(2)-dety00+RotAxisOffset SamposW(3)-detz00]');
            dety = round(-dety22/pixelysize+dety0); % [pixel]
            detz = round(-detz22/pixelzsize+detz0); % [pixel]
            
            dis=sqrt((dety-dety_exp).^2+(detz-detz_exp).^2); % [pixel]
            Err=[Err;dis]; % [mm]
%         end
%     end
end
ErrMean=mean(Err); % minimize the Euclidian distance [mm]
% ErrMean=sum(Err);
