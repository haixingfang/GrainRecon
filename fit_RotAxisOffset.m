% fit RotAxisOffset
function [RotAxisOffset_fitted,fval]=fit_RotAxisOffset(RotAxisOffset,hittedSpots_pair,S,B,P0y,P0z, ...
                        pixelysize,pixelzsize,dety0,detz0,Lsam2sou,Lsam2det,dety00, ...
                        detz00,tilt_x,tilt_y,tilt_z)
                    
x0=RotAxisOffset;
ErrMean=@(x)fit_RotAxisOffset_fun(x,hittedSpots_pair,S,B,P0y,P0z, ...
                        pixelysize,pixelzsize,dety0,detz0,Lsam2sou,Lsam2det,dety00, ...
                        detz00,tilt_x,tilt_y,tilt_z);
stop_fitting=0;
while stop_fitting~=1
    LB=x0-3*pixelysize/(Lsam2det/Lsam2sou+1);
    UB=x0(1)+3*pixelysize/(Lsam2det/Lsam2sou+1);
%         opts=optimset('Display','iter','Algorithm','trust-region-reflective', ...
%             'MaxFunEvals',2000,'MaxIter',500,'TolFun',1e-10,'PlotFcns','optimplotfval');
    opts=optimset('Display','iter','Algorithm','trust-region-reflective', ...
        'MaxFunEvals',2000,'MaxIter',500,'TolFun',1e-10);
    % opts=optimset('Display','iter','Algorithm','sqp','UseParallel',true, ...
    %     'MaxFunEvals',2000,'MaxIter',2000,'TolFun',1e-10,'PlotFcns','optimplotfval');
%             opts=optimset('Display','iter','Algorithm',{'levenberg-marquardt',0.005}, ...
%                 'MaxFunEvals',2000,'MaxIter',2000,'TolFun',1e-10,'PlotFcns','optimplotfval');
    [x,fval,exitflag,output] = fminsearchbnd(ErrMean,x0,LB,UB,opts);
    if all(abs(x-LB)>0.01*pixelysize/(Lsam2det/Lsam2sou+1)) && all(abs(x-UB)>0.01*pixelysize/(Lsam2det/Lsam2sou+1))
        stop_fitting=1;
    else
        x0=x;
    end
end
sprintf('The resulted residual is: %0.5f',fval)
RotAxisOffset_fitted=x;
end
        
function ErrMean=fit_RotAxisOffset_fun(x,hittedSpots_pair,S,B,P0y,P0z, ...
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
    dety_exp=hittedSpots_pair(i,16); % [pixel]
    detz_exp=hittedSpots_pair(i,17); % [pixel]   
    
    omega=hittedSpots_pair(i,12)*pi/180; % [rad]
    Omega=[cos(omega) -sin(omega) 0;sin(omega) cos(omega) 0;0 0 1];
    hkl=[hittedSpots_pair(i,5) hittedSpots_pair(i,6) hittedSpots_pair(i,7)]';

    pos=hittedSpots_pair(i,2:4);
    pos=pos+trans_pos;
    pos(:,2)=pos(:,2)-RotAxisOffset;
    U=euler2u(hittedSpots_pair(i,8)*pi/180,hittedSpots_pair(i,9)*pi/180,hittedSpots_pair(i,10)*pi/180);
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
end
