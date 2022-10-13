% 1st step: fit the Lss, Lsd and detector shift
% 2nd step: fit the detector tilts
% Oct 28, 2021
function [ParaFit,fval]=L_shift_tilt_fitting(SpotsPair,S,B,P0y,P0z,RotAxisOffset,pixelysize,pixelzsize, ...
    dety0,detz0,Lsam2sou,Lsam2det,dety00,detz00,tilt_x,tilt_y,tilt_z,FitOption)

    switch FitOption
    case 'FitWithouTilts'
        % fit with Lss, Lsd, dety00, detz00 and no tilt
        clear x;
        x0=[Lsam2sou Lsam2det dety00 detz00];
        ErrMean=@(x)geo_fitting(x,SpotsPair,S,B,P0y,P0z,RotAxisOffset,pixelysize,pixelzsize,dety0,detz0,tilt_x,tilt_y,tilt_z);
        LB=[max([x0(1)-5 0]) max([x0(2)-5 0]) -pixelysize*150 -pixelzsize*150];
        UB=[x0(1)+5 x0(2)+5 pixelysize*150 pixelzsize*150];
        opts=optimset('Display','iter','Algorithm','trust-region-reflective', ...
            'MaxFunEvals',2000,'MaxIter',500,'TolFun',1e-10,'PlotFcns','optimplotfval');
        % opts=optimset('Display','iter','Algorithm','sqp','UseParallel',true, ...
        %     'MaxFunEvals',2000,'MaxIter',2000,'TolFun',1e-10,'PlotFcns','optimplotfval');
        [x,fval,exitflag,output] = fminsearchbnd(ErrMean,x0,LB,UB,opts);
        sprintf('The resulted residual is: %0.5f',fval)
        ParaFit=[x tilt_x tilt_y tilt_z];

    case 'FitAllOnce'
        % fit with Lss, Lsd, dety00, detz00 and tilt angles
        clear x;
        x0=[Lsam2sou Lsam2det dety00 detz00 tilt_x tilt_y tilt_z];
        ErrMean=@(x)geo_fitting(x,SpotsPair,S,B,P0y,P0z,RotAxisOffset,pixelysize,pixelzsize,dety0,detz0);
        stop_fitting=0;
        while stop_fitting~=1
            LB=[x0(1)-0.5 x0(2)-0.5 x0(3)-pixelysize*50 x0(4)-pixelysize*50 x0(5)-0.5 x0(6)-0.5 x0(7)-0.5];
            UB=[x0(1)+0.5 x0(2)+0.5 x0(3)+pixelysize*50 x0(4)+pixelysize*50 x0(5)+0.5 x0(6)+0.5 x0(7)+0.5];
%             LB=[max([x0(1)-1 0]) max([x0(2)-1 0]) dety00-pixelysize*100 detz00-pixelzsize*100 tilt_x-0.5 tilt_y-0.5 tilt_z-0.5];
%             UB=[x0(1)+1 x0(2)+1 dety00+pixelysize*100 detz00+pixelzsize*100 tilt_x+0.5 tilt_y+0.5 tilt_z+0.5];
    %         opts=optimset('Display','iter','Algorithm','trust-region-reflective', ...
    %             'MaxFunEvals',2000,'MaxIter',500,'TolFun',1e-10,'PlotFcns','optimplotfval');
            opts=optimset('Display','iter','Algorithm','trust-region-reflective', ...
                'MaxFunEvals',2000,'MaxIter',500,'TolFun',1e-10);
            % opts=optimset('Display','iter','Algorithm','sqp','UseParallel',true, ...
            %     'MaxFunEvals',2000,'MaxIter',2000,'TolFun',1e-10,'PlotFcns','optimplotfval');
%             opts=optimset('Display','iter','Algorithm',{'levenberg-marquardt',0.005}, ...
%                 'MaxFunEvals',2000,'MaxIter',2000,'TolFun',1e-10,'PlotFcns','optimplotfval');
            [x,fval,exitflag,output] = fminsearchbnd(ErrMean,x0,LB,UB,opts);
            if all(abs(x-LB)>0.1) && all(abs(x-UB)>0.1)
                stop_fitting=1;
            else
                x0=x;
            end
        end
        fprintf('The resulted residual is: %0.5f\n',fval)
        ParaFit=x;
    case 'FitAll2Steps'
        x0=[Lsam2sou Lsam2det P0y P0z dety00 detz00];
        ErrMean=@(x)L_shift_fitting(x,SpotsPair,S,B, ...
            pixelysize,pixelzsize,dety0,detz0,tilt_x,tilt_y,tilt_z,RotAxisOffset)
        LB=[max([x0(1)-5 0]) max([x0(2)-5 0]) x0(3)-0.02 x0(4)-0.02 -pixelysize*150 -pixelzsize*150];
        UB=[x0(1)+5 x0(2)+5 x0(3)+0.02 x0(4)+0.02 pixelysize*150 pixelzsize*150];
        opts=optimset('Display','iter','Algorithm','trust-region-reflective', ...
            'MaxFunEvals',2000,'MaxIter',500,'TolFun',1e-10,'PlotFcns','optimplotfval');
        % opts=optimset('Display','iter','Algorithm','sqp','UseParallel',true, ...
        %     'MaxFunEvals',2000,'MaxIter',2000,'TolFun',1e-10,'PlotFcns','optimplotfval');
        [x,fval1,exitflag,output] = fminsearchbnd(ErrMean,x0,LB,UB,opts);
        sprintf('The resulted residual is: %0.5f',fval1)
        ParaFit1=x

        clear x;
        x0=[tilt_x tilt_y tilt_z];
        ErrMean=@(x)det_tilt_fitting(x,SpotsPair,S,B, ...
            pixelysize,pixelzsize,dety0,detz0,ParaFit1(1),ParaFit1(2), ...
            ParaFit1(3),ParaFit1(4),RotAxisOffset,ParaFit1(5),ParaFit1(6));
        LB=[-5 -5 -5];
        UB=[5 5 5];
        opts=optimset('Display','iter','Algorithm','trust-region-reflective', ...
            'MaxFunEvals',2000,'MaxIter',500,'TolFun',1e-10,'PlotFcns','optimplotfval');
        % opts=optimset('Display','iter','Algorithm','sqp','UseParallel',true, ...
        %     'MaxFunEvals',2000,'MaxIter',2000,'TolFun',1e-10,'PlotFcns','optimplotfval');
        [x,fval,exitflag,output] = fminsearchbnd(ErrMean,x0,LB,UB,opts);
        sprintf('The resulted residual is: %0.5f',fval)
        ParaFit2=x;

        ParaFit=[ParaFit1 ParaFit2];
    end
end


% maximize the completeness to derive Lss and Lsd
% fit Lss, Lsd, dety00, detz00, tilt angles (optional)
% July 6, 2021
% updated on September 9, 2021
function ErrMean=L_shift_fitting(x,hittedSpots_pair,S,B, ...
    pixelysize,pixelzsize,dety0,detz0,tilt_x,tilt_y,tilt_z,RotAxisOffset)

Lsam2sou=x(1);
Lsam2det=x(2);
P0y=x(3)-RotAxisOffset;
P0z=x(4);
dety00=x(5)-RotAxisOffset;
detz00=x(6);
RotDet=get_det_R(tilt_x,tilt_y,tilt_z);

L=Lsam2sou+Lsam2det; % [mm]
Err=[];
for i=1:length(hittedSpots_pair(:,1)) 
    dety_exp=hittedSpots_pair(i,16); % [pixel]
    detz_exp=hittedSpots_pair(i,17); % [pixel]   
    
    omega=hittedSpots_pair(i,12)*pi/180; % [rad]
    Omega=[cos(omega) -sin(omega) 0;sin(omega) cos(omega) 0;0 0 1];
    hkl=[hittedSpots_pair(i,5) hittedSpots_pair(i,6) hittedSpots_pair(i,7)]';
    
    pos=hittedSpots_pair(i,2:4);
    pos(:,2)=pos(:,2)-RotAxisOffset;
    U=euler2u(hittedSpots_pair(i,8)*pi/180,hittedSpots_pair(i,9)*pi/180,hittedSpots_pair(i,10)*pi/180);

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
            t = (RotDet(1,1)*(Lsam2det-SamposW(1))+RotDet(2,1)*(dety00-SamposW(2))+RotDet(3,1)*(detz00-SamposW(3)))./ ...
                (RotDet(1,1)*K_out_unit(1)+RotDet(2,1)*K_out_unit(2)+RotDet(3,1)*K_out_unit(3));
            dety22 = [RotDet(1,2) RotDet(2,2) RotDet(3,2)]*(t*K_out_unit+[SamposW(1)-Lsam2det SamposW(2)-dety00 SamposW(3)-detz00]');
            detz22 = [RotDet(1,3) RotDet(2,3) RotDet(3,3)]*(t*K_out_unit+[SamposW(1)-Lsam2det SamposW(2)-dety00 SamposW(3)-detz00]');
%             dety = -round(dety22/pixelysize-0.5)+dety0; % [pixel]
%             detz = -round(detz22/pixelzsize-0.5)+detz0; % [pixel]
            dety = round(-dety22/pixelysize+dety0); % [pixel]
            detz = round(-detz22/pixelzsize+detz0); % [pixel]
            
            dis=sqrt((dety-dety_exp).^2+(detz-detz_exp).^2); % [mm]
            Err=[Err;dis]; % [mm]
%         end
%     end
end
ErrMean=mean(Err); % minimize the Euclidian distance [mm]
% ErrMean=sum(Err);
end


% maximize the completeness to derive Lss and Lsd
% fit Lss, Lsd, dety00, detz00, tilt angles (optional)
% July 6, 2021
% updated on September 9, 2021
function ErrMean=det_tilt_fitting(x,hittedSpots_pair,S,B, ...
    pixelysize,pixelzsize,dety0,detz0,Lsam2sou,Lsam2det, ...
    P0y,P0z,RotAxisOffset,dety00,detz00)

tilt_x=x(1);         % detector tilt counterclockwise around lab x axis [deg] 
tilt_y=x(2);         % detector tilt counterclockwise around lab y axis [deg] 
tilt_z=x(3);         % detector tilt counterclockwise around lab z axis [deg]
RotDet=get_det_R(tilt_x,tilt_y,tilt_z);

P0y=P0y-RotAxisOffset;
dety00=dety00-RotAxisOffset;

L=Lsam2sou+Lsam2det; % [mm]
Err=[];
for i=1:length(hittedSpots_pair(:,1)) 
    dety_exp=hittedSpots_pair(i,16); % [pixel]
    detz_exp=hittedSpots_pair(i,17); % [pixel]   
    
    omega=hittedSpots_pair(i,12)*pi/180; % [rad]
    Omega=[cos(omega) -sin(omega) 0;sin(omega) cos(omega) 0;0 0 1];
    hkl=[hittedSpots_pair(i,5) hittedSpots_pair(i,6) hittedSpots_pair(i,7)]';
    
    pos=hittedSpots_pair(i,2:4);
    pos(:,2)=pos(:,2)-RotAxisOffset;
    U=euler2u(hittedSpots_pair(i,8)*pi/180,hittedSpots_pair(i,9)*pi/180,hittedSpots_pair(i,10)*pi/180);

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
            t = (RotDet(1,1)*(Lsam2det-SamposW(1))+RotDet(2,1)*(dety00-SamposW(2))+RotDet(3,1)*(detz00-SamposW(3)))./ ...
                (RotDet(1,1)*K_out_unit(1)+RotDet(2,1)*K_out_unit(2)+RotDet(3,1)*K_out_unit(3));
            dety22 = [RotDet(1,2) RotDet(2,2) RotDet(3,2)]*(t*K_out_unit+[SamposW(1)-Lsam2det SamposW(2)-dety00 SamposW(3)-detz00]');
            detz22 = [RotDet(1,3) RotDet(2,3) RotDet(3,3)]*(t*K_out_unit+[SamposW(1)-Lsam2det SamposW(2)-dety00 SamposW(3)-detz00]');
            dety = round(-dety22/pixelysize+dety0); % [pixel]
            detz = round(-detz22/pixelzsize+detz0); % [pixel]
            
            dis=sqrt((dety-dety_exp).^2+(detz-detz_exp).^2); % [mm]
            Err=[Err;dis]; % [mm]
%         end
%     end
end
ErrMean=mean(Err); % minimize the Euclidian distance [mm]
% ErrMean=sum(Err);
end
