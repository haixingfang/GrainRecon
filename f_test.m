
function Err1sum=f_test(x,U,proj_bin_bw,rot_angles,S,B,RotDet,lambda_min,lambda_max, ...
    Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
    pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ)

pos=x;
% hkl=[0 0 2]';
% hkl=[0 2 -2]';
% hkl=[-1 -1 -1]';

hkl=[0 0 2;0 2 -2;-1 -1 -1]';
exp_spot=[744 201;603 1641;1716 1344];

Nr_simu=0;
Nr_intersect=0;
Err1=[];
Err2=[];
L=Lsam2sou+Lsam2det;

Gw = S*U*B*hkl;
i=1;
% for i=1:length(rot_angles)    
    omega=rot_angles(i)*pi/180; % [rad]
    Omega(:,:,i)=[cos(omega) -sin(omega) 0;sin(omega) cos(omega) 0;0 0 1];
    
    SamposW=Omega(:,:,i)*S*pos';
    center = [L, (SamposW(2)-P0y)*L/(Lsam2sou+SamposW(1)), ...
        (SamposW(3)-P0z)*L/(Lsam2sou+SamposW(1))]; % sample center projected to the position of the detector
    alpha = atan(sqrt((SamposW(2)-P0y)^2+(SamposW(3)-P0z)^2)/(Lsam2sou+SamposW(1)));
    grainpos = [Lsam2sou+SamposW(1) SamposW(2)-P0y SamposW(3)-P0z];
    
    Gt=Omega(:,:,i)*Gw;
    v1 = [zeros(1,size(hkl,2));Gt(2,:);Gt(3,:)];
    Glen = (Gt(1,:).^2 + Gt(2,:).^2 + Gt(3,:).^2).^0.5;
    beta = acos(dot((ones(size(hkl))'.*grainpos/norm(grainpos))',Gt./Glen)); % [rad]
    
%     select1=find(beta > pi/2 & beta < (90+thetamax*4)/180*pi);
    
    theta = beta-pi/2;
    sintth = sin(2*theta);
    costth = cos(2*theta);
    d = 1./Glen*2*pi;
    lambdahkl = 2 * d .*sin(theta);
    Energy_hkl=12.398./lambdahkl; % [keV]
    
%     select2=find(lambdahkl > lambda_min & lambdahkl < lambda_max);
    
    phix = acos(dot(v1./sqrt(v1(1,:).^2+v1(2,:).^2+v1(3,:).^2),(ones(size(hkl))'.*center/norm(center))'));
    phiy = phix-2*theta;
    L2 = (Lsam2det-SamposW(1))/cos(alpha);
    diffvec = L2*sintth./sin(phiy); % [mm]
    konst=sqrt(Gt(2,:).^2+Gt(3,:).^2);
    
    dety22 = (center(2)+ (diffvec.*Gt(2,:)./konst)); % dety [mm]
    detz22 = (center(3)+ (diffvec.*Gt(3,:)./konst)); % detz [mm]
    
    K_out_unit = ([ones(1,size(hkl,2))*Lsam2det;dety22;detz22]-ones(1,size(hkl,2)).*SamposW) ...
        ./((Lsam2det-SamposW(1)).^2+(dety22-SamposW(2)).^2+(detz22-SamposW(3)).^2).^(1/2);
    t = (RotDet(1,1)*(Lsam2det-SamposW(1))+RotDet(2,1)*(dety00-SamposW(2))+RotDet(3,1)*(detz00-SamposW(3)))./ ...
        (RotDet(1,1)*K_out_unit(1,:)+RotDet(2,1)*K_out_unit(2,:)+RotDet(3,1)*K_out_unit(3,:));

    dety22 = [RotDet(1,2) RotDet(2,2) RotDet(3,2)]*(t.*K_out_unit+[SamposW(1)-Lsam2det SamposW(2)-dety00 SamposW(3)-detz00]');
    detz22 = [RotDet(1,3) RotDet(2,3) RotDet(3,3)]*(t.*K_out_unit+[SamposW(1)-Lsam2det SamposW(2)-dety00 SamposW(3)-detz00]');
    dety = round(-dety22/pixelysize+dety0); % [pixel]
    detz = round(-detz22/pixelzsize+detz0); % [pixel]

    select3=find(beta > pi/2 & ...
    lambdahkl > lambda_min & lambdahkl < lambda_max & ...
    dety>=1 & dety<=detysize & detz>=1 & detz<=detzsize & ...
    ~(dety>=BeamStopY(1) & dety<=BeamStopY(2) & ...
        detz>=BeamStopZ(1) & detz<=BeamStopZ(2)));
    select3=1:length(dety);
    if ~isempty(select3)
        Nr_simu=Nr_simu+length(select3);
        for j=1:length(select3)
            dis = proj_bin_bw(detz(select3(j)),dety(select3(j)),i); % euclidian [pixel]
            if dis <= minEucDis/mean([pixelysize pixelzsize])
                Nr_intersect=Nr_intersect+1;
            end
            Err2=[Err2;dis]; % [mm]
            
            dis0=sqrt((dety(j)-exp_spot(j,1))^2+(detz(j)-exp_spot(j,2))^2);
            Err1=[Err1;dis0];
        end
    end
    Err1sum=sum(Err1);
% end


