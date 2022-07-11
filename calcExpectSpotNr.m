% calculate the number of spots expected on the detector
function NrSpotExpect=calcExpectSpotNr(Gt,hkl,grainpos,center,alpha,RotDet, ...
                SamposW,thetamax,lambda_min,lambda_max,Lsam2det,dety00,detz00, ...
                pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ)

v1 = [zeros(1,size(hkl,2));Gt(2,:);Gt(3,:)];
Glen = (Gt(1,:).^2 + Gt(2,:).^2 + Gt(3,:).^2).^0.5;
beta = acos(dot((ones(size(hkl))'.*grainpos/norm(grainpos))',Gt./Glen)); % [rad]

theta = beta-pi/2;
sintth = sin(2*theta);
costth = cos(2*theta);
d = 1./Glen*2*pi;
lambdahkl = 2 * d .*sin(theta);
Energy_hkl=12.398./lambdahkl; % [keV]

phix = acos(dot(v1./sqrt(v1(1,:).^2+v1(2,:).^2+v1(3,:).^2),(ones(size(hkl))'.*center/norm(center))'));
phiy = phix-2*theta;
L2 = (Lsam2det-SamposW(1))/cos(alpha);
diffvec = L2*sintth./sin(phiy); % [mm]
konst=sqrt(Gt(2,:).^2+Gt(3,:).^2);

dety22 = (center(2)+ (diffvec.*Gt(2,:)./konst)); % dety [mm]
detz22 = (center(3)+ (diffvec.*Gt(3,:)./konst)); % detz [mm]

K_out_unit = ([ones(1,size(hkl,2))*Lsam2det;dety22;detz22]-ones(1,size(hkl,2)).*SamposW) ...
    ./((Lsam2det-SamposW(1)).^2+(dety22-SamposW(2)).^2+(detz22-SamposW(3)).^2).^(1/2);
t = (RotDet(1,1)*(Lsam2det-SamposW(1))-RotDet(2,1)*(dety00-SamposW(2))-RotDet(3,1)*(detz00-SamposW(3)))./ ...
(RotDet(1,1)*K_out_unit(1,:)+RotDet(2,1)*K_out_unit(2,:)+RotDet(3,1)*K_out_unit(3,:));

dety22 = [RotDet(1,2) RotDet(2,2) RotDet(3,2)]*(t.*K_out_unit+[SamposW(1)-Lsam2det SamposW(2)-dety00 SamposW(3)-detz00]');
detz22 = [RotDet(1,3) RotDet(2,3) RotDet(3,3)]*(t.*K_out_unit+[SamposW(1)-Lsam2det SamposW(2)-dety00 SamposW(3)-detz00]');
% dety = -round(dety22/pixelysize-0.5)+dety0; % [pixel]
% detz = -round(detz22/pixelzsize-0.5)+detz0; % [pixel]
dety = -round(dety22/pixelysize)+dety0; % [pixel]
detz = -round(detz22/pixelzsize)+detz0; % [pixel]

select3=find(beta > pi/2 & beta < (90+thetamax*4)/180*pi & ...
    lambdahkl > lambda_min & lambdahkl < lambda_max & ...
    dety>=1 & dety<=detysize & detz>=1 & detz<=detzsize & ...
    ~(dety>=BeamStopY(1) & dety<=BeamStopY(2) & ...
        detz>=BeamStopZ(1) & detz<=BeamStopZ(2)));
if ~isempty(select3)
    NrSpotExpect=length(select3);
else
    NrSpotExpect=0;
end
    
        
        