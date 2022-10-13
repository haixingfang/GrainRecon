% grain verification
% June 14, 2021
% calculate completeness and Icorr given pos and U
function [Nr_simu,Nr_intersect,Icorr,Iexp_sum,dis_median,SimuSpots,HittedSpots]=index_verify_grey(U,proj,proj_bin_bw,pos,rot_angles,S,B,Ahkl,nrhkl, ...
                RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
                RotAxisOffset,pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ, ...
                K1,V,I0E,Energy,ExpTime,atomparam,Transmission,rou,Rsample,CsI,Swank,lambda,VoxSize)
% % for testing
% pos=SubGrain{grainno}(67,2:4); % [mm]

L=Lsam2sou+Lsam2det;
Nr_simu=0;
Nr_intersect=0;
dis_mean=0;
HittedSpots=[];
SimuSpots=[];
Imulti_sum=0;
Icalc_sum2=0;
Iexp_sum2=0;
Iexp_sum=0;

hkl = [Ahkl(1:nrhkl,1) Ahkl(1:nrhkl,2) Ahkl(1:nrhkl,3)]';
Gw = S*U*B*hkl;
pos(:,2)=pos(:,2)-RotAxisOffset;
for i=1:length(rot_angles)    
    omega=rot_angles(i)*pi/180; % [rad]
    Omega(:,:,i)=[cos(omega) -sin(omega) 0;sin(omega) cos(omega) 0;0 0 1];
    
    SamposW=Omega(:,:,i)*S*pos';
    center = [L, (SamposW(2)-P0y+RotAxisOffset)*L/(Lsam2sou+SamposW(1)), ...
        (SamposW(3)-P0z)*L/(Lsam2sou+SamposW(1))]; % sample center projected to the position of the detector
    alpha = atan(sqrt((SamposW(2)-P0y+RotAxisOffset)^2+(SamposW(3)-P0z)^2)/(Lsam2sou+SamposW(1)));
    grainpos = [Lsam2sou+SamposW(1) SamposW(2)-P0y+RotAxisOffset SamposW(3)-P0z];
    
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
    t = (RotDet(1,1)*(Lsam2det-SamposW(1))+RotDet(2,1)*(dety00-RotAxisOffset-SamposW(2))+RotDet(3,1)*(detz00-SamposW(3)))./ ...
        (RotDet(1,1)*K_out_unit(1,:)+RotDet(2,1)*K_out_unit(2,:)+RotDet(3,1)*K_out_unit(3,:));

    dety22 = [RotDet(1,2) RotDet(2,2) RotDet(3,2)]*(t.*K_out_unit+[SamposW(1)-Lsam2det SamposW(2)-dety00+RotAxisOffset SamposW(3)-detz00]');
    detz22 = [RotDet(1,3) RotDet(2,3) RotDet(3,3)]*(t.*K_out_unit+[SamposW(1)-Lsam2det SamposW(2)-dety00+RotAxisOffset SamposW(3)-detz00]');
%     dety = -round(dety22/pixelysize-0.5)+dety0; % [pixel]
%     detz = -round(detz22/pixelzsize-0.5)+detz0; % [pixel]
    dety = round(-dety22/pixelysize+dety0); % [pixel]
    detz = round(-detz22/pixelzsize+detz0); % [pixel]
    
    select3=find(beta > pi/2 & beta < (90+thetamax*4)/180*pi &...
        lambdahkl > lambda_min & lambdahkl < lambda_max & ...
        dety>=1 & dety<=detysize & detz>=1 & detz<=detzsize & ...
        ~(dety>=BeamStopY(1) & dety<=BeamStopY(2) & ...
            detz>=BeamStopZ(1) & detz<=BeamStopZ(2)));
    if ~isempty(select3)
        Nr_simu=Nr_simu+length(select3);
        for j=1:length(select3)
            % calculate the intensity
            Lorentz=1./(sin(2*theta(select3(j)))); % Lorentz factor
            P=(1+costth(select3(j))^2)/2; % polarization factor
            ee=min(find(Energy>(Energy_hkl(select3(j))-(Energy(2)-Energy(1))) & Energy<(Energy_hkl(select3(j))+(Energy(2)-Energy(1)))));
            [A_Ehkl, L_total]=beam_attenuation(SamposW,Lsam2sou,Lsam2det,dety(select3(j)),detz(select3(j)), ...
                atomparam,Transmission,rou,Energy_hkl(select3(j)),Rsample);       % attenuation intensity factor, June 22, 2020
            [DQE_Ehkl]=Detector_efficiency(CsI,Swank,Energy_hkl(select3(j))); % DQE, June 22, 2020 
            K2 = lambda(ee)^3*VoxSize^3*10^21/V^2;
            K2 = A_Ehkl*DQE_Ehkl*K2; % consider attenuation and detector efficiency
            Icalc=K1*K2*abs(I0E(ee))*Lorentz*P*Ahkl(select3(j),5)*ExpTime; % calculated pixel intensity [photons]
            Iexp=double(proj{i}(detz(select3(j)),dety(select3(j)))); % experimental pixel intensity
            Imulti_sum=Imulti_sum+Icalc*Iexp;           
            Icalc_sum2=Icalc_sum2+Icalc^2;
            Iexp_sum2=Iexp_sum2+Iexp^2;
            Iexp_sum=Iexp_sum+Iexp; % sum of experimental pixel intensities
            
            dis = proj_bin_bw(detz(select3(j)),dety(select3(j)),i); % euclidian [pixel]
            dis_mean=dis_mean+dis;
            if dis <= minEucDis/mean([pixelysize pixelzsize])
                Nr_intersect=Nr_intersect+1;
                HittedSpots=[HittedSpots;hkl(:,select3(j))' Energy_hkl(select3(j)) ...
                    rot_angles(i) dety22(select3(j)) detz22(select3(j)) dis dety(select3(j)) detz(select3(j))];
            elseif dis <= 2*minEucDis/mean([pixelysize pixelzsize])
                HittedSpots=[HittedSpots;hkl(:,select3(j))' Energy_hkl(select3(j)) ...
                    rot_angles(i) dety22(select3(j)) detz22(select3(j)) dis dety(select3(j)) detz(select3(j))];
            end
            SimuSpots=[SimuSpots;hkl(:,select3(j))' Energy_hkl(select3(j)) rot_angles(i) ...
                dety22(select3(j)) detz22(select3(j)) dis dety(select3(j)) detz(select3(j))];
        end
    end
end
% calculate the intensity correlation, the same as E. Rauch, Z. Kristallogr. 2009
Icorr=Imulti_sum/(sqrt(Icalc_sum2)*sqrt(Iexp_sum2));

% dis_mean=dis_mean/Nr_intersect;
% dis_mean=dis_mean/Nr_simu;
if ~isempty(SimuSpots)
    dis_median=median(SimuSpots(:,8));
else
    dis_median=2000;
end
      