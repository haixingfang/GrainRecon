% grain verification
% June 14, 2021
function [Nr_simu,Nr_intersect,dis_median,SimuSpots,HittedSpots]=index_verify_v3(U,proj_bin_bw,Spots,pos,rot_angles,S,B,Ahkl,nrhkl,hklnumber,hkl_square, ...
                RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
                pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ)

L=Lsam2sou+Lsam2det;
Nr_simu=0;
Nr_intersect=0;
dis_mean=0;
HittedSpots=[];
SimuSpots=[];

hkl = [Ahkl(1:nrhkl,1) Ahkl(1:nrhkl,2) Ahkl(1:nrhkl,3)]';
Gw = S*U*B*hkl;
for i=1:length(rot_angles)    
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
%     dety = -round(dety22/pixelysize-0.5)+dety0; % [pixel]
%     detz = -round(detz22/pixelzsize-0.5)+detz0; % [pixel]
    dety = -round(dety22/pixelysize)+dety0; % [pixel]
    detz = -round(detz22/pixelzsize)+detz0; % [pixel]

    select3=find(beta > pi/2 & beta < (90+thetamax*4)/180*pi &...
        lambdahkl > lambda_min & lambdahkl < lambda_max & ...
        dety>=1 & dety<=detysize & detz>=1 & detz<=detzsize & ...
        ~(dety>=BeamStopY(1) & dety<=BeamStopY(2) & ...
            detz>=BeamStopZ(1) & detz<=BeamStopZ(2)));
    if ~isempty(select3)
        Nr_simu=Nr_simu+length(select3);
        for j=1:length(select3)
            dis = proj_bin_bw(detz(select3(j)),dety(select3(j)),i); % euclidian [pixel]
            dis_mean=dis_mean+dis;
            if dis <= minEucDis/mean([pixelysize pixelzsize])
                Nr_intersect=Nr_intersect+1;
                % pairing the spots
                [minDis,spotpair_exp_ID]=min(sqrt((Spots{i}(:,9)-dety(select3(j))).^2+ ...
                    (Spots{i}(:,10)-detz(select3(j))).^2)); % [pixel]
                if minDis<min([sqrt(Spots{i}(spotpair_exp_ID,11)) 15])
                    HittedSpots=[HittedSpots;reshape(U,1,9) pos hkl(:,select3(j))' Energy_hkl(select3(j)) ...
                        rot_angles(i) Spots{i}(spotpair_exp_ID,9:10) minDis dety(select3(j)) detz(select3(j)) dis];
                end
            end
            SimuSpots=[SimuSpots;hkl(:,select3(j))' Energy_hkl(select3(j)) rot_angles(i) ...
                dety22(select3(j)) detz22(select3(j)) dis dety(select3(j)) detz(select3(j))];
        end
    end
end
% dis_mean=dis_mean/Nr_intersect;
% dis_mean=dis_mean/Nr_simu;
if ~isempty(SimuSpots)
    dis_median=median(SimuSpots(:,8));
else
    dis_median=2000;
end
            
       
            
%{
for i=1:length(rot_angles)
    omega=rot_angles(i); % [rad]
    Omega=[cos(omega) -sin(omega) 0;sin(omega) cos(omega) 0;0 0 1];
    for j=1:nrhkl
    %     for j=find(Ahkl(:,1)==-1 & Ahkl(:,2)==1 & Ahkl(:,3)==3)
        if (Ahkl(j,1)^2+Ahkl(j,2)^2+Ahkl(j,3)^2)<=hkl_square(hklnumber)
            hkl = [Ahkl(j,1) Ahkl(j,2) Ahkl(j,3)]';
%             dhkl=sqrt(1./sum(hkl(1:3)'.^2./cell(1:3).^2)); % for orthogonal [angstrom]
    
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
            if beta > pi/2 && beta < (90+thetamax*4)/180*pi
                theta = beta-pi/2;
                sintth = sin(2*theta);
                costth = cos(2*theta);
                d = 1/Glen*2*pi;
                lambdahkl = 2 * d *sin(theta);
                Energy_hkl=12.398/lambdahkl; % [keV]
                if lambdahkl > lambda_min && lambdahkl < lambda_max
                    phix = acos(dot(v1/norm(v1),center/norm(center)));
                    phiy = phix-2*theta;
                    L2 = (Lsam2det-SamposW(1))/cos(alpha);
                    diffvec = L2*sintth/sin(phiy); % [mm]
                    konst = norm([0 Gt(2) Gt(3)]);
                    dety22 = (center(2)+ (diffvec*Gt(2)/konst)); % dety [mm]
                    detz22 = (center(3)+ (diffvec*Gt(3)/konst)); % detz [mm]
                    
                    %%% include detector tilt
                    K_out_unit = ([Lsam2det dety22 detz22]'-SamposW)./norm([Lsam2det dety22 detz22]'-SamposW);
                    t = (RotDet(1,1)*(Lsam2det-SamposW(1))-RotDet(2,1)*(dety00-SamposW(2))-RotDet(3,1)*(detz00-SamposW(3)))./ ...
                    (RotDet(1,1)*K_out_unit(1)+RotDet(2,1)*K_out_unit(2)+RotDet(3,1)*K_out_unit(3));
                    dety22 = [RotDet(1,2) RotDet(2,2) RotDet(3,2)]*(t*K_out_unit+[SamposW(1)-Lsam2det SamposW(2)-dety00 SamposW(3)-detz00]');
                    detz22 = [RotDet(1,3) RotDet(2,3) RotDet(3,3)]*(t*K_out_unit+[SamposW(1)-Lsam2det SamposW(2)-dety00 SamposW(3)-detz00]');
                    dety = -round(dety22/pixelysize+0.5)+dety0; % [pixel]
                    detz = -round(detz22/pixelzsize+0.5)+detz0; % [pixel]
                    
%                     if (dety22>=-(detysize-dety0)*pixelysize && dety22<=(detysize-dety0)*pixelysize) && ...
%                             (detz22>=-(detzsize-detz0)*pixelzsize && detz22<=(detzsize-detz0)*pixelzsize) && ...
%                             ~((dety22>=(BeamStopY(1)-dety0)*pixelysize && dety22<=(BeamStopY(2)-dety0)*pixelysize) && ...
%                             (detz22>=(BeamStopZ(1)-detz0)*pixelzsize && detz22<=(BeamStopZ(2)-detz0)*pixelzsize))
                    if (dety>=1 && dety<=detysize && detz>=1 && detz<=detzsize) ...
                       && ~(dety>=dety0-BeamStopY(1) && dety<=dety0+BeamStopY(2) && ...
                            detz>=detz0-BeamStopZ(1) && detz<=detz0+BeamStopZ(2))                        
                        Nr_simu=Nr_simu+1;
                        dis = proj_bin_bw(detz,dety,i); % euclidian [pixel]

                        dis_mean=dis_mean+dis;
                        if dis <= minEucDis/mean([pixelysize pixelzsize])
                            Nr_intersect=Nr_intersect+1;
                            HittedSpots=[HittedSpots;hkl' Energy_hkl rot_angles(i) dety22 detz22 dis];
                        elseif dis <= 2*minEucDis/mean([pixelysize pixelzsize])
                            HittedSpots=[HittedSpots;hkl' Energy_hkl rot_angles(i) dety22 detz22 dis];
                        end
                        SimuSpots=[SimuSpots;hkl' Energy_hkl rot_angles(i) dety22 detz22 dis];
                    end
                end
            end
        end
    end % end of nrhkl
end
% dis_mean=dis_mean/Nr_intersect;
% dis_mean=dis_mean/Nr_simu;
dis_mean=median(SimuSpots(:,8));
%}

