% discretize the quaternion space
% perform forward simulation
% check the frequency of intersection
% July 12, 2021

function [Output,SimuSpots,HittedSpots]=ForwardCalc_OR_space_comp(q_b,q_c,q_d,dim_OR,RotDet, ...
    proj_bin_bw,pos,rot_angles,S,B,Ahkl,nrhkl,hklnumber,hkl_square,thetamax,lambda_min,lambda_max, ...
    Lsam2sou,Lsam2det,minEucDis,dety00,detz00,P0y,P0z, ...
    pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ)

% % for testing
% q_b=OR.q_b;
% q_c=OR.q_c;
% q_d=OR.q_d;
% dim_OR=size(OR.q_b,1);
% pos=pos_indexing;

    L=Lsam2sou+Lsam2det; % [mm]
%     parfor i=1:dim(1)*dim(2)*dim(3)
% %         clear indices;
%         [ind1,ind2,ind3]=ind2sub(dim,i);
    if dim_OR>2000
        parfor i=1:dim_OR
            x = [q_b(i) q_c(i) q_d(i)];
            % x is represented by the last three quantities of the quaternions (a, b, c, d)
            % Extract the values from Q
            if (x(1)^2+x(2)^2+x(3)^2)<=1
                q0 = sqrt(1-(x(1)^2+x(2)^2+x(3)^2));
                q1 = x(1);
                q2 = x(2);
                q3 = x(3);

                % First row of the rotation matrix
                r00 = 2 * (q0 * q0 + q1 * q1) - 1;
                r01 = 2 * (q1 * q2 - q0 * q3);
                r02 = 2 * (q1 * q3 + q0 * q2);

                % Second row of the rotation matrix
                r10 = 2 * (q1 * q2 + q0 * q3);
                r11 = 2 * (q0 * q0 + q2 * q2) - 1;
                r12 = 2 * (q2 * q3 - q0 * q1);

                % Third row of the rotation matrix
                r20 = 2 * (q1 * q3 - q0 * q2);
                r21 = 2 * (q2 * q3 + q0 * q1);
                r22 = 2 * (q0 * q0 + q3 * q3) - 1;

                % 3x3 rotation matrix
                U = [r00 r01 r02;r10 r11 r12;r20 r21 r22];
                [Nr_simu,Nr_intersect,dis_median,SimuSpots{i},HittedSpots{i}]=calcForward(U,proj_bin_bw,pos,rot_angles,S,B,Ahkl,nrhkl,hklnumber,hkl_square, ...
                    RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,L,minEucDis,dety00,detz00,P0y,P0z, ...
                    pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ);
                Output(i,:)=[Nr_simu Nr_intersect Nr_intersect/Nr_simu dis_median q0 q1 q2 q3];
            end
            i;
        end
    else
        for i=1:dim_OR
        x = [q_b(i) q_c(i) q_d(i)];
        % x is represented by the last three quantities of the quaternions (a, b, c, d)
        % Extract the values from Q
        if (x(1)^2+x(2)^2+x(3)^2)<=1
            q0 = sqrt(1-(x(1)^2+x(2)^2+x(3)^2));
            q1 = x(1);
            q2 = x(2);
            q3 = x(3);

            % First row of the rotation matrix
            r00 = 2 * (q0 * q0 + q1 * q1) - 1;
            r01 = 2 * (q1 * q2 - q0 * q3);
            r02 = 2 * (q1 * q3 + q0 * q2);

            % Second row of the rotation matrix
            r10 = 2 * (q1 * q2 + q0 * q3);
            r11 = 2 * (q0 * q0 + q2 * q2) - 1;
            r12 = 2 * (q2 * q3 - q0 * q1);

            % Third row of the rotation matrix
            r20 = 2 * (q1 * q3 - q0 * q2);
            r21 = 2 * (q2 * q3 + q0 * q1);
            r22 = 2 * (q0 * q0 + q3 * q3) - 1;

            % 3x3 rotation matrix
            U = [r00 r01 r02;r10 r11 r12;r20 r21 r22];
            [Nr_simu,Nr_intersect,dis_median,SimuSpots{i},HittedSpots{i}]=calcForward(U,proj_bin_bw,pos,rot_angles,S,B,Ahkl,nrhkl,hklnumber,hkl_square, ...
                RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,L,minEucDis,dety00,detz00,P0y,P0z, ...
                pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ);
            Output(i,:)=[Nr_simu Nr_intersect Nr_intersect/Nr_simu dis_median q0 q1 q2 q3];
        end
        i;
        end
    end
end

function [Nr_simu,Nr_intersect,dis_median,SimuSpots,HittedSpots]=calcForward(U,proj_bin_bw,pos,rot_angles,S,B,Ahkl,nrhkl,hklnumber,hkl_square, ...
                RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,L,minEucDis,dety00,detz00,P0y,P0z, ...
                pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ)
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
    t = (RotDet(1,1)*(Lsam2det-SamposW(1))-RotDet(2,1)*(dety00-SamposW(2))-RotDet(3,1)*(detz00-SamposW(3)))./ ...
    (RotDet(1,1)*K_out_unit(1,:)+RotDet(2,1)*K_out_unit(2,:)+RotDet(3,1)*K_out_unit(3,:));

    dety22 = [RotDet(1,2) RotDet(2,2) RotDet(3,2)]*(t.*K_out_unit+[SamposW(1)-Lsam2det SamposW(2)-dety00 SamposW(3)-detz00]');
    detz22 = [RotDet(1,3) RotDet(2,3) RotDet(3,3)]*(t.*K_out_unit+[SamposW(1)-Lsam2det SamposW(2)-dety00 SamposW(3)-detz00]');
%     dety = -round(dety22/pixelysize-0.5)+dety0; % [pixel]
%     detz = -round(detz22/pixelzsize-0.5)+detz0; % [pixel]
    dety = -round(dety22/pixelysize)+dety0; % [pixel]
    detz = -round(detz22/pixelzsize)+detz0; % [pixel]
    
    select3=find(beta > pi/2 & beta < (90+thetamax*4)/180*pi & ...
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
                HittedSpots=[HittedSpots;hkl(:,select3(j))' Energy_hkl(select3(j)) ...
                    rot_angles(i) dety22(select3(j)) detz22(select3(j)) dis];
            end
            SimuSpots=[SimuSpots;hkl(:,select3(j))' Energy_hkl(select3(j)) rot_angles(i) ...
                dety22(select3(j)) detz22(select3(j)) dis];
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
end


% this is about 30% slower
%{
function [Nr_simu,Nr_intersect,dis_median,SimuSpots,HittedSpots]=calcForward(U,proj_bin_bw,pos,rot_angles,S,B,Ahkl,nrhkl,hklnumber,hkl_square, ...
                RotDet,thetamax,lambda_min,lambda_max,Lsam2sou,Lsam2det,L,minEucDis,dety00,detz00,P0y,P0z, ...
                pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ)
Nr_simu=0;
Nr_intersect=0;
dis_mean=0;
HittedSpots=[];
SimuSpots=[];

for i=1:length(rot_angles)    
    omega=rot_angles(i)*pi/180; % [rad]
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
            if beta > pi/2 %&& beta < (90+thetamax*4)/180*pi
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
                    
                    %%% include detector tilt and offset
                    K_out_unit = ([Lsam2det dety22 detz22]'-SamposW)./norm([Lsam2det dety22 detz22]'-SamposW);
                    t = (RotDet(1,1)*(Lsam2det-SamposW(1))-RotDet(2,1)*(dety00-SamposW(2))-RotDet(3,1)*(detz00-SamposW(3)))./ ...
                    (RotDet(1,1)*K_out_unit(1)+RotDet(2,1)*K_out_unit(2)+RotDet(3,1)*K_out_unit(3));
                    dety22 = [RotDet(1,2) RotDet(2,2) RotDet(3,2)]*(t*K_out_unit+[SamposW(1)-Lsam2det SamposW(2)-dety00 SamposW(3)-detz00]');
                    detz22 = [RotDet(1,3) RotDet(2,3) RotDet(3,3)]*(t*K_out_unit+[SamposW(1)-Lsam2det SamposW(2)-dety00 SamposW(3)-detz00]');
                    dety = -round(dety22/pixelysize-0.5)+dety0; % [pixel]
                    detz = -round(detz22/pixelzsize-0.5)+detz0; % [pixel]                    
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

dis_median=median(SimuSpots(:,8));
end
%}
