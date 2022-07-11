% discretize the quaternion space
% perform forward simulation
% check the matched number of G vectors
% Sep 24, 2021

function [Output,Gt_matched_all]=ForwardCalc_OR_space_Gt(q_b,q_c,q_d,dim_OR,RotDet, ...
    Spots,pos,rot_angles,rot_start,rot_step,S,B,Ahkl,nrhkl,hkl_family,hkl_family_square,d_possible,Glen_possible, ...
    Lsam2sou,Lsam2det,dety00,detz00,P0y,P0z,pixelysize,pixelzsize,dety0,detz0, ...
    thetamax,lambda_min,lambda_max,detysize,detzsize,BeamStopY,BeamStopZ)

% % for testing
% q_b=OR.q_b;
% q_c=OR.q_c;
% q_d=OR.q_d;
% dim_OR=size(OR.q_b,1);
% pos=pos_indexing;

    L=Lsam2sou+Lsam2det; % [mm]
    Output=zeros(dim_OR,8);
    Nr_angles=length(rot_angles);
%     Nr_angles=90;
    if length(rot_angles)>Nr_angles
        rot_angles_select=rot_angles(randperm(length(rot_angles))); % make a random shuffle
        rot_angles_select=rot_angles_select(1:Nr_angles);
    else
        rot_angles_select=rot_angles;
    end
    
    % Jan 6, 2021
    % select spots which are not located close to rotation axis
    SpotsSelect_flag=0;
    if SpotsSelect_flag==1
         for i=1:length(Spots)
            SpotsSelect{i}=Spots{i}(find((Spots{i}(:,4)>10 & Spots{i}(:,4)<170) | (Spots{i}(:,4)>190 & Spots{i}(:,4)<350)),:);
        end
    else
        SpotsSelect=Spots;
    end

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
                [Gt_matched_all{i},Nr_match,~,Gt_matched_cost,NrSpotExpectAll]=calc_Gt_match(U,pos,SpotsSelect, ...
                    rot_angles_select,rot_start,rot_step,S,B,Ahkl,nrhkl,RotDet, ...
                    hkl_family,hkl_family_square,d_possible,Glen_possible,Lsam2sou,Lsam2det,dety00,detz00,P0y,P0z, ...
                    pixelysize,pixelzsize,dety0,detz0,thetamax,lambda_min,lambda_max,detysize,detzsize,BeamStopY,BeamStopZ);
                Output(i,:)=[i Nr_match Gt_matched_cost Nr_match/NrSpotExpectAll q0 q1 q2 q3];
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
            [Gt_matched_all{i},Nr_match,~,Gt_matched_cost,NrSpotExpectAll]=calc_Gt_match(U,pos,SpotsSelect, ...
                rot_angles_select,rot_start,rot_step,S,B,Ahkl,nrhkl,RotDet, ...
                hkl_family,hkl_family_square,d_possible,Glen_possible,Lsam2sou,Lsam2det,dety00,detz00,P0y,P0z, ...
                pixelysize,pixelzsize,dety0,detz0,thetamax,lambda_min,lambda_max,detysize,detzsize,BeamStopY,BeamStopZ);            
            Output(i,:)=[i Nr_match Gt_matched_cost Nr_match/NrSpotExpectAll q0 q1 q2 q3];
        end
        i;
        end
    end
end

% calculate the matched cases of G vectors
% Sep 23, 2021
function [Gt_matched_all,Nr_match,Gt_matched_all_mean,Gt_matched_cost,NrSpotExpectAll]=calc_Gt_match(U,pos,Spots, ...
    rot_angles,rot_start,rot_step,S,B,Ahkl,nrhkl,RotDet, ...
    hkl_family,hkl_family_square,d_possible,Glen_possible,Lsam2sou,Lsam2det,dety00,detz00,P0y,P0z, ...
    pixelysize,pixelzsize,dety0,detz0,thetamax,lambda_min,lambda_max,detysize,detzsize,BeamStopY,BeamStopZ)
    
    L=Lsam2sou+Lsam2det;
    hkl = [Ahkl(1:nrhkl,1) Ahkl(1:nrhkl,2) Ahkl(1:nrhkl,3)]';
    Gw = S*U*B*hkl;
    Gt_matched_all=[];
    NrSpotExpectAll=0;
    
    % consider the detector tilts for calculating the diffraction vector
    d_tr=RotDet'*[Lsam2det dety00 detz00]';
    Lsam2det_tr=d_tr(1)/RotDet(1,1);
    dety00_tr=RotDet(1,2)*Lsam2det_tr-d_tr(2);
    detz00_tr=RotDet(1,3)*Lsam2det_tr-d_tr(3);
    L_tr=Lsam2sou+Lsam2det_tr;
    
    for rot=rot_angles
        rot_number=(rot-rot_start)/rot_step+1;
        dety=Spots{rot_number}(:,9);
        detz=Spots{rot_number}(:,10);

        omega=rot*pi/180; % [rad]
        Omega=[cos(omega) -sin(omega) 0;sin(omega) cos(omega) 0;0 0 1];
        SamposW=Omega*S*pos';
        center = [L_tr, (SamposW(2)-P0y)*L_tr/(Lsam2sou+SamposW(1)), ...
            (SamposW(3)-P0z)*L_tr/(Lsam2sou+SamposW(1))]; % sample center projected to the position of the detector
        alpha = atan(sqrt((SamposW(2)-P0y)^2+(SamposW(3)-P0z)^2)/(Lsam2sou+SamposW(1)));
        grainpos = [Lsam2sou+SamposW(1) SamposW(2)-P0y SamposW(3)-P0z];

        % unit vectors along incoming and diffracted beams
        Kin_unit=[Lsam2sou+SamposW(1),SamposW(2)-P0y,SamposW(3)-P0z]./ ...
            norm([Lsam2sou+SamposW(1),SamposW(2)-P0y,SamposW(3)-P0z]); % unit vector along the incoming beam 
%         dety22=(dety0-dety+0.5)*pixelysize; % [mm]
%         detz22=(detz0-detz+0.5)*pixelzsize; % [mm]
        dety22=(dety0-dety)*pixelysize; % [mm]
        detz22=(detz0-detz)*pixelzsize; % [mm]
        
        Kout_unit=[repmat(Lsam2det-SamposW(1),length(dety),1),dety22-SamposW(2)-dety00,detz22-SamposW(3)-detz00]./ ...
            sqrt((Lsam2det-SamposW(1)).^2+(dety22-SamposW(2)-dety00).^2+(detz22-SamposW(3)-detz00).^2); % unit vector along the diffracted beam
        
        % consider the detector tilts for calculating the diffraction vector
%         d_tr=RotDet'*[Lsam2det dety00 detz00]';
%         Lsam2det_tr=d_tr(1)/RotDet(1,1);
%         dety00_tr=RotDet(1,2)*Lsam2det_tr-d_tr(2);
%         detz00_tr=RotDet(1,3)*Lsam2det_tr-d_tr(3);
        dety22_tr=dety22-dety00_tr;
        detz22_tr=detz22-detz00_tr;
        Kout_unit_tr=[repmat(Lsam2det_tr-SamposW(1),length(dety),1),dety22_tr-SamposW(2),detz22_tr-SamposW(3)]./ ...
            sqrt((Lsam2det_tr-SamposW(1)).^2+(dety22_tr-SamposW(2)).^2+(detz22_tr-SamposW(3)).^2);        

        ttheta_spot=acos(sum(repmat(Kin_unit,length(Kout_unit_tr(:,1)),1).*Kout_unit_tr,2))*180/pi; % 2-theta [deg]
        lambdahkl_possible = 2 * d_possible *sind(ttheta_spot'/2); % [A]
        Gt_spot_factor=Kout_unit_tr-Kin_unit;  % [mm]
    %     Gt_spot{rot_number}=[];
    %     for j=1:length(hkl_family(:,1))
    %         Gt_spot{rot_number}(:,:,j)=(2*pi./lambdahkl_possible(j,:)').*Gt_spot_factor;
    %     end
        for j=1:length(hkl_family(:,1))
            Gt_spot{j}=(2*pi./lambdahkl_possible(j,:)').*Gt_spot_factor;
        end    
        Gt=Omega*Gw;
        Gt=Gt';
        NrSpotExpect=calcExpectSpotNr(Gt',hkl,grainpos,center,alpha,RotDet, ...
        SamposW,thetamax,lambda_min,lambda_max,Lsam2det,dety00,detz00, ...
        pixelysize,pixelzsize,dety0,detz0,detysize,detzsize,BeamStopY,BeamStopZ);
        NrSpotExpectAll=NrSpotExpectAll+NrSpotExpect;

        % find matched diffraction vectors
        [~,Gt_matched]=find_Gt_match(Gt_spot,Gt,hkl,hkl_family_square,Glen_possible,rot,dety,detz);
        Gt_matched_all=[Gt_matched_all;Gt_matched];
    end
    if ~isempty(Gt_matched_all)
        Nr_match=length(Gt_matched_all(:,1));
        Gt_matched_all_mean=mean(Gt_matched_all(:,7));
        Gt_matched_all_median=median(Gt_matched_all(:,7));
        Gt_matched_cost=sum((Gt_matched_all(:,1)-Gt_matched_all(:,8)).^2)+ ...
            sum((Gt_matched_all(:,2)-Gt_matched_all(:,9)).^2)+ ...
            sum((Gt_matched_all(:,3)-Gt_matched_all(:,10)).^2);
        Gt_matched_cost=Gt_matched_cost/length(Gt_matched_all(:,1));
    else
        Nr_match=0;
        Gt_matched_all_mean=NaN;
        Gt_matched_all_median=NaN;
        Gt_matched_cost=100;
    end
end


function [Gt_match,Gt_matched]=find_Gt_match(Gt_spot,Gt,hkl,hkl_family_square,Glen_possible,rot,dety,detz)
% Gt_spot: list of Glab vectors calculated from the spot peaks;
% Gt: list of Glab vectors calculated from the forward simulation given a specific orientation
% Sep 24, 2021

for i=1:length(Gt_spot)
    Gt_spot_norm{i}=Gt_spot{i}./sqrt((Gt_spot{i}(:,1)).^2+(Gt_spot{i}(:,2)).^2+(Gt_spot{i}(:,3)).^2);
end
Gt_match=zeros(length(Gt(:,1)),17);
Gt_match(:,1:3)=Gt./sqrt(Gt(:,1).^2+Gt(:,2).^2+Gt(:,3).^2);
Gt_match(:,4)=sqrt(Gt(:,1).^2+Gt(:,2).^2+Gt(:,3).^2);
% angle_threshold=2.7; % minimum angle to accept as correct match, should be consistent with the mis-angle of discretized OR
angle_threshold=1; % Nov 9, 2021
for i=1:length(Gt_match(:,1))
    switch sum(hkl(:,i).^2)
        case hkl_family_square(1)
            ang=acos(sum(repmat(Gt_match(i,1:3),length(Gt_spot_norm{1}(:,1)),1).*Gt_spot_norm{1},2))*180/pi; % angle betwee [deg]
            [ang_min,MatchSpotID]=min(ang);
            if ang_min<angle_threshold && abs(Gt_match(i,4)-Glen_possible(1))<0.1
                Gt_match(i,5)=1; % find a matched experimental diffraction vector
                Gt_match(i,6)=MatchSpotID; % matched spot ID
                Gt_match(i,7)=ang_min;
                Gt_match(i,8:10)=Gt_spot_norm{1}(MatchSpotID,:); % Gt_spot
                Gt_match(i,11)=Glen_possible(1);  % Gt_spot length
                Gt_match(i,12:14)=hkl(:,i)';  % (hkl)
                Gt_match(i,15)=rot;  % [deg]
                Gt_match(i,16:17)=[dety(MatchSpotID) detz(MatchSpotID)];  % [dety, detz] [pixel]
            end
        case hkl_family_square(2)
            ang=acos(sum(repmat(Gt_match(i,1:3),length(Gt_spot_norm{2}(:,1)),1).*Gt_spot_norm{2},2))*180/pi;
            [ang_min,MatchSpotID]=min(ang);
            if ang_min<angle_threshold && abs(Gt_match(i,4)-Glen_possible(2))<0.1
                Gt_match(i,5)=1; % find a matched experimental diffraction vector
                Gt_match(i,6)=MatchSpotID; % matched spot ID
                Gt_match(i,7)=ang_min;
                Gt_match(i,8:10)=Gt_spot_norm{2}(MatchSpotID,:); % Gt_spot
                Gt_match(i,11)=Glen_possible(2);  % Gt_spot length
                Gt_match(i,12:14)=hkl(:,i)';  % (hkl)
                Gt_match(i,15)=rot;  % [deg]
                Gt_match(i,16:17)=[dety(MatchSpotID) detz(MatchSpotID)];  % [dety, detz] [pixel]
            end
        case hkl_family_square(3)
            ang=acos(sum(repmat(Gt_match(i,1:3),length(Gt_spot_norm{3}(:,1)),1).*Gt_spot_norm{3},2))*180/pi;
            [ang_min,MatchSpotID]=min(ang);
            if ang_min<angle_threshold && abs(Gt_match(i,4)-Glen_possible(3))<0.1
                Gt_match(i,5)=1; % find a matched experimental diffraction vector
                Gt_match(i,6)=MatchSpotID; % matched spot ID
                Gt_match(i,7)=ang_min;
                Gt_match(i,8:10)=Gt_spot_norm{3}(MatchSpotID,:); % Gt_spot
                Gt_match(i,11)=Glen_possible(3);  % Gt_spot length
                Gt_match(i,12:14)=hkl(:,i)';  % (hkl)
                Gt_match(i,15)=rot;  % [deg]
                Gt_match(i,16:17)=[dety(MatchSpotID) detz(MatchSpotID)];  % [dety, detz] [pixel]
            end
        case hkl_family_square(4)
            ang=acos(sum(repmat(Gt_match(i,1:3),length(Gt_spot_norm{4}(:,1)),1).*Gt_spot_norm{4},2))*180/pi;
            [ang_min,MatchSpotID]=min(ang);
            if ang_min<angle_threshold && abs(Gt_match(i,4)-Glen_possible(4))<0.1
                Gt_match(i,5)=1; % find a matched experimental diffraction vector
                Gt_match(i,6)=MatchSpotID; % matched spot ID
                Gt_match(i,7)=ang_min;
                Gt_match(i,8:10)=Gt_spot_norm{4}(MatchSpotID,:); % Gt_spot
                Gt_match(i,11)=Glen_possible(4);  % Gt_spot length
                Gt_match(i,12:14)=hkl(:,i)';  % (hkl)
                Gt_match(i,15)=rot;  % [deg]
                Gt_match(i,16:17)=[dety(MatchSpotID) detz(MatchSpotID)];  % [dety, detz] [pixel]
            end
    end
end
Gt_matched=Gt_match(Gt_match(:,5)==1,:);
end
    


