% calculate the beam path length along the sample
% incoming beam length + outcoming diffracting beam length
% June 20, 2020
% length: [mm]
% angle: [rad]
% diffraction occurring at (x,y,z)

function [A_Ehkl L_total]=beam_attenuation(SubGrain_posW,Lsam2sou,Lsam2det,dety22,detz22, ...
    atomparam,Transmission,rou,Energy_hkl,Rsample)

if nargin<10
    Rsample=0.35; % [mm]
end

% length of incoming beam path
if (Rsample^2*(Lsam2sou+SubGrain_posW(1))^2+ ...
    SubGrain_posW(2)^2*(Rsample^2-Lsam2sou^2))>=0
    t1=(Lsam2sou*(Lsam2sou+SubGrain_posW(1))-sqrt(Rsample^2*(Lsam2sou+SubGrain_posW(1))^2+ ...
        SubGrain_posW(2)^2*(Rsample^2-Lsam2sou^2)))/((Lsam2sou+SubGrain_posW(1))^2+SubGrain_posW(2)^2);
else
    t1=1-Rsample/Lsam2sou; % approximate solution is used when Rsample is not accurately estimated
end
xn=-Lsam2sou+t1*(Lsam2sou+SubGrain_posW(1));
yn=t1*SubGrain_posW(2);
zn=t1*SubGrain_posW(3);
L_NM=sqrt((xn-SubGrain_posW(1))^2+(yn-SubGrain_posW(2))^2+(zn-SubGrain_posW(3))^2);

% length of diffracted beam path
if (2*SubGrain_posW(1)*SubGrain_posW(2)*(Lsam2det-SubGrain_posW(1))*(dety22-SubGrain_posW(2))+ ...
    Rsample^2*((Lsam2det-SubGrain_posW(1))^2+(dety22-SubGrain_posW(2))^2)- ...
    SubGrain_posW(1)^2*(dety22-SubGrain_posW(2))^2-SubGrain_posW(2)^2*(Lsam2det-SubGrain_posW(1))^2)>=0
    t2=(-SubGrain_posW(1)*(Lsam2det-SubGrain_posW(1))-SubGrain_posW(2)*(dety22-SubGrain_posW(2))+ ...
        sqrt(2*SubGrain_posW(1)*SubGrain_posW(2)*(Lsam2det-SubGrain_posW(1))*(dety22-SubGrain_posW(2))+ ...
        Rsample^2*((Lsam2det-SubGrain_posW(1))^2+(dety22-SubGrain_posW(2))^2)- ...
        SubGrain_posW(1)^2*(dety22-SubGrain_posW(2))^2-SubGrain_posW(2)^2*(Lsam2det-SubGrain_posW(1))^2)) ...
        /((Lsam2det-SubGrain_posW(1))^2+(dety22-SubGrain_posW(2))^2);
else
    t2=(-SubGrain_posW(1)+sqrt(Rsample^2-SubGrain_posW(2)^2))/Lsam2det; % approximate solution is used when Rsample is not accurately estimated
end
xq1=SubGrain_posW(1)+t2*(Lsam2det-SubGrain_posW(1));
yq1=SubGrain_posW(2)+t2*(dety22-SubGrain_posW(2));
zq1=SubGrain_posW(3)+t2*(detz22-SubGrain_posW(3));
L_MQ1=sqrt((xq1-SubGrain_posW(1))^2+(yq1-SubGrain_posW(2))^2+(zq1-SubGrain_posW(3))^2);

L_total=L_NM+L_MQ1;

% % The following is wrong!!!!!!!!!
% % length of incoming beam path
% L_SO1=sqrt(Lsam2sou^2+SubGrain_posW(3)^2);
% L_SM=sqrt((Lsam2sou+SubGrain_posW(1))^2+SubGrain_posW(2)^2+SubGrain_posW(3)^2);
% L_SN=L_SO1*cos(alpha)-sqrt(L_SO1^2*(cos(alpha)^2-1)+Rsample^2); % cosine law in triangle SNO1
% L_NM=L_SM-L_SN; % length of incoming beam path
% 
% % length of outcoming diffracting beam path
% L_SP1=L_SO1*cos(alpha)+sqrt(L_SO1^2*(cos(alpha)^2-1)+Rsample^2);
% L_MP1=L_SP1-L_SM;
% L_MO1=sqrt(SubGrain_posW(1)^2+SubGrain_posW(2)^2);
% cos_theta1=(L_MP1^2+L_MO1^2-Rsample^2)/(2*L_MP1*L_MO1); % cosine law in triangle O1MP1
% theta1=acos(cos_theta1);
% L_MQ1=L_MO1*cos(theta1+2*theta)+sqrt(L_MO1^2*(cos(theta1+2*theta)^2-1)+Rsample^2);% length of outcoming diffracting beam path
% L_total=L_NM+L_MQ1;
% L_MO1*cos(theta1-2*theta)+sqrt(L_MO1^2*(cos(theta1-2*theta)^2-1)+Rsample^2)
    
for j=1:length(Transmission)-1
    if Transmission(j,1)<=Energy_hkl/1e3 && Transmission(j+1,1)>=Energy_hkl/1e3
        Eindex=[j j+1];
    end
end
A_Ehkl=exp(-(Transmission(Eindex(1),2)+(Energy_hkl/1e3-(Transmission(Eindex(1),1))) ...
    *(Transmission(Eindex(2),2)-Transmission(Eindex(1),2))/ ...
    (Transmission(Eindex(2),1)-Transmission(Eindex(1),1)))*L_total*1e-1*rou);










