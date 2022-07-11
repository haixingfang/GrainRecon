% DQE: detective quantum efficiency
% quantum efficiency assumed to be 1
% scientillation efficiency assumed to be proportional to the absorption by
% the scintillator CsI
% data from https://physics.nist.gov/PhysRefData/Xcom/html/xcom1.html
% June 23, 2020

function [DQE_Ehkl]=Detector_efficiency(CsI,Swank,Energy_hkl)

for j=1:length(CsI(:,1))-1
    if CsI(j,1)<=Energy_hkl/1e3 && CsI(j+1,1)>=Energy_hkl/1e3
        Eindex=[j j+1];
    elseif CsI(1,1)>Energy_hkl/1e3
        Eindex=[1 1+1];
    end
end
thickness=150; % thickness of scintillator [um]
density=4.51; % density of CsI [g/cm^3]
aqe=1-exp(-(CsI(Eindex(1),8)+(Energy_hkl/1e3-(CsI(Eindex(1),1))) ...
    *(CsI(Eindex(2),8)-CsI(Eindex(1),8))/ ...
    (CsI(Eindex(2),1)-CsI(Eindex(1),1)))*thickness*1e-4*density); % absorption quantum efficiency

% Swank statistical factor at zero spatial frequency
% data from Hajdok, G., Battista, J.J. & Cunningham, I.A (2008). Med. Phys. 35, 3194-3204.
for j=1:length(Swank(:,1))-1
    if Swank(j,1)<=Energy_hkl && Swank(j+1,1)>=Energy_hkl
        Eindex=[j j+1];
    end
end
Swank_factor=(Swank(Eindex(1),2)+(Energy_hkl-(Swank(Eindex(1),1))) ...
    *(Swank(Eindex(2),2)-Swank(Eindex(1),2))/ ...
    (Swank(Eindex(2),1)-Swank(Eindex(1),1)));

DQE_Ehkl=aqe*Swank_factor;

%{
%%%% for testing
Energy0=5:1:140;
thickness=[50 100 150];
for i=1:length(thickness)
    for k=1:length(Energy0)
        Energy_hkl=Energy0(k);
        for j=1:length(CsI(:,1))-1
            if CsI(j,1)<=Energy_hkl/1e3 && CsI(j+1,1)>=Energy_hkl/1e3
                Eindex=[j j+1];
            end
        end
        density=4.51; % density of CsI [g/cm^3]
        aqe=1-exp(-(CsI(Eindex(1),8)+(Energy_hkl/1e3-(CsI(Eindex(1),1))) ...
            *(CsI(Eindex(2),8)-CsI(Eindex(1),8))/ ...
            (CsI(Eindex(2),1)-CsI(Eindex(1),1)))*thickness(i)*1e-4*density);
        
        % Swank statistical factor at zero spatial frequency
        % data from Hajdok, G., Battista, J.J. & Cunningham, I.A (2008). Med. Phys. 35, 3194-3204.
        for j=1:length(Swank(:,1))-1
            if Swank(j,1)<=Energy_hkl && Swank(j+1,1)>=Energy_hkl
                Eindex=[j j+1];
            end
        end
        Swank_factor=(Swank(Eindex(1),2)+(Energy_hkl-(Swank(Eindex(1),1))) ...
            *(Swank(Eindex(2),2)-Swank(Eindex(1),2))/ ...
            (Swank(Eindex(2),1)-Swank(Eindex(1),1)));        

        DQE_Ehkl=aqe*Swank_factor;
        DQE(i,k)=DQE_Ehkl;
    end
end
figure;
hold all;
plot(Energy0,DQE(1,:),'ro-');
plot(Energy0,DQE(2,:),'bx-');
plot(Energy0,DQE(3,:),'m^-');
% legend('t = 25 \mum','t = 50 \mum','t = 100 \mum');
legend('t = 50 \mum','t = 100 \mum','t = 150 \mum');
dlmwrite('DQE.txt',[Energy0' DQE'],'delimiter',' ');
%}




