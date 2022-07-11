% This function generates X-ray spectrum for W anode material
% X-ray source at SIMaP
% Nano_80kv_Large_11uA_calbration_spectrum

function [Ebin, I0E]=Xray_Spectrum_W_simap(I0,flag)
I0 = 1e6;       % Beam flux from a lab source 10^8 (photons/s/mm2)

Xray=textread('Nano_80kv_Large_11uA_calbration_spectrum.txt');
E=Xray(:,1); % [keV]
PhotonsPerE=Xray(:,2)/60; % photon flux per energy bin [photons/s/mm^2/keV]

Ebin = linspace(2,145,400);
cs = spline(E,PhotonsPerE,Ebin);
I0E=cs*I0*(Ebin(2)-Ebin(1)); % photon flux for a specific energy [photons/s/mm^2]
I0E=abs(I0E);

if nargin<=1
    figure('Name','X-ray spectra');
    subplot(1,2,1);
    plot(E,PhotonsPerE,'bo',Ebin,cs,'r-');
    xlabel('Energy (keV)');
    ylabel('Photons per energy bin');

    subplot(1,2,2);
    plot(Ebin,I0E,'ro-');
    xlabel('Energy (keV)');
    ylabel('Photon flux (photons/s/mm^{2})');
end

