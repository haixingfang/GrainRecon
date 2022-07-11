% This function generates X-ray spectrum for W anode material
% based on John M. Boonea and J. Anthony Seibert, 1997, Am. Assoc. Med. Phys.
% Empirical fitting

function [Ebin, I0E]=Xray_Spectrum_W(I0,flag)
% I0 = 1e8;       % Beam flux from a lab source 10^8 (photons/s/mm2)
E=[12.50 14.900222 16.762749 19.24612 20.798225 22.039911 23.90244 24.833702 ...
    26.385809 28.248337 30.421286 33.215076 37.871395 42.527718 47.80488 ...
    51.529934 54.013306 55.56541 56.807095 59.600887 60.84257 ...
    60.94257 61.152992 61.773834 62.7051 65.49889 66.14058 66.43015 67.36142 68.6031 ...
    68.91353 69.84479 70.776054 73.25942 78.53658 83.81374 90.64301 99.02439 105.543236 ...
    114.23503 120.44346 128.204 131.92905 138.13747]; % [keV]

PhotonsPerE=[0.00015 0.000414907 0.001659629 0.003319257 0.004875159 0.006742242 0.008505597 0.010061499 ...
    0.011824854 0.013380756 0.014936658 0.016492561 0.017011194 0.016492561 0.015870198 ...
    0.014729205 0.014106843 0.013691936 0.021575173 0.036096922 0.021990081 0.018670823 ...
    0.014832932 0.012965849 0.011824854 ...
    0.010580133 0.012758396 0.014936658 0.017737282 0.015870198 0.012654669 0.010165226 ...
    0.00840187 0.006638514 0.006119881 0.00518634 0.004252798 0.003526711 0.002696897 ...
    0.002178263 0.001555902 0.001244722 0.000829814 0.000207454]; % photon flux per energy bin [photons/s/mm^2/keV]
Ebin = linspace(13.0,138,50);
% Ebin = linspace(14.9,138,50);
cs = spline(E,PhotonsPerE,Ebin);
I0E=cs*I0*(Ebin(2)-Ebin(1)); % photon flux for a specific energy [photons/s/mm^2]

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





