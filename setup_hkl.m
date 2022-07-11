% set up X-ray spectra, range of wave vectors and Miller indices hkl for indexing

% basic information of X-ray spectra and range of wave vector
mono_beam=0; % monochoromatic photon beam or not, 0-polychromatic; 1-monochromatic. By default it is polychromatic
% I0 = 1e14; % Beam flux from a synchrotron source 10^14-15 (photons/s/mm2)
I0 = 5e8;    % Beam flux from a lab source 10^8 (photons/s/mm2)
if mono_beam==1
    Energy = 59.31;  % photon energy [keV], K alpha1 line of W anode material
else
    if simap_data_flag~=1
        [Energy, I0E]=Xray_Spectrum_W(I0,0);
    else
        [Energy, I0E]=Xray_Spectrum_W_simap(I0,0);
    end
    I0E=abs(I0E);
end
tthetamax = acos(dot([Lsam2det 0 0]./norm([Lsam2det 0 0]), ...
    [Lsam2det 0.5*detysize*pixelysize 0.5*detzsize*pixelzsize]./norm([Lsam2det 0.5*detysize*pixelysize 0.5*detzsize*pixelzsize])));
tthetamax = tthetamax*180/pi; % two-theta [deg]
thetamax=tthetamax/2;
lambda_max = 12.39818746/min(Energy);   % [Angstrom]
lambda_min = 12.39818746/max(Energy); % [Angstrom]
Kmax = 1/lambda_min;
Kmin = 1/lambda_max;
lambda = 12.398./Energy;
% sintlmax = sin(thetamax*pi/180)/(12.398/Energy(find(I0E==max(I0E)))); % sin(theta)/lambda [A^-1], consider the characteristic wavelength which corresponds to highest flux
sintlmax = sin(thetamax*pi/180)/(12.398/55);
Ki_max = [-2*pi/lambda_min 0 0];
Klen_max = -Ki_max(1);
Ki_min = [-2*pi/lambda_max 0 0];
Klen_min = -Ki_min(1);
emass =9.1093826e-31;
echarge = 1.60217653e-19;
pi4eps0 = 1.11265e-10;
c = 299792458 ;
K1 =  (echarge^2/(pi4eps0*emass*c*c)*1000)^2; % square of Thomson scattering length r0^2, [mm^2]

sg = sglib(space_group_IT_number); % get sysconditions for specific element from the sglib.m
sysconditions=sg.sysconditions;

readhkl = 0;
structfact = 1; % do calculate structure factors
if readhkl == 0
    % Generate Miller indices for reflections within a certain resolution  
    % only compute the first several hkl families
    Ahkl0  = genhkl(cell,sysconditions,1.5*sintlmax);
    hkl_square=Ahkl0(:,1).^2+Ahkl0(:,2).^2+Ahkl0(:,3).^2;
    hkl_square=sortrows(unique(hkl_square));
    Ahkl=[];
%     hklnumber=3; % maximum is 10, recommended be at least >= 3, 4
    if length(hkl_square(:,1))>=hklnumber
        hkl2_max=hkl_square(hklnumber);
    else
        hkl2_max=hkl_square(end);
    end
    for i=1:length(Ahkl0(:,1))
        if (Ahkl0(i,1).^2+Ahkl0(i,2).^2+Ahkl0(i,3).^2)<=hkl2_max
            Ahkl=[Ahkl;Ahkl0(i,:)];
        end
    end
      
    % Initialize Ahkl
    nrhkl = size(Ahkl,1);
    if structfact == 1
        disp('Calculating Structure factors');
        atomlib;
        hkl = [0 0 0];
        for i=1:nrhkl
            hkl2 = [Ahkl(i,1) Ahkl(i,2) Ahkl(i,3)];            
            if all(hkl2 == -1*hkl) %Only calculate F^2 if not Friedel mate 
                hkl = hkl2;
            else
                hkl = hkl2;
                [Freal, Fimg] = structure_factor(hkl,cell,atomparam,sg,atom);
                int = Freal^2 + Fimg^2;
            end
            Ahkl(i,5) = int;
        end
        %disp('Finished Calculating Structure factors'); disp(' ')
        
        
%         %%%%%%%%%%%%%% select the most strongest hkl for indexing June 28, 2021
%         Ahkl=sortrows(Ahkl,5,'descend');
%         Ahkl(:,6)=Ahkl(:,1).^2+Ahkl(:,2).^2+Ahkl(:,3).^2;
% %         figure;
% %         plot(Ahkl(:,6),Ahkl(:,5),'ro','MarkerSize',8,'LineWidth',1.5);
% %         xlabel('h^{2}+k^{2}+l^{2}');
% %         ylabel('Structure factor');
% %         set(gca,'LineWidth',1.5);
% %         set(gca,'FontSize',18);
%         hklnumber=3;
%         hkl_square_temp=[];
%         Ahkl_temp=[];
%         count=0;
%         for i=1:length(Ahkl(:,1))
%             if isempty(find(hkl_square_temp==Ahkl(i,6)))
%                 count=count+1;
%                 hkl_square_temp=[hkl_square_temp;Ahkl(i,6)];
%                 Ahkl_temp=[Ahkl_temp;Ahkl(find(Ahkl(:,6)==Ahkl(i,6)),:)];
%             end
%             if count==hklnumber
%                 break;
%             end
%         end
%         hkl_square=hkl_square_temp;
%         Ahkl=Ahkl_temp;
%         nrhkl = size(Ahkl,1);

%         %%% specify hkl for indexing
%         hkl_square_temp=[3 8 11]; % specify the hkls
%         hklnumber=length(hkl_square_temp);
%         Ahkl(:,6)=Ahkl(:,1).^2+Ahkl(:,2).^2+Ahkl(:,3).^2;
%         Ahkl_temp=[];
%         for i=1:length(hkl_square_temp)
%             Ahkl_temp=[Ahkl_temp;Ahkl(find(Ahkl(:,6)==hkl_square_temp(i)),:)];          
%         end
%         hkl_square=hkl_square_temp;
%         Ahkl=Ahkl_temp;
%         nrhkl = size(Ahkl,1);
    else
        for i=1:nrhkl
            Ahkl(i,5) = 32768; % half of 2^16
        end
    end
%     save('Ahkl_Al_fcc_4hkl.mat','Ahkl','hklnumber','nrhkl','hkl_square');
else
    disp('Please load the pre-calculated list of hkl reflections *.mat from folder');
    filename = uigetfile({'*.mat'}, 'Select mat file'); 
    load(filename)
end

% add on Sep 24, 2021
hkl_family=unique(sort(abs(Ahkl(:,1:3)),2),'rows'); % family of hkl indices for indexing
hkl_family(:,4)=hkl_family(:,1).^2+hkl_family(:,2).^2+hkl_family(:,3).^2;
hkl_family=sortrows(hkl_family,4);
hkl_family=hkl_family(:,1:3);
for i=1:length(hkl_family(:,1))
    hkl_family_square(i)=sum(hkl_family(i,:).^2);
end
d_possible=sqrt(1./(hkl_family(:,1).^2./cell(1)^2+ ...
    hkl_family(:,2).^2./cell(2)^2+hkl_family(:,3).^2./cell(3)^2)); % general expression for Orthorhombic [A]
Glen_possible=2*pi./d_possible;

