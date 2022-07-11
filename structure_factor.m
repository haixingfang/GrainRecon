% Calculation of the structure factor of reflection hkl
%
%     [Freal Fimg] = structure_factor(hkl,unit_cell,atomparam,sg,atom)
%
% INPUT : hkl =       [h k l] 
%         unit_cell = [a b c alpha beta gamma] 
%         atomparam:  structural parameters
%         sg:         read from sglib(sgno) (space group library)
%         atom:       read from atomlib (library of form factors)
% OUTPUT: The real and imaginary parts of the the structure factor
%

function [Freal Fimg] = structure_factor(hkl,cell,atomparam,sg,atom)

stl = sintl(cell,hkl);
noatoms = size(atomparam,2);
Freal = 0.0;
Fimg = 0.0;


for i = 1:noatoms
    sfreal = 0.0;
    sfimg = 0.0;
    
    %Check whether isotrop or anisotropic displacements 
    Uelements = size(atomparam(i).adp,2);

    if Uelements == 1
        U = atomparam(i).adp;
    elseif  Uelements == 6
        % transform Uij to betaij
        betaij = Uij2betaij(atomparam(i).adp,cell);
    else
        error(['wrong no of elements in atomparam(',num2str(i),').adp'])
    end

    

    for j = 1:sg.nosymop
        % atomic displacement factor
        if Uelements == 1
            expij=exp(-2*pi^2*U*stl^2);
        else
            betaijrot = sg.symop(j).rot*betaij*sg.symop(j).rot';
            expij=exp(-hkl*betaijrot*hkl');
        end
        % exponent for phase factor
        r = sg.symop(j).rot*atomparam(i).pos' + sg.symop(j).trans';
        exponent = 2*pi*hkl*r;

        %forming the real and imaginary parts of F
        sfreal = sfreal + expij*cos(exponent);
        sfimg = sfimg + expij*sin(exponent);
    end

    % Including the atomic formfactor
    formfac = formfactor(atomparam(i).atomno,atom,stl)*atomparam(i).occ/atomparam(i).symmulti ;

    Freal = Freal + formfac*sfreal ;
    Fimg = Fimg + formfac*sfimg ;
end

