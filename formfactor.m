% Calculation of the atomic form factor at a specified sin(theta)/lambda
% using the analytic fit to the Direc form factors from Int. Tab. Cryst Sect. C
% INPUT : atom_no
%         atomlib
%         sin(theta)/lambda 
% OUTPUT: formfactor of atom with atom_no at sin(theta)/lambda
%

function formfac = formfactor(atomno,atom,stl)

formfac = 0.0;
for i=1:4
  formfac = formfac + atom(atomno).a(i)*exp(-atom(atomno).b(i)*(stl)^2);
end
formfac = formfac + atom(atomno).c;
