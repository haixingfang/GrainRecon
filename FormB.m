% calculate B matrix of (Gcart = B Ghkl) following eq. 3.4 in 
%   H.F. Poulsen.
%   Three-dimensional X-ray diffraction microscopy. 
%   Mapping polycrystals and their dynamics. 
%   (Springer Tracts in Modern Physics, v. 205), (Springer, Berlin, 2004).
%
%
% B = FormB(unit_cell)
%
% INPUT : unit_cell = [a b c alpha beta gamma] 
% OUTPUT: B [3x3]
%
% Henning Osholm Sørensen, June 23, 2006.

function B = FormB(ucell)

     a   = ucell(1);
     b   = ucell(2);
     c   = ucell(3);
     calp = cos(ucell(4)*pi/180.);
     cbet = cos(ucell(5)*pi/180.);
     cgam = cos(ucell(6)*pi/180.);
     salp = sin(ucell(4)*pi/180.);
     sbet = sin(ucell(5)*pi/180.);
     sgam = sin(ucell(6)*pi/180.);

  angular = sqrt(1 - calp^2 - cbet^2 - cgam^2 + 2*calp*cbet*cgam);
  V = a*b*c*angular ;                                            %Volume of unit cell


  %  Calculate reciprocal lattice parameters: NOTICE PHYSICIST DEFINITION of recip axes with 2*pi
  astar = 2*pi*b*c*salp/V                     ;   
  bstar = 2*pi*a*c*sbet/V                     ;   
  cstar = 2*pi*a*b*sgam/V                     ;   
  salpstar = V/(a*b*c*sbet*sgam)              ;   
  sbetstar = V/(a*b*c*salp*sgam)              ;   
  sgamstar = V/(a*b*c*salp*sbet)              ;   
  calpstar = (cbet*cgam-calp)/(sbet*sgam)     ;   
  cbetstar = (calp*cgam-cbet)/(salp*sgam)     ;   
  cgamstar = (calp*cbet-cgam)/(salp*sbet)     ;   

% Form B matrix following eq. 3.4 in H.F Poulsen
  B(1,1) = astar                              ;
  B(1,2) = bstar*cgamstar                     ;
  B(1,3) = cstar*cbetstar                     ;
  B(2,1) = 0                                  ;
  B(2,2) = bstar*sgamstar                     ;
  B(2,3) = -cstar*sbetstar*calp               ;
  B(3,1) = 0                                  ; 
  B(3,2) = 0                                  ;
  B(3,3) = cstar*sbetstar*salp                ;



