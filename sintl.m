% sintl calculate sin(theta)/lambda of the reflection "H" given the unit cell "cell" 
%
%    sintlc = sintl(unit_cell,hkl)
%
% INPUT:  unit_cell = [a b c alpha beta gamma]
%         hkl = [h k l]
% OUTPUT: sin(theta)/lambda

function sintlc = sintl(ucell,HKL)

     a   = ucell(1);
     b   = ucell(2);
     c   = ucell(3);
     calp = cos(ucell(4)*pi/180.);
     cbet = cos(ucell(5)*pi/180.);
     cgam = cos(ucell(6)*pi/180.);

     h = HKL(1);
     k = HKL(2);
     l = HKL(3);


  PART1 =  h^2/a^2 * (1-calp^2) + k^2/b^2 * (1-cbet^2) + l^2/c^2 * (1-cgam^2) +...
     2*h*k*(calp*cbet-cgam)/(a*b) + 2*h*l*(calp*cgam-cbet)/(a*c) + 2*k*l*(cbet*cgam-calp)/(b*c);

  PART2 = 1 - (calp^2 + cbet^2 + cgam^2) + 2*calp*cbet*cgam;

  sintlc = sqrt(PART1) / (2*sqrt(PART2));
