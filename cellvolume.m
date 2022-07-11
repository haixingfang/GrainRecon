function V = cellvolume(ucell)

a   = ucell(1);
b   = ucell(2);
c   = ucell(3);
calp = cos(ucell(4)*pi/180.);
cbet = cos(ucell(5)*pi/180.);
cgam = cos(ucell(6)*pi/180.);

angular = sqrt(1 - calp^2 - cbet^2 - cgam^2 + 2*calp*cbet*cgam);
V = a*b*c*angular;                                           %Volume of unit cell
