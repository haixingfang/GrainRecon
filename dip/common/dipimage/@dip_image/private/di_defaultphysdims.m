%OUT_PHYS = DI_DEFAULTPHYSDIMS(N)

% (C) Copyright 1999-2014               Pattern Recognition Group
%     All rights reserved               Faculty of Applied Physics
%                                       Delft University of Technology
%                                       Lorentzweg 1
%                                       2628 CJ Delft
%                                       The Netherlands
%
% Cris Luengo 2008.
% 29 October 2014:  Significant speedup.

function out_type = di_defaultphysdims(n)
%out_type.PixelSize = repmat(1,1,n);
out_type.PixelSize = ones(1,n);
%out_type.PixelUnits = repmat({'px'},1,n);
pu = cell(1,n);
pu(:) = {'px'};
out_type.PixelUnits = pu;
