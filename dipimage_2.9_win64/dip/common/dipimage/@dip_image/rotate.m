%ROTATE   Rotate a vector image around an axis
%   ROTATE(V,AXIS,ANGLE) rotates the 3D vector image V about ANGLE
%   around the AXIS given by a second vector image or a vector.
%
%   ROTATE(V,ANGLE) rotates the 2D vector image V about ANGLE.
%
%   The ANGLE should be given between [0,2pi]
%
%   LITERATURE: Computer Graphics, D. Hearn and M.P. Baker, Prentice
%               Hall, p.408-419
%
%   SEE ALSO: rotation, rotation3d

% (C) Copyright 1999-2015               Pattern Recognition Group
%     All rights reserved               Faculty of Applied Physics
%                                       Delft University of Technology
%                                       Lorentzweg 1
%                                       2628 CJ Delft
%                                       The Netherlands
%
% Bernd Rieger, Nov 2000
% 26 July 2007:    Better testing for wrong input.
% 21 January 2015: Now also for 2D vectors (CL).

function out=rotate(v,axis,theta)
if ~builtin('isa',v,'dip_image')
   error('First input argument should be a vector image.');
end
if ~isvector(v)
   error('First input argument should be a vector image.');
end
s = builtin('numel',v);
switch(s)
   case 2
      if nargin~=2
         error('Not enough input arguments.');
      end
      theta = axis;
      if ~isnumeric(theta)
         error('Theta input argument must be numeric');
      end
      if theta<0 | theta>2*pi
         error('Angle out of range: 0,2pi');
      end
      s = sin(theta);
      c = cos(theta);
      M = [c,-s;s,c];
      out = M*v;
      
   case 3
      if nargin~=3
         error('Not enough input arguments.');
      end
      if ~isnumeric(theta) | numel(theta)~=1
         error('Theta input argument must be numeric scalar');
      end
      if theta<0 | theta>2*pi
         error('Angle out of range: 0,2pi');
      end
      if isnumeric(axis)
         s1 = size(axis);
      else
         s1 = imarsize(axis);
      end
      if prod(s1)~=3
         error('Axis needs to be 3D vector.');
      end
      M = dip_image('array',[3,3]);
      axis = axis./norm(axis).*sin(theta/2);
      a = axis(1);
      b = axis(2);
      c = axis(3);
      s = cos(theta/2);
      M(1,1) = 1-2*b*b-2*c*c;
      M(1,2) = 2*a*b-2*s*c;
      M(1,3) = 2*a*c+2*s*b;
      M(2,1) = 2*a*b+2*s*c;
      M(2,2) = 1-2*a*a-2*c*c;
      M(2,3) = 2*b*c-2*s*a;
      M(3,1) = 2*a*c-2*s*b;
      M(3,2) = 2*b*c+2*s*a;
      M(3,3) = 1-2*a*a-2*b*b;
      out = M*v;
      
   otherwise
      error('First input argument should be a 2D or 3D vector image.');
end



