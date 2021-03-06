%GETLINKPIXEL   Get link-pixels from skeleton
%
% SYNOPSIS:
%  image_out = getlinkpixel(image_in)

% (C) Copyright 1999-2000               Pattern Recognition Group
%     All rights reserved               Faculty of Applied Physics
%                                       Delft University of Technology
%                                       Lorentzweg 1
%                                       2628 CJ Delft
%                                       The Netherlands
%
% Lucas van Vliet, May 2000.
% Cris Luengo, August 2001: Using COUNTNEIGHBOURS.

function image_out = getlinkpixel(varargin)

d = struct('menu','Binary Filters',...
           'display','Get link-pixels from skeleton',...
           'inparams',struct('name',       {'image_in'},...
                             'description',{'Input image'},...
                             'type',       {'image'},...
                             'dim_check',  {0},...
                             'range_check',{'bin'},...
                             'required',   {1},...
                             'default',    {'a'}...
                            ),...
           'outparams',struct('name',{'image_out'},...
                              'description',{'Output image'},...
                              'type',{'image'}...
                              )...
          );
if nargin == 1
   s = varargin{1};
   if ischar(s) & strcmp(s,'DIP_GetParamList')
      image_out = d;
      return
   end
end
try
   [image_in] = getparams(d,varargin{:});
catch
   if ~isempty(paramerror)
      error(paramerror)
   else
      error(firsterr)
   end
end

image_out = (countneighbours(image_in,1,ndims(image_in))==2);
