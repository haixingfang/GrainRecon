%RADIALSUM   Computes the sum as a function of the R-coordinate
%
% SYNOPSIS:
%  image_out = radialsum(image_in,[mask],binSize,innerRadius)
%
% PARAMETERS:
%  image_in:    the input image
%  mask:        binary mask image (OPTIONAL)
%  binSize:     the size of the radial bins
%  innerRadius: the maximum radius to use: the smallest or largest radius
%               that fits the image
%
% DEFAULTS:
%  binSize = 1
%  innerRadius = 0
%
% NOTE:
%  The center of the image, around which the measurements are done, is
%  the same one as defined by the Fourier Transform. That is, on even
%  size, it is to the right of the true center. This is also the default
%  for functions like RR.

% (C) Copyright 1999-2009               Pattern Recognition Group
%     All rights reserved               Faculty of Applied Physics
%                                       Delft University of Technology
%                                       Lorentzweg 1
%                                       2628 CJ Delft
%                                       The Netherlands
%
% Cris Luengo, August 2001.
% 18 August 2009:    Added optional MASK parameter.
% 24 February 2016:  Added fast 2D implementation plus some comments (EB).

function image_out = radialsum(varargin)

d = struct('menu','Statistics',...
           'display','Radial sum',...
           'inparams',struct('name',       {'image_in',   'mask',      'binSize',        'innerRadius'},...
                             'description',{'Input image','Mask image','Radial bin size','Use inner radius'},...
                             'type',       {'image',      'image',     'array',          'boolean'},...
                             'dim_check',  {0,            0,           0,                []},...
                             'range_check',{[],           [],          'R+',             []},...
                             'required',   {1,            0,           0,                0},...
                             'default',    {'a',          '[]',        1,                0}...
                              ),...
           'outparams',struct('name',{'image_out'},...
                              'description',{'Output image'},...
                              'type',{'image'}...
                              )...
          );

% Check whether just the input-variables are being asked
if nargin == 1
   s = varargin{1};
   if ischar(s) & strcmp(s,'DIP_GetParamList')
      image_out = d;
      return
   end
end

% If possible, perform the faster radialsum (is approximate a factor 5 faster)
% We skip getparams() because that is a slow function, costs roughly the same
% as this method
if nargin == 1 && ndims(varargin{1}) == 2	    
  % Note that at this point image_in is not always an dipimage
  image_out = fast2D_radialsum(varargin{1});
  % To force output to be identical to the C-implementation (always a dipimage):
  image_out = mat2im(image_out);
  return
end

% The fast 2D radialsum is not applicable and it seems we have to run the
% C-implementation after parsing the input.

% Input parsing:
if nargin>=2
   s = varargin{2};
   if isnumeric(s) & prod(size(s))==1
      % This looks like the BINSIZE parameter, the user skipped the MASK parameter
      varargin = [varargin(1),{[]},varargin(2:end)];
   end
end
try
   [image_in,mask,binSize,innerRadius] = getparams(d,varargin{:});
catch
   if ~isempty(paramerror)
      error(paramerror)
   else
      error(firsterr)
   end
end

% Run the C-implementation
image_out = dip_radialsum(image_in,mask,[],binSize,innerRadius,[]);

end

% This method computes the radialsum for a 2D matrix without a mask or bin-size.
% It works by creating a labels matrix and then one call to accumarray().
function [ y ] = fast2D_radialsum( A )

  % Often radialsum() is called in sucsession for different matrixes of the
  % same size. To reduce computation time we use persistent variables to store
  % the labels matrix between different calls.
  persistent lookup_size lookup_M

  % Make sure the input is a matrix and not accidently a dipimage.
  if ~isnumeric(A), A = im2mat(A); end

  % Check whether the previous computed labels matrix corresponds to the size
  % of the input matrix.
  if ~isempty(lookup_size)
    if any(size(A) ~= lookup_size)
      lookup_size = [];
    end
  end

  % If we have to recompute the labels matrix
  if isempty(lookup_size)

    lookup_size = size(A);
    lookup_M = zeros(size(A));
    M = ceil(size(A) / 2);

    % Check if we can use more symmetry properties
    if size(A, 1) == size(A, 2)
      s = size(A, 1);
      for i = 1:ceil(s/2)
        for j = i:ceil(s/2)
          m = floor(norm([i;j] - M')) + 1;
          lookup_M(i, j) = m;
          lookup_M(i, s-j+1) = m;
          lookup_M(j, i) = m;
          lookup_M(j, s-i+1) = m;
          lookup_M(s-i+1, j) = m;
          lookup_M(s-i+1, s-j+1) = m;
          lookup_M(s-j+1, i) = m;
          lookup_M(s-j+1, s-i+1) = m;
        end
      end
    else
      for i = 1:ceil(size(A, 1)/2)
        for j = 1:ceil(size(A, 2)/2)
          m = floor(norm([i;j] - M')) + 1;
          lookup_M(i, j) = m;
          lookup_M(size(A, 1)-i+1, j) = m;
          lookup_M(i, size(A, 2)-j+1) = m;
          lookup_M(size(A, 1)-i+1, size(A, 2)-j+1) = m;
        end
      end
    end

    % Store it in a way that we can use it with accumarray
    lookup_M = uint32(lookup_M(:));

  end

  % Prepare the input and remove NaN elements
  A = A(:);
  A(isnan(A)) = 0;

  % Let Matlab perform its magic
  y = accumarray(lookup_M, A);

end
