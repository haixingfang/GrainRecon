%SUBSASGN   Overloaded operator for a{i}(j)=b.
%   Additionally defines:
%      a.pixelsize(j)  = b % sets the physical dimensions of the pixels.
%      a.pixelunits(j) = b % sets the units for the values above.
%      a.whitepoint    = b % sets the whitepoint for the XYZ color space.
%
%   Speed/memory usage:
%   Due to the way matlab handles custom implementations of subscripted
%   assignments, speed and memory usage are awful. A full copy is made of
%   the image being subscripted, it is then modified and finally the copy
%   replaces the original
%
%   dipsetpref('FastSubscriptedAssignment',1);
%
%   enables a much improved version, available since DIPimage v2.9.0. It is not
%   enabled by default, as a pure matlab version cannot avoid the copy, but
%   a mex file implementation can. However, since this implies a
%   reimplementation of a significant part of the original code paths,
%   there is a serious potential for bugs and regressions. The code has been
%   thoroughly tested, and the above dipsetpref does therefore come warmly
%   recommended. Still, it was decided to enforce a conscious decision to 
%   enable it as any bug could result in hard-to-track down errors.
%
%   i1=readim;
%   i1(100:199,100:199)=rr([100,100])
%
%   does not create a copy. However,
%
%   i1=readim;
%   i2=i1;
%   i1(100:199,100:199)=rr([100,100])
%
%   does create a copy: the initial "copy" "i2=i1;" does not actually copy
%   the image data, and as soon as we attempt to modify either i1 or i2,
%   the copy is made at that point, i.e. the copy is deferred until really
%   necessary. In fact, an easy way to enforce the copy is:
%
%   i2(0)=i2(0);

% (C) Copyright 1999-2014               Pattern Recognition Group
%     All rights reserved               Faculty of Applied Physics
%                                       Delft University of Technology
%                                       Lorentzweg 1
%                                       2628 CJ Delft
%                                       The Netherlands
%
% Cris Luengo, July 2000.
% 28 March 2001:    Fixed bug that allowed image to be resized upon assignment.
%                   It caused incorrectly defined dimensionality.
% 16 April 2001:    Not converting the data type of target any more.
% 14 August 2001:   Fixed bug with trailing singleton dimensions.
% 18 December 2001: Color information is not as volatile anymore.
% 6 February 2001:  Introduced same changes BR added to SUBSREF.
% 15 November 2002: Fixed binary images to work in MATLAB 6.5 (R13)
% 27 November 2003: Added changing value of single pixel in tensor images.
% February 2008:    Adding pixel dimensions and units to dip_image. (BR)
% 5 March 2008:     Allowing pixel dimensions to be 0. (CL)
% 6 March 2008:     Allowing slightly more complex syntax with property indexing.
%                   Fixed bug when A is an image array. (CL)
% 10 March 2008:    Added whitepoint property. (CL)
% 29 October 2014:  The profiler said IMARSIZE is nearly 30% of the execution time of SUBSREF. (CL)

function a = subsasgn(a,s,b)

if ~di_isdipimobj(a)
   % This only happens when a doesn't exist yet: a(:,:,1)=dip_image(b)
   error('Illegal syntax: please create image first.')
end
N = length(s);
arrayindex = 0;
imageindex = 0;
propindex = 0;

%arsz = imarsize(a); % This apparently has a large overhead? Why?
arsz = builtin('size',a);
for ii=1:N
   if strcmp(s(ii).type,'{}')
      if arrayindex | propindex
         error('Illegal indexing.')
      end
      arrayindex = ii;
      idx = s(arrayindex).subs;
      if length(idx) == 1
         if isnumeric(idx{1}) & (any(idx{1} > prod(arsz)) | any(idx{1} < 1))
            error('Index exceeds array dimensions.')
         end
      else
         if length(idx) ~= length(arsz)
            error('Number of array indices not the same as image array dimensionality.')
         end
         for jj=1:length(idx)
            if isnumeric(idx{jj}) & (any(idx{jj} > arsz(jj)) | any(idx{jj} < 1))
               error('Index exceeds array dimensions.')
            end
         end
      end
   elseif strcmp(s(ii).type,'()')
      if imageindex
         error('Illegal indexing on image array.')
      end
      imageindex = ii;
   else   % Dot-indexing
      % Allowed:
      %  a.p = b
      %  a.p(i) = b
      % Meaning: propindex [ imageindex ]
      if arrayindex | imageindex | propindex
         error('Illegal indexing.')
      end
      propindex = ii;
   end
end
if propindex
   % Set properties
   if imageindex
      idx = s(imageindex).subs;
      if length(idx)~=1
         error('Illegal indexing.')
      end
      idx = idx{1};
      if isequal(idx,':')
         idx = [];
      end
   else
      idx = [];
   end
   switch lower(s(propindex).subs)
   case 'pixelsize'
      if ~isnumeric(b) | any(b<0)
         error('Pixel size must be positive and numeric.');
      end
      for ii=1:prod(arsz)
         v = b;
         if isempty(idx)
            if length(v)==1
               v = repmat(v,1,a(ii).dims);
            elseif length(v)~=a(ii).dims;
               error('Pixel size array must match image dimesions.');
            end
            a(ii).physDims.PixelSize = v(:)';
         else
            if any(idx<=0) | any(idx>a(ii).dims)
               error('Index exceeds matrix dimensions.')
            end
            if length(v)~=1 & length(v)~=length(idx)
               error('Pixel size array must match image dimesions.');
            end
            a(ii).physDims.PixelSize(idx) = v(:)';
         end
      end
   case 'pixelunits'
      if isstr(b)
         b = {b};
      elseif ~iscellstr(b)
         error('Units array must be a string cell array.');
      end
      for ii=1:prod(arsz)
         v = b;
         if isempty(idx)
            if length(v)==1
               v = repmat(v,1,a(ii).dims);
            elseif length(v)~=a(ii).dims;
               error('Pixel size array must match image dimesions.');
            end
            a(ii).physDims.PixelUnits = v(:)';
         else
            if any(idx<=0) | any(idx>a(ii).dims)
               error('Index exceeds matrix dimensions.')
            end
            if length(v)~=1 & length(v)~=length(idx)
               error('Pixel size array must match image dimesions.');
            end
            a(ii).physDims.PixelUnits(idx) = v(:)';
         end
      end
   case 'whitepoint'
      if isempty(a(1).color)
         error('The image is not a color image');
      end
      if ~isnumeric(b) | prod(size(b))~=3
         error('The whitepoint must be a 1x3 numeric array.')
      end
      b = b(:)';
      for ii=1:prod(arsz)
         a(ii).color.xyz = b;
      end
   otherwise
      error('Illegal indexing.')
   end
else
   if ~di_isdipimobj(b)
      b = dip_image(b);
   end
   if iscolor(a)
      col = a(1).color;
   else
      col = '';
   end
   if ~isempty(col) & ~isempty(b(1).color)
      if ~strcmp(b(1).color.space,col.space)
         warning('Color spaces do not match: removing color space information.')
         col = '';
      end
   end
   % Assign images into array.
   if arrayindex & ~imageindex
      try
         a = builtin('subsasgn',a,substruct('()',s(arrayindex).subs),b);
      catch
         error(di_firsterr)
      end
      % Set color information.
      for ii=1:prod(arsz)
         a(ii).color = col;
      end
   % Assign pixels into images (possibly in array).
   elseif imageindex
      if arrayindex
         tmpa = builtin('subsref',a,substruct('()',s(arrayindex).subs));
         tmparsz = builtin('size',tmpa);
      else
         tmpa = a;
         tmparsz = arsz;
      end
      N = prod(tmparsz);
      if N > 1
         if ~istensor(tmpa)
            error('Cannot index by pixel in a dip_image_array that is not a tensor.')
         end
      end
      if ~istensor(b)
         error('Cannot assign into a tensor a dip_image_array that is not a tensor.')
      end
      M = prod(builtin('size',b));
      if M ~= 1 & M ~= N
         error('Image tensor sizes do not match.')
      end
      singlepixel = 0; % if M==1 and singlepixel==1, then size(b)==N
      sz = size(tmpa(1).data);
      sindx = s(imageindex);
      if length(sindx.subs) == 1
         if di_isdipimobj(sindx.subs{1})
            if ndims(squeeze(sindx.subs{1}))==0 % added for Piet
               sindx.subs{1}=double(sindx.subs{1})+1;
            elseif ~islogical(sindx.subs{1})
               error('Only binary images can be used to index.')
            else
               sindx.subs{1} = logical(sindx.subs{1}.data);
            end
         elseif ~ischar(sindx.subs{1}) & ~islogical(sindx.subs{1})
            sindx.subs{1} = sindx.subs{1}+1;
            if any(sindx.subs{1} > prod(sz))
               error('Index exceeds image dimensions.')
            end
            if prod(size(sindx.subs{1})) == 1
               singlepixel = 1;
            end
         end
         if islogical(sindx.subs{1})
            if ~isequal(size(sindx.subs{1}),sz)
               error('Mask image must match image size when indexing.')
            end
         end
      elseif length(sindx.subs) == tmpa(1).dims
         tmp = sindx.subs(2);
         sindx.subs(2) = sindx.subs(1);
         sindx.subs(1) = tmp;
         for ii=1:length(sz)
            if di_isdipimobj(sindx.subs{ii})
               if ndims(squeeze(sindx.subs{ii}))==0 %added for Piet
                  sindx.subs{ii}=double(sindx.subs{ii})+1;
               else
                  error('Illegal indexing on image.')
               end
            elseif islogical(sindx.subs{ii})
               error('Illegal indexing on image.')
            elseif ~ischar(sindx.subs{ii})
               sindx.subs{ii} = sindx.subs{ii}+1;
               if any(sindx.subs{ii} > sz(ii)) | any(sindx.subs{ii} < 1)
                  error('Index exceeds image dimensions.')
               end
            end
         end
         for ii=(length(sz)+1):length(sindx.subs) % indexing into the singleton dimensions
            if di_isdipimobj(sindx.subs{ii}) | islogical(sindx.subs{ii})
               error('Illegal indexing on image.')
            elseif ~ischar(sindx.subs{ii})
               if any(sindx.subs{ii} ~= 0)
                  error('Index exceeds image dimensions.')
               end
            end
         end
         if all(cellfun('length',sindx.subs)==1) & all(cellfun('isclass',sindx.subs,'double'))
            singlepixel = 1;
         end
         sindx.subs = sindx.subs(1:length(sz));
      else
         error('Number of indices not the same as image dimensionality.')
      end
      if M == 1 & N > 1
         % prepare
         tmpb = b;
         if singlepixel
            sz1 = tmparsz;
            sz2 = size(b.data);
            if isequal(sz1,sz2) | (prod(sz1)==max(sz1) & prod(sz2)==max(sz2) & prod(sz1)==prod(sz2))
               % assigning single tensor into tensor image, OR
               % vector into vector image (ignore vector orientation).
               b = dip_image('array',N);
               for ii=1:N
                  b(ii).data = tmpb.data(ii);
                  b(ii).dims = 0;
                  b(ii).dip_type = tmpb.dip_type;
               end
               M = N;
            else
               error('Image tensor sizes do not match.')
            end
         end
      end
      for ii=1:N
         if M == N
            out_type = tmpa(ii).dip_type;
            tmpb = dip_image(b(ii),out_type);
         else % M == 1
            out_type = tmpa(ii).dip_type;
            tmpb = dip_image(tmpb,out_type);
         end
         try
            tmpa(ii).data = subsasgn(tmpa(ii).data,sindx,tmpb.data);
         catch
            error(di_firsterr)
         end
         tmpa(ii).color = col;
      end
      if arrayindex
         a = builtin('subsasgn',a,substruct('()',s(arrayindex).subs),tmpa);
      else
         a = tmpa;
      end
   end
end
