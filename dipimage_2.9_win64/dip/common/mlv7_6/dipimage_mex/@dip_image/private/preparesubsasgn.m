%SUBSASGN   Overloaded operator for a{i}(j)=b.
%   Additionally defines:
%      a.pixelsize(j)  = b % sets the physical dimensions of the pixels.
%      a.pixelunits(j) = b % sets the units for the values above.
%      a.whitepoint    = b % sets the whitepoint for the XYZ color space.

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
%
% July 2017:        New hybrid mex/matlab implementation to allow in-place
%                   subscripted assignment.
% 12-8-2017 - MvG - Fix for binary images. Binary and complex are the most
%                   likely candidates for trouble.
% 19-8-2017 - MvG - Missed the use case im(:)= completely. Two subcases
%                   exist: im(:)=pl(:) and im(:)=<scalar>.

% This file containing the matlab code is a refactored version of the
% original  subsasgn.m. When a subscripted assignment is invoked, the
% subsasgn.mex file initially takes control. This is the only way to
% prevent spurious references to the input data which trigger unwanted
% copies of the data. Then control is passed to this m-file which performs
% the vast  majority of checking the parameters, type of indexing etc.
% Control is then  passed back to the mex-file which performs the actual
% image data copy if required.
%
% Communication between the different functions is performed by passing
% the "idi" (short for index information) structure around.
%
% Critical warning: the mex code does not do much checking. It expects
% given types for given fields. Be very careful to cross check against
% the mex code when you modify this code.
%
% Nomenclature:
%
%   {} value subscripting:
%      indexes into the value (colour/tensor) dimensions
%   () spatial subscripting
%   .  property subscripting, provides access to
%      pixelsize/color properties
%
%   payload : the data that is copied into the dip_image, i.e. the
%             right-hand side of the assignment
%
% Code walk-through:
%
%   validate_subscript_types_and_check_value_subscripts
%     checks whether the index types are compatible and in the correct
%     order. Permissible are:
%     single subscript actions  {}, () or .
%     two subscripting actions: {}() or .()
%
%     This is enforced by the three checks, summarised by the following
%     matrix:
%
%                 second         +: allowed
%                () {}  .        M: forbidden by the matlab interpreter
%             ()  2  M  3        1:     "      " rule 1 
%     first   {}  +  1  3        2:     "      " rule 2
%             .   +  1  3        3:     "      " rule 3
%
%     idi.si_value, idi.si_spatial and idi.si_property store which
%     element of the subs array stores the subscripts for that particular
%     type of indexing. Value indexes can be checked for validity
%     without further information and this check is performed here.
%     (si_ for subscript index)
%
%   the code then proceeds to deal with each of three types of
%   subscripting. The first two can be dealt with entirely here in the
%   matlab code. The result is passed on to the subsasgn mex file which
%   passes it on without further processing. It is instructed to do so
%   using the idi.pass_on field. The three types are:
%
%   1. property subscripting, performed by <process_property_subscripting>
%   2. pure value subscripting (no spatial subscripting)
%   3. subscripting involving spatial subscripting
%
%   for types 1&2 control is now passed back to the mex file. The
%   reminder of the code performs checks, and prepares for the copying
%   of the image data that will be done by the mex file.
%
%   normalise_value_subscripts_and_add_spatial_info
%     The value subscript(s) can be none, one-dimensional or
%     multi-dimensional. By adapting a linear value index all these
%     cases can be dealt with in the same way. The fields added to idi are:
%
%     .li_value     the linear index into the value dimensions
%     .li_value_N   the length of the linear index array (for easy access
%                   by the mex code)
%     .sbvl_dims    the dimensions of the dip_image along the value
%                   dimensions after subscripting (for checking against
%                   the payload dimensions)
%     .sp_dims      the spatial dimensions
%     .sp_dm        spatial dimensionality
%     .sp_N         number of elements in the spatial dimensions
%
%   check_value_assignment_against_payload
%     Checks wether the payload value dimensions are compatible with
%     the value dimension subscripts. I did not change this code, it
%     tests only for the number of elements, not the layout... not sure
%     if this was the intention.
%
%     .li_payload_value_N       number of value dimension elements of the payload
%
%   verify_and_adapt_spatial_subscripts
%     Three styles of subscripting are supported for the spatial
%     dimensions:
%
%     - an array (or ':') for each dimension    (sp_style=0)
%     - a single linear addressing array        (sp_style=1)
%     - a mask image                            (sp_style=2)
%
%     the code identifies which style it is and stores it in the
%     idi.sp_style field. Then it verifies that the subscripts are
%     compatible with the image dimensions and where necessary adapts them
%     to support some special cases.
%
%     Partly through the helper function <normalise_spatial_subscripts>,
%     it further does the following:
%
%     - it determines idi.assign_to_single_pixel
%     - if sp_style=0:
%       + determines idi.first_sp_sub_contiguous
%       + any ':' subscript is converted to an array 0:(size(dim)-1)
%       + any logical array is converted to an array of indices
%         (through find)
%       + any non-double numerical array is converted to double
%       (all of these are to simplify the mex code)
%
%     assign_to_single_pixel is later used to support the special case
%     of assigning a spatial payload into a single multi-valued pixel.
%
%     first_sp_sub_contiguous reflects whether the indices along the first
%     axis (i.e. actually the y-axis due to the x,y dimension swapping
%     thing dip_image's do) are contiguous. This allows the mex file to
%     use memcpy.
%
%   deal_with_special_tensor_and_vector_assignments
%     creates a mock-up version of the payload such that the following
%     special case requires no further separate attention in the remainder
%     of the code:
%
%     Example:
%       vldims=[3,2];
%       spdims=[100,100];
%       tim=dip_image(repmat({newim(spdims)},vldims));
%       payload=newim(vldims([2,1]));  % spatial dimensions get swapped,
%                                      % but not the value dimensions!
%       tim(50,40)=payload;
%
%   create_linear_index_into_payload_and_prep_data_type_conversion
%     to match the linear indexing into the assignee value dimensions,
%     a linear index into the payload is set up. Before copying a
%     payload .data field, it must be converted to match the data type
%     of the corresponding assignee field. To facilitate passing this
%     information around, it is stored in .payload_convert_to. The
%     conversion itself is performed by a call to payloadsubsasgn() from
%     the mex code. Finally, to copy the data we do not care about the
%     actual data type, only its size in bytes.
%
%     idi.li_payload_value     linear index into payload value dimensions
%     idi.payload_convert_to   list of data types
%     idi.sizeof               sizeof each of the above data types
%     idi.payload_is_scalar    if true, the payload can be "broadcast"
%
%   check_spatial_dimensions
%     check the spatial dimensions of the assignee and payload are
%     compatible.


function [oa,b,idi] = preparesubsasgn(a,s,b)
   mcode = 0;
   if ~dipgetpref('FastSubscriptedAssignment')
      mcode = 1;
   else
      cv = matlabver_ge([]);
      % MATLAB 9.2 	R2017a 	r37   last tested version
      if vngt(cv,[9,2])
         disp('Dear developer, please check whether fast subscripted');
         disp('assignment works with your shiny new version of matlab');
         disp('and update dipimage/xfiles/@dip_image/private/preparesubsasgn.m');
         disp('accordingly');
         mcode = 1;
      end
   end
   if mcode
      oa = mcode_subsasgn(a,s,b);
      idi.pass_on = 1;
      return;
   end
   
   % in-place notation a=preparesubsasgn(a,s,b) spells disaster. matlab
   % gets confused and starts needlessly copying data. Hence the explicitly
   % different name "oa" for the output parameter corresponding to the image
   oa=[];

   if ~di_isdipimobj_local(a)
      % This only happens when a doesn't exist yet: a(:,:,1)=dip_image(b)
      error('Illegal syntax: please create image first.')
   end
   %arsz = imarsize(a); % This apparently has a large overhead? Why?
   idi = [];
   idi.p_subs =s;
   idi.vl_size = builtin('size',a);
   idi.pass_on = 0;
   idi.large_dims = strcmp(class(di_cast_to_mwindex(1)),'int64');
   
   idi=validate_subscript_types_and_check_value_subscripts(idi);
   % First category of assignment: property assignment
   if idi.si_property
      oa=process_property_subscripting(a,b,idi);
      idi.pass_on = 1;
      return;
   end

   if ~di_isdipimobj_local(b)
      b = dip_image(b);
   end
   idi.col = determine_color_space(a,b);
   % Second category of assignment: images into (cell/tensor) array
   if idi.si_value & ~idi.si_spatial
      oa = process_pure_value_subscripting(a,b,idi);
      idi.pass_on = 1;
      return;
   end
   % Third and final category: assign sub-(tensor)images into (tensor) images
   idi = normalise_value_subscripts_and_add_spatial_info(a,idi);
   idi = check_value_assignment_against_payload(a,b,idi);  % throws errors if necessary
   idi = verify_and_adapt_spatial_subscripts(idi);
   % the following sets up a mock-up b if necessary to accomodate two special cases
   [b,idi] = deal_with_special_tensor_and_vector_assignments(b,idi);
   idi = create_linear_index_into_payload_and_prep_data_type_conversion(a,b,idi);
   check_spatial_dimensions(b,idi);
   % convert to 0-based indexing where appropriate - do not use in matlab scope
   % from now on
   idi.li_value = idi.li_value-1;
   idi.li_payload_value = idi.li_payload_value-1;
end

% +---------------------------------------------------------------------------+
% | Above we handed back control to the mex-file. Below are helper functions  |
% | for mcode-handled cases and further on helper functions to prepare for    |
% | the image copy (indicated by a marker block like this one)                |
% +---------------------------------------------------------------------------+

function idi = validate_subscript_types_and_check_value_subscripts(idi)
   s = idi.p_subs;
   arsz = idi.vl_size;
   N = length(s);
   idi.si_value = 0;
   idi.si_spatial = 0;
   idi.si_property = 0;
   for ii=1:N
      if strcmp(s(ii).type,'{}')
         if idi.si_value | idi.si_property
            % Seems to be forbidden at the interpreter level anyway - MvG
            error('Illegal indexing.')
         end
         idi.si_value = ii;
         idx = s(idi.si_value).subs;
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
         if idi.si_spatial
            % Seems to be forbidden at the interpreter level anyway - MvG
            error('Illegal indexing on image array.')
         end
         idi.si_spatial = ii;
      else   % Dot-indexing
         % Allowed:
         %  a.p = b
         %  a.p(i) = b
         % Meaning: idi.si_property [ idi.si_spatial ]
         if idi.si_value | idi.si_spatial | idi.si_property
            error('Illegal indexing.')
         end
         idi.si_property = ii;
      end
   end
end


function a = process_property_subscripting(a,b,idi)
   % Set properties
   s = idi.p_subs;
   arsz = idi.vl_size;
   if idi.si_spatial
      idx = s(idi.si_spatial).subs;
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
   switch lower(s(idi.si_property).subs)
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
end

function col = determine_color_space(a,b)
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
end

function oa = process_pure_value_subscripting(a,b,idi);
   try
      oa = builtin('subsasgn',a,substruct('()',idi.p_subs(idi.si_value).subs),b);
   catch
      error(di_firsterr)
   end
   % Set color information.
   for ii=1:prod(idi.vl_size)
      oa(ii).color = idi.col;
   end
end

% +---------------------------------------------------------------------------+
% | Below here the helper functions to prepare for the actual image data copy |
% +---------------------------------------------------------------------------+

function idi=verify_and_adapt_spatial_subscripts(idi)
   idi.assign_to_single_pixel = 0; % if M==1 and singlepixel==1, then size(b)==N
   idi.full_linear_copy=logical(0);
   idi.first_sp_sub_contiguous = logical(0);
   sz = idi.sp_dims;
   sindx = idi.p_subs(idi.si_spatial);
   idi.sp_style=0;  % make life easier for the mex-file, 0="simple indexing"
   if length(sindx.subs) == 1
      idi.sp_style=-1;
      if di_isdipimobj_local(sindx.subs{1})
         if ndims(squeeze(sindx.subs{1}))==0 % added for Piet
            sindx.subs{1}=double(sindx.subs{1});
         elseif ~islogical(sindx.subs{1})
            error('Only binary images can be used to index.')
         else
            sindx.subs{1} = logical(sindx.subs{1}.data);
         end
         idi.sp_style=1;   % indexing with logicals (i.e. binary)
      elseif ~ischar(sindx.subs{1}) & ~islogical(sindx.subs{1})
         if any(sindx.subs{1} >= prod(sz))
            error('Index exceeds image dimensions.')
         end
         if prod(size(sindx.subs{1})) == 1
            idi.assign_to_single_pixel = 1;
         end
         idi.sp_style=2;  % linear array indexing
         if ~strcmp(class(sindx.subs{1}),'double')
            % for mex/copy_using_linear_index()
            sindx.subs{1}=double(sindx.subs{1});
         end
      end
      if islogical(sindx.subs{1})
         if ~isequal(size(sindx.subs{1}),sz)
            error('Mask image must match image size when indexing.')
         end
      end
      idi.im_subs=sindx.subs;
      if ischar(idi.im_subs{1})
        if ~strcmp(idi.im_subs{1},':')
          error('Unhandled case for string index ~= :');
        end
        idi.full_linear_copy=logical(1);
        idi.sp_style=2;
      end
      if idi.sp_style==-1
        error('Unhandled indexing style');
      end
   elseif length(sindx.subs) == idi.sp_dm
      tmp = sindx.subs(2);
      sindx.subs(2) = sindx.subs(1);
      sindx.subs(1) = tmp;
      for ii=1:length(sz)
         if di_isdipimobj_local(sindx.subs{ii})
            if ndims(squeeze(sindx.subs{ii}))==0 %added for Piet
               sindx.subs{ii}=double(sindx.subs{ii});
            else
               error('Illegal indexing on image.')
            end
         elseif islogical(sindx.subs{ii})
            error('Illegal indexing on image.')
         elseif ~ischar(sindx.subs{ii})
            sindx.subs{ii} = sindx.subs{ii};
            if any(sindx.subs{ii} >= sz(ii)) | any(sindx.subs{ii} < 0)
               error('Index exceeds image dimensions.')
            end
         end
      end
      for ii=(length(sz)+1):length(sindx.subs) % indexing into the singleton dimensions
         if di_isdipimobj_local(sindx.subs{ii}) | islogical(sindx.subs{ii})
            error('Illegal indexing on image.')
         elseif ~ischar(sindx.subs{ii})
            if any(sindx.subs{ii} ~= 0)
               error('Index exceeds image dimensions.')
            end
         end
      end
      if all(cellfun('length',sindx.subs)==1) & all(cellfun('isclass',sindx.subs,'double'))
         idi.assign_to_single_pixel = 1;
      end
      sindx.subs = sindx.subs(1:length(sz));
      idi=normalise_spatial_subscripts(idi,sindx);
   else
      error('Number of indices not the same as image dimensionality.')
   end
   idi.im_subs_N=length(idi.im_subs);
   idi.im_subs_dims=di_cast_to_mwsize(cellfun('length',idi.im_subs));
end

function idi=normalise_spatial_subscripts(idi,sindx)
   for ii=1:length(sindx.subs)
      switch class(sindx.subs{ii})
         case 'char'
            if ~strcmp(sindx.subs{ii},':')
               error('expected ":", got %s',sindx.subs{ii});
            end
            sindx.subs{ii}=di_cast_to_mwindex(0:(idi.sp_dims(ii)-1));
         case 'logical'
            sindx.subs{ii}=find(sindx.subs{ii})-1;
         otherwise
            % for mex/copy_contiguous_line() and mex/copy_regular_subscripts()
            sindx.subs{ii}=di_cast_to_mwindex(sindx.subs{ii});
      end
   end
   % if the first dimension is contiguous, use memcpy to speed things up
   if strcmp(class(sindx.subs{1}),'double') && all(diff(sindx.subs{1})==1)
      idi.first_sp_sub_contiguous=logical(1);
   end
   idi.im_subs=sindx.subs;
end


function idi=normalise_value_subscripts_and_add_spatial_info(a,idi)
   % set up a regular array in the same shape as the (potentially tensor) image a
   model = zeros(idi.vl_size);
   model(:) = 1:prod(idi.vl_size);
   if idi.si_value
     model=model(idi.p_subs(idi.si_value).subs{:});
   end
   % now we can trivially obtain a list of linear indices into the original array
   idi.li_value = model(:);
   idi.li_value_N = length(idi.li_value);
   idi.sbvl_dims = size(model);
   idi.sp_dims = size(a(idi.li_value(1)).data);
   idi.sp_dm = a(idi.li_value(1)).dims;
   idi.sp_N = prod(idi.sp_dims);
   idi.sp_strides = di_cast_to_mwindex( cumprod([1,idi.sp_dims(1:end-1)]));
end

function idi=check_value_assignment_against_payload(a,b,idi)
   if idi.li_value_N > 1
      if ~is_tensor_subarray(a,idi)
         error('Cannot index by pixel in a dip_image_array that is not a tensor.')
      end
   end
   if ~istensor(b)
      error('Cannot assign into a tensor a dip_image_array that is not a tensor.')
   end
   idi.li_payload_value_N = prod(builtin('size',b));
   if idi.li_payload_value_N ~= 1 & idi.li_payload_value_N ~= idi.li_value_N
      error('Image tensor sizes do not match.')
   end
end

% this reflects the original subsasgn code: the original is allowed to be
% non-tensor as long as the subarray (indicated by s(arrayindex) ) _is_
% a tensor
function istensor=is_tensor_subarray(a,idi)
   istensor=logical(1);
   refsize=size(a(idi.li_value(1)).data);
   for ii=idi.li_value(2:end)';
      if ~isequal(refsize,size(a(ii).data))
         istensor = logical(0);
         return;
      end
   end
end

function [b,idi]=deal_with_special_tensor_and_vector_assignments(b,idi)
   if idi.li_payload_value_N == 1 & idi.li_value_N > 1 & idi.assign_to_single_pixel
      sz1 = idi.sbvl_dims;
      sz2 = size(b.data);
      if isequal(sz1,sz2) | (prod(sz1)==max(sz1) & prod(sz2)==max(sz2) & prod(sz1)==prod(sz2))
         % assigning single tensor into tensor image, OR
         % vector into vector image (ignore vector orientation).
         bdata = b.data;
         btype = b.dip_type;
         b = dip_image('array',idi.li_value_N);
         for ii=1:idi.li_value_N
            b(ii).data = bdata(ii);
            b(ii).dims = 0;
            b(ii).dip_type = btype;
         end
         idi.li_payload_value_N = idi.li_value_N;
      else
         error('Image tensor sizes do not match.')
      end
   end
end

function idi=create_linear_index_into_payload_and_prep_data_type_conversion(a,b,idi)
   if idi.li_payload_value_N==1
     idi.li_payload_value=ones(1,idi.li_value_N);
   else
     idi.li_payload_value=1:idi.li_value_N;
   end
   idi.payload_convert_to={a(idi.li_value).dip_type};
   types=idi.payload_convert_to;
   % 12-8-2017 - MvG - fix for binary images
   binmask=strcmp(types,'bin');
   if any(binmask)
      for dtix=find(binmask)
         types{dtix}=class(a(idi.li_value(dtix)).data);
      end
   end
   idi.payload_is_scalar=prod(size(b(1).data))==1;
   typeinfo=di_dtinfo(types);
   % call again to obtain sizeof=4,8 for scomplex,dcomplex
   typeinfo=di_dtinfo({typeinfo.diplib_real});
   idi.sizeof=int32([typeinfo.sizeof]);
end

function check_spatial_dimensions(b,idi)
   % in the original subsasgn.m these checks were implicitly done
   % calling subsasgn() on the .data field.
   if idi.sp_style==0 && idi.payload_is_scalar
      return;
   end
   if idi.sp_style==0
      plsz=size(b(1).data);
      sbsz=cellfun('length',idi.im_subs);
      % remove singleton dimensions from the comparison
      plsz(plsz==1)=[];
      sbsz(sbsz==1)=[];
      if ~isequal(plsz,sbsz)
         error('Subscripted assignment dimension mismatch');
      end
   end
end

function b = di_isdipimobj_local(a)
  %#function isa
  b = builtin('isa',a,'dip_image');
end

function r = vngt(a,b) % >=
   r = 0;
   if ( a(1)>b(1) ) | ( a(1)==b(1) & a(2)>b(2) )
      r = 1;
   end
end
