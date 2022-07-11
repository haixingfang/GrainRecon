%DATATYPE  Extracts the data type string from a dip_image.
%   A = DTINFO(B) returns info about the data type of B (image or as string)
%   DIPlib convention of data type names.
%

% (C) Copyright 1999-2007               Pattern Recognition Group
%     All rights reserved               Faculty of Applied Physics
%                                       Delft University of Technology
%                                       Lorentzweg 1
%                                       2628 CJ Delft
%                                       The Netherlands
%
% Michael van Ginkel, July 2017

function di = dtinfo(varargin)
  persistent dtinfo_table;  % to speed up things
  persistent dtinfo_keys;

  di=[];
  if nargin == 1
    dtype = varargin{1};
    if ischar(dtype) && strcmp(dtype,'DIP_GetParamList')
      di = [];
      return
    end
  end

  if isempty(dtinfo_table)
    dtinfo_table=       struct('diplib','uint8', 'sizeof',1,'classid','mxUINT8_CLASS', 'matlab','uint8', 'mex','uint8_T', 'c','unsigned char');
    dtinfo_table(end+1)=struct('diplib','uint16','sizeof',2,'classid','mxUINT16_CLASS','matlab','uint16','mex','uint16_T','c','unsigned short');
    dtinfo_table(end+1)=struct('diplib','uint32','sizeof',4,'classid','mxUINT32_CLASS','matlab','uint32','mex','uint32_T','c','unsigned int');
    dtinfo_table(end+1)=struct('diplib','sint8', 'sizeof',1,'classid','mxINT8_CLASS',  'matlab','int8',  'mex','uint8_T', 'c','signed char');
    dtinfo_table(end+1)=struct('diplib','sint16','sizeof',2,'classid','mxINT16_CLASS', 'matlab','int16', 'mex','uint16_T','c','signed short');
    dtinfo_table(end+1)=struct('diplib','sint32','sizeof',4,'classid','mxINT32_CLASS', 'matlab','int32', 'mex','uint32_T','c','signed int');
    dtinfo_table(end+1)=struct('diplib','sfloat','sizeof',4,'classid','mxSINGLE_CLASS','matlab','single','mex','float',   'c','float');
    dtinfo_table(end+1)=struct('diplib','dfloat','sizeof',8,'classid','mxDOUBLE_CLASS','matlab','double','mex','double',  'c','double');
    [dtinfo_table.iscomplex]=deal(logical(0));
    [dtinfo_table.diplib_real]=deal(dtinfo_table.diplib);
    dtinfo_table(end+1)=struct('diplib','scomplex','sizeof',[],'classid',[],'matlab',[],'mex',[],'c',[],'iscomplex',logical(1),'diplib_real','sfloat');
    dtinfo_table(end+1)=struct('diplib','dcomplex','sizeof',[],'classid',[],'matlab',[],'mex',[],'c',[],'iscomplex',logical(1),'diplib_real','dfloat');
    [dtinfo_table.real_equivalent]=deal([]);
    keys=[{dtinfo_table.diplib};num2cell(1:length(dtinfo_table))];
    keys=keys(:);
    dtinfo_keys=struct(keys{:});
    dtinfo_table(dtinfo_keys.scomplex).real_equivalent=...
        dtinfo_table(dtinfo_keys.sfloat);
    dtinfo_table(dtinfo_keys.dcomplex).real_equivalent=...
        dtinfo_table(dtinfo_keys.dfloat);
  end

  if nargin==0
    di=dtinfo_table;
    return;
  end
  if nargin~=1
    error('single argument expected');
  end
  if isempty(dtype)
    di=dtinfo_table;
    return;
  end
  if isa(dtype,'dip_image')
    dtype=datatype(dtype);
  end
  if isa(dtype,'char')
    dtype={dtype};
  end
  if ~(iscell(dtype) && all(cellfun('isclass',dtype,'char')) )
    error('argument must be a dip_image or a (cell array of) string(s) containing a data type');
  end
  % fetch the correct elements from dtinfo_table
  di=cellfun(@(x) dtinfo_table(dtinfo_keys.(x)),dtype,'uniformoutput',true);
end







