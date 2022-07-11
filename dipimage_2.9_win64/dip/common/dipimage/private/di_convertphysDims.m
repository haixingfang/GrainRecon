%DI_CONVERTPHYSDIMS    Changes the power on units.
%   DI_CONVERTPHYSDIMS(IN,OP,N) applies the operation OP with
%   parameter N to the string IN.
%
%   The string IN represents units of measurement, but are expected
%   to be simple, such as 'mm^2', 'micron^-1', 'km', etc. The expected
%   format is '<base>^<exp>', with <exp> any number (doubles) and
%   <base> any string not containing a '^' character. However, if the
%   base represents the multiplication of two or more units, then the
%   result will not be correct.
%
%   OP is any of the strings 'invert', '+' or '-'. N is a number.
%
%   Input or output <base> can be '1' or 'undefined', handled correctly.
%
%  EXAMPLES:
%   OP '+', N=1 : m -> m^2  , m^2 -> m^3  , m^-1 -> 1
%   OP '+', N=2 : m -> m^3  , m^2 -> m^4  , m^-1 -> m
%   OP '-', N=1 : m -> 1    , m^2 -> m^1  , m^-1 -> m^-2
%   OP '-', N=2 : m -> m^-1 , m^2 -> 1    , m^-1 -> m^-3
%   OP 'invert' : m -> m^-1 , m^2 -> m^-2 , m^-1 -> m
%   OP 'base'   : m -> m    , m^2 -> m    , m^-1 -> m

% (C) Copyright 2015 by Cris Luengo.
% 19 May 2015.
% Centre for Image Analysis, Uppsala, Sweden.

function out = di_convertphysDims(in,op,n)

if nargin<3
   n = 0;
elseif nargin<2
   error('More input arguments needed!')
end

if iscell(in)
   out = in;
   for ii=1:length(in)
      out{ii} = dcpd(in{ii},op,n);
   end
else
   out = dcpd(in,op,n);
end



function out = dcpd(in,op,n)

ii = find(in=='^');
if length(ii)==0
   base = in;
   exp = 1;
elseif length(ii)==1
   base = in(1:ii-1);
   exp = str2double(in(ii+1:end));
else
   exp = NaN;
end

if isnan(exp)
   warning('Physical dimensions units have unexpected format.')
   out = in;
   return
end
if strcmp(base,'1') || strcmp(base,'undefined')
   out = base;
   return
end

switch lower(op)
   case 'invert'
      exp = -exp;
   case '+'
      exp = exp+n;
   case '-'
      exp = exp-n;
   otherwise
      error('Unkown option.');
end

if exp==1
   out = base;
elseif exp==0
   out = '1';
else
   out = sprintf('%s^%d',base,exp);
end
