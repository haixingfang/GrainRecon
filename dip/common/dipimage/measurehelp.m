%MEASUREHELP   Provides help on the measurement features
%
% SYNOPSIS:
%  measurehelp
%
% This function provides help on the usage of the MEASUREMENT function.

% (C) Copyright 1999-2015               Pattern Recognition Group
%     All rights reserved               Faculty of Applied Physics
%                                       Delft University of Technology
%                                       Lorentzweg 1
%                                       2628 CJ Delft
%                                       The Netherlands
%
% Cris Luengo, February 2002.
% 12 Feb 2007:    Better way of removing Mike's private measurement functions.
% 5 March 2008:   Also listing Bernd's new derived measurements. (CL)
% 19 May 2015:    Removed duplicate code by calling MEASURE('-getfeatures'). (CL)

function out = measurehelp(in)

if nargin == 1
   if ischar(in) & strcmp(in,'DIP_GetParamList')
      out = struct('menu','Analysis',...
                   'display','Help on measurement features'...
                   );
      return
   end
end

msmts = measure('-getfeatures');

disp(' ')
disp('Measurement Features to use in MEASURE:')
disp(' ')
fprintf('%20s   %s\n','Name','Description')
fprintf('%20s---%s\n','-------------','--------------------------------------------------------')
for ii=1:length(msmts)
   fprintf('%20s - %s\n',msmts(ii).name,msmts(ii).description)
end
disp(' ')
disp('Measurements marked with an * require a grey-value input image.')
disp(' ')
