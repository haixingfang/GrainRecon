%MEASURE   Do measurements on objects in an image
%
% SYNOPSIS:
%  msr = measure(object_in,gray_in,measurmentIDs,objectIDs,...
%                 connectivity,minSize,maxSize)
%
% PARAMETERS:
%  object_in: binary or labelled image holding the objects.
%  gray_in: (original) gray value image of object_in. It is needed for
%           several types of measurements. Otherwise you can use [].
%  measurementIDs: measurements to be performed, either a single string
%                  or a cell array with strings (e.g.: {'Size','Perimeter'} ).
%                  See MEASUREHELP for a full list of possible
%                  measurements.
%  objectIDs: labels of objects to be measured. Use [] to measure all
%             objects
%  connectivity: defines which pixels are considered neighbours: up to
%     'connectivity' coordinates can differ by maximally 1. Thus:
%     * A connectivity of 1 indicates 4-connected neighbours in a 2D image
%       and 6-connected in 3D.
%     * A connectivity of 2 indicates 8-connected neighbourhood in 2D, and
%       18 in 3D.
%     * A connectivity of 3 indicates a 26-connected neighbourhood in 3D.
%     Connectivity can never be larger than the image dimensionality.
%     Setting the connectivity to Inf (the default) makes it equal to the image
%     image dimensionality.
%  minSize, maxSize: minimum and maximum size of objects to be measured.
%
% DEFAULTS:
%  measurementIDs = 'size'
%  objectIDs = []
%  connectivity = inf
%  minSize = 0
%  maxSize = 0
%
% RETURNS:
%  msr: a dip_measurement object containing the results.
%
% EXAMPLE:
%  img = readim('cermet')
%  msr = measure(img<100, img, ({'size', 'perimeter','mean'}), [], ...
%                1, 1000, 0)
%
% NOTE:
%  The function MEASUREHELP provides help on the measurement
%  features available in this function.
%
% NOTE:
%  Several measures use the boundary chain code (i.e. 'Feret', 'Perimeter',
%  'CCBendingEnergy', 'P2A' and 'PodczeckShapes'). These measures will fail
%  if the object is not compact (one chain code must represent the whole
%  object). If the object is not compact under the connectivity chosen, only
%  one part of the object will be measured. Make sure the connectivity in
%  MEASURE matches that used in LABEL!

% 'UN-DOCUMENTED' SYNTAX:
%   measure -getfeatures
%   m = measure('-getfeatures')
%      Returns a struct with name and description of all recognized features.

% (C) Copyright 1999-2015               Pattern Recognition Group
%     All rights reserved               Faculty of Applied Physics
%                                       Delft University of Technology
%                                       Lorentzweg 1
%                                       2628 CJ Delft
%                                       The Netherlands
%
% Cris Luengo, August 2000.
% September 2000: Added 'optionarray' for measurementIDs.
% September 2000: Added MINSIZE and MAXSIZE parameters.
% December 2000:  Fixed bug when using MINSIZE or MAXSIZE with new dip_measurement.
% January 2002:   Returning empty structure if MINSIZE is larger than any object.
% February 2002:  Using the new function 'dip_getmeasurefeatures' - removed static lists.
% February 2004:  Improved help (KvW).
% February 2004:  Fixed bug involving min and max size. Default measurement is 'size'.
% 18 April 2005:  Exit gracefully in absence of any objects when a size
%                    restriction is imposed (MvG)
% 12 Feb 2007:    Better way of removing Mike's private measurement functions.
% February 2008:  Added strucutre for non-diplib measurements;
%                    now: extent of 2/3D cubes and elliposids. (BR)
% 6 March 2008:   Simplified removing of small objects from object ID list.
%                 Using new overloaded method RMFIELD.
% 20 Oct 2008:    Added 'MajorAxes' measurement (BR)
% 7 Dec 2009:     Added a note to the help about chain codes and connectivity.
% 14 Sep 2010:    Changed default connectivity to match label default (BR)
% 1 October 2010: Moved sub-function MAPALIASES to PRIVATE/DI_MAPALIASES.
% 19 May 2015:    Moved DI_DERIVEDMEASUREMENTS into here, simplified by using function
%                     handles, added private '-getfeatures' option, added 'axes' and
%                     'units' to non-diplib measurements.

function data = measure(varargin)

msmts = dip_getmeasurefeatures;

% Remove private elements from list - This is Michael's useless creation.
[tmp,I] = intersect({msmts.name},{'BendingEnergy','CCLongestRun','Orientation2D','Anisotropy2D'});
if ~isempty(I)
   msmts(I) = [];
end
% Add non-diplib measurements
add_msmts = di_derivedmeasurements;
msmts = cat(2,msmts,rmfield(add_msmts,{'featureID','conv_func'})');

d = struct('menu','Analysis',...
           'display','Measure',...
           'inparams',struct('name',       {'object_in',   'gray_in',         'measurementID','objectIDs', 'connectivity','minSize',            'maxSize'},...
                             'description',{'Object image','Grey-value image','Measurement',  'Object IDs','Connectivity','Minumum object size','Maximum object size'},...
                             'type',       {'image',       'image',           'optionarray',  'array',     'array',       'array',              'array'},...
                             'dim_check',  {0,             0,                 1,              -1,          0,             0,                    0},...
                             'range_check',{[],            [],                msmts,          'N+',        'N+',          'N',                  'N'},...
                             'required',   {1,             0,                 0,              0,           0,             0,                    0},...
                             'default',    {'a',           '[]',              'size',         [],          inf,             0,                    0}...
                            ),...
           'outparams',struct('name',{'msr'},...
                              'description',{'Output measurement data'},...
                              'type',{'measurement'}...
                             )...
           );
if nargin == 1
   s = varargin{1};
   if ischar(s)
      switch(s)
         case 'DIP_GetParamList'
            data = d;
            return
         case '-getfeatures'
            data =  msmts;
            return
      end
   end
end
% Aliases for elements in the 'msmts' list (backwards compatability).
if nargin>=3
   if ischar(varargin{3})
      varargin{3} = {varargin{3}};
   end
   if iscellstr(varargin{3})
      for ii=1:prod(size(varargin{3}))
         varargin{3}{ii} = di_mapaliases(varargin{3}{ii});
      end
   end
end

try
   [object_in,gray_in,measurementID,objectIDs,connectivity,minSize,maxSize] = getparams(d,varargin{:});
catch
   if ~isempty(paramerror)
      error(paramerror)
   else
      error(firsterr)
   end
end
if isinf(connectivity)
   connectivity = ndims(object_in);
end

% Check for non-DIPlib measurements
if ~iscell(measurementID)
   measurementID = {measurementID};
end
orgmeasurementID = measurementID;
added_msrID = add_msmts([]);
for ii=1:size(add_msmts,1)
   jj = find(strcmpi(measurementID,add_msmts(ii).name));
   if jj
      measurementID{jj} = add_msmts(ii).featureID;
      added_msrID(end+1) = add_msmts(ii);
   end
end
measurementID = unique(lower(measurementID));

if islogical(object_in)
   % meaning it is a binary image...
   object_in = dip_label(object_in,connectivity,'threshold_on_size',...
      minSize,maxSize,'');
   if max(object_in) < 1
      data = dip_measurement;
      return
   end
elseif minSize ~= 0 | maxSize ~= 0
   % measure object size and select objects in range
   % We don't need to do this if we just labelled the image
   if minSize ~= 0 & maxSize ~= 0 & minSize > maxSize
      error('maxSize must be larger than minSize.')
   end
   m = dip_measure(object_in,gray_in,'size',objectIDs,connectivity);
   % 18-04-2005 - MvG - Exit gracefully in absence of any objects
   if isempty( m )
      data = dip_measurement;
      return
   end
   sz = m.size;
   I = sz==sz;   % same as true(size(sz)), but works on older MATLABs.
   if minSize ~= 0
      I(sz<minSize) = 0;
   end
   if maxSize ~= 0
      I(sz>maxSize) = 0;
   end
   objectIDs = m.id(I);
   if length(objectIDs) < 1
      data = dip_measurement;
      return
   end
elseif max(object_in) < 1
   data = dip_measurement;
   return
end
data = dip_measure(object_in,gray_in,measurementID,objectIDs,connectivity);

% Was a non-DIPlib measurements requested??
% if yes, add it to the dip_measurement object
if ~isempty(added_msrID)
   for ii=1:length(added_msrID)
      %data.(added_msrID(ii).name) = added_msrID(ii).conv_func(data);
      [d,a,u] = added_msrID(ii).conv_func(data);
      tmp = [];
      tmp.id = data.id;
      tmp.data = {d};
      tmp.names = {added_msrID(ii).name};
      tmp.axes = {a};
      tmp.units = {u};
      tmp = dip_measurement('trust_me',tmp);
      data = [data,tmp];
   end
   % remove added measures from list if neccessary
   tmp = setdiff(lower(measurementID),lower(orgmeasurementID));
   for ii=1:length(tmp)
      data = rmfield(data,tmp{ii});
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = di_derivedmeasurements
out = { ...
   'DimensionsCube','extent along the principal axes of a cube','Inertia',@feat_dimensionscube;...
   'GreyDimensionsCube','extent along the principal axes of a cube (grey-weighted) *','GreyInertia',@feat_dimensionscube_grey;...
   'DimensionsEllipsoid','extent along the principal axes of an ellipsoid','Inertia',@feat_dimensionsellipsoid;...
   'GreyDimensionsEllipsoid','extent along the principal axes of an elliposid (grey-weighted)*','GreyInertia',@feat_dimensionsellipsoid_grey;...
   'MajorAxes','principal axes of an object','Mu',@feat_majoraxes;...
   'GreyMajorAxes','principal axes of an object (grey-weighted) *','GreyMu',@feat_majoraxes_grey};
out = cell2struct(out,{'name','description','featureID','conv_func'},2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [data,axes,units] = feat_dimensionsellipsoid(data)
J = flipud(data.Inertia);
units = di_convertphysDims(data.units.Inertia,'-',1);
switch size(J,1)
    case 2
        m = [0 16; 16 0];
    case 3
        m = 10.*[-1 1 1; 1 -1 1; 1 1 -1];
    otherwise
        error('Not implemented for higher dimensions than 3.');
end
data = sqrt(m*J);
axes = {'axis1';'axis2';'axis3'};
axes = axes(1:size(data,1));

function [data,axes,units] = feat_dimensionsellipsoid_grey(data)
J = flipud(data.GreyInertia);
units = di_convertphysDims(data.units.GreyInertia,'-',1);
switch size(J,1)
    case 2
        m = [0 16; 16 0];
    case 3
        m = 10.*[-1 1 1; 1 -1 1; 1 1 -1];
    otherwise
        error('Not implemented for higher dimensions than 3.');
end
data = sqrt(m*J);
axes = {'axis1';'axis2';'axis3'};
axes = axes(1:size(data,1));

function [data,axes,units] = feat_dimensionscube(data)
J = flipud(data.Inertia);
units = di_convertphysDims(data.units.Inertia,'-',1);
switch size(J,1)
    case 2
        m = [0 12; 12 0];
    case 3
        m = 6.*[-1 1 1; 1 -1 1; 1 1 -1];
    otherwise
        error('Not implemented for higher dimensions than 3.');
end
data = sqrt(m*J);
axes = {'axis1';'axis2';'axis3'};
axes = axes(1:size(data,1));

function [data,axes,units] = feat_dimensionscube_grey(data)
J = flipud(data.GreyInertia);
units = di_convertphysDims(data.units.GreyInertia,'-',1);
switch size(J,1)
    case 2
        m = [0 12; 12 0];
    case 3
        m = 6.*[-1 1 1; 1 -1 1; 1 1 -1];
    otherwise
        error('Not implemented for higher dimensions than 3.');
end
data = sqrt(m*J);
axes = {'axis1';'axis2';'axis3'};
axes = axes(1:size(data,1));

function [data,axes,units]=feat_majoraxes_grey(data)
J = data.GreyMu;
units = di_convertphysDims(data.units.GreyMu(1),'-',1);
switch size(J,1)
    case 3 %2D
        N = 2;
        axes = {'V1_x';'V1_y';'V2_x';'V2_y'};
    case 6 %3D
        N = 3;
        axes = {'V1_x';'V1_y';'V1_z';'V2_x';'V2_y';'V2_z';'V3_x';'V3_y';'V3_z'};
    otherwise
        error('should not happen.')
end
data = zeros(N*N,size(J,2));
for objectnumber = 1:size(J,2)
    j = zeros(N,N);
    l=0;
    for ii = 1:N
        for jj=ii:N
            l=l+1;
            j(ii,jj) = J(l,objectnumber);
            j(jj,ii) = J(l,objectnumber);
        end
    end
    %if any(isnan(j)) || any(isinf(j)); continue;end
    [tmp,bla] = eig(j);
    data(:,objectnumber) = tmp(:);
end
units = repmat(units,size(data,1),1);

function [data,axes,units]=feat_majoraxes(data)
J = data.Mu;
units = di_convertphysDims(data.units.Mu(1),'-',1);
switch size(J,1)
    case 3 %2D
        N = 2;
        axes = {'V1_x';'V1_y';'V2_x';'V2_y'};
    case 6 %3D
        N = 3;
        axes = {'V1_x';'V1_y';'V1_z';'V2_x';'V2_y';'V2_z';'V3_x';'V3_y';'V3_z'};
    otherwise
        error('should not happen.')
end
data = zeros(N*N,size(J,2));
for objectnumber = 1:size(J,2)
    j = zeros(N,N);
    l=0;
    for ii = 1:N
        for jj=ii:N
            l=l+1;
            j(ii,jj) = J(l,objectnumber);
            j(jj,ii) = J(l,objectnumber);
        end
    end
    %if any(isnan(j)) || any(isinf(j)); continue;end
    [tmp,bla] = eig(j);
    data(:,objectnumber) = tmp(:);
end
units = repmat(units,size(data,1),1);
