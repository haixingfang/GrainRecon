% DISTANCEBETWEENPOINTSETS  Compute distance and statistics thereof between sets of points
%
%   di = distancebetweenpointsets(p1,p2,<keywords>);
%
%   The coordinate arrays must be in the format NxD where N is the number
%   of coordinates and D the dimensionality of the points. The main advantage
%   of this function over directly computing the distance matrix is memory
%   usage if the distance matrix itself is not needed.
%
%   - p2 can have a different number of points to p1
%   - p2 must have the same dimensionality as p1
%   - p2 can be [] in which case p1 is compared against itself
%
%   The structure of di is as follows:
%     di.dim        with dim as below
%     di.distance   the distance matrix d(i,j)
%
%   dim is a struct array of length 2 if both p1 and p2 are given, and of
%   length 1 if only p1 is given (self-comparison). dim has the fields:
%     dim.max       : max(i) = max_j(d(i,j))    max(j) = max_i(d(i,j))
%     dim.index_max : argmax_j(d(i,j))          argmax_i(d(i,j))
%     dim.count_max : sum_j(max(i)==d(i,j))     sum_i(max(j)==d(i,j))
%     dim.min, dim.index_min, dim.count_min  as above for "min"
%     dim.mean      : mean_j(d(i,j))            mean_i(d(i,j))
%
%   comments: d(i,j) is the distance between the ith point of p1 and
%             the jth point of p2.
%             index_max/min holds the index of the first maximum/minimum
%             in case there is a tie. count_max/min holds the number of ties
%             for the first place.
%             max/min values are compared through strict equality through
%             the == operator so the usual caveats about floating point
%             comparisons do apply.
%
%   which results are actually computed is controlled through the following
%   keywords: "distance", "max", "min", "mean", "count", "index"
%   "count" and "index" are meaningless unless at least one of "max" or "min"
%   is requested as well.
%
% WARNING WARNING WARNING
%   this function can be passed both matlab arrays and dip_image's. In both
%   case the first index is taken to represent the ith point, the second index
%   the dimension j (0,1,2=x,y,z). Notably if im is a dip_image these are
%   NOT equivalent:
%     di = distancebetweenpointsets(im,[])
%     di = distancebetweenpointsets(double(im),[])
%   Also, the indices reported by DISTANCEBETWEENPOINTSETS are in the
%   awkward matlab convention of starting at one. The reason for this is that
%   in the typical usage case, the results are more likely to be used in
%   matlab functions rather than dipimage functions. The underlying DIPlib
%   function does use proper indices.
%
% Periodic space
%   There is support for points in a periodic space. Use the keywords "dims"
%   followed by an array of the period along each dimension. For example:
%   "dims", [4,2,10] for a three-dimensional space with period 4 along the
%   x-axis, 2 along the y-axis and 10 along the z-axis.
%   The distance along each axis is computed as:
%   dist=(p1-p2)-round((p1-p2)/period)*period;
%   The Euclidean distance is then computed as normal from the per-axis
%   distance.

function di=distancebetweenpointsets(p1,p2,varargin)
  if ~strcmp(class(p1),'dip_image')
    p1=p1';
  end
  if ~isempty(p2) && ~strcmp(class(p2),'dip_image')
    p2=p2';
  end
  di=dip_distancebetweenpointsets(p1,p2,varargin{:});
  fns=fieldnames(di.dim);
  for ii=1:length(di.dim)
    for fn=fns'
      di.dim(ii).(fn{1})=di.dim(ii).(fn{1})';
      if ~isempty(strfind(fn{1},'index'))
        di.dim(ii).(fn{1})=di.dim(ii).(fn{1})+1;
      end
    end
  end
  di.distance=di.distance';
end
