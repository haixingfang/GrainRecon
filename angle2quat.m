function q = angle2quat( r1, r2, r3, varargin )
%  ANGLE2QUAT Convert rotation angles to quaternion.
%   Q = ANGLE2QUAT( R1, R2, R3 ) calculates the quaternion, Q, for given,
%   R1, R2, R3.   R1 is an M array of first rotation angles.  R2 is an M
%   array of second rotation angles.  R3 is an M array of third rotation
%   angles.  Q returns an M-by-4 matrix containing M quaternions. Q has its
%   scalar number as the first column.  Rotation angles are input in radians.    
%
%   Q = ANGLE2QUAT( R1, R2, R3, S ) calculates the quaternion, Q, for a
%   given set of rotation angles, R1, R2, R3, and a specified rotation
%   sequence, S.  
%
%   The default rotation sequence is 'ZYX' where the order of rotation
%   angles for the default rotation are R1 = Z Axis Rotation, R2 = Y Axis
%   Rotation, and R3 = X Axis Rotation. 
%
%   All rotation sequences, S, are supported: 'ZYX', 'ZYZ', 'ZXY', 'ZXZ',
%   'YXZ', 'YXY', 'YZX', 'YZY', 'XYZ', 'XYX', 'XZY', and 'XZX'.
%
%   Examples:
%
%   Determine the quaternion from rotation angles:
%      yaw = 0.7854; 
%      pitch = 0.1; 
%      roll = 0;
%      q = angle2quat( yaw, pitch, roll )
%
%   Determine the quaternions from multiple rotation angles:
%      yaw = [0.7854 0.5]; 
%      pitch = [0.1 0.3]; 
%      roll = [0 0.1];
%      q = angle2quat( pitch, roll, yaw, 'YXZ' )
%
%   See also DCM2QUAT, QUAT2DCM, QUAT2ANGLE.

%   Copyright 2000-2016 The MathWorks, Inc.

narginchk(3,4);

if any(~isreal(r1) || ~isnumeric(r1))
    error(message('aero:angle2quat:isNotReal1'));
end

if any(~isreal(r2) || ~isnumeric(r2))
    error(message('aero:angle2quat:isNotReal2'));
end

if any(~isreal(r3) || ~isnumeric(r3))
    error(message('aero:angle2quat:isNotReal3'));
end

if (length(r1) ~= length(r2)) || (length(r1) ~= length(r3))
    error(message('aero:angle2quat:wrongDimension'));
end

if nargin == 3
    type = 'zyx';
else
    if ischar( varargin{1} ) || isstring(varargin{1})
        type = char(varargin{1});
    else
        error(message('aero:angle2quat:notChar'));
    end
end

angles = [r1(:) r2(:) r3(:)];

cang = cos( angles/2 );
sang = sin( angles/2 );

switch lower( type )
    case 'zyx'
        q = [ cang(:,1).*cang(:,2).*cang(:,3) + sang(:,1).*sang(:,2).*sang(:,3), ...
            cang(:,1).*cang(:,2).*sang(:,3) - sang(:,1).*sang(:,2).*cang(:,3), ...
            cang(:,1).*sang(:,2).*cang(:,3) + sang(:,1).*cang(:,2).*sang(:,3), ...
            sang(:,1).*cang(:,2).*cang(:,3) - cang(:,1).*sang(:,2).*sang(:,3)];
    case 'zyz'
        q = [ cang(:,1).*cang(:,2).*cang(:,3) - sang(:,1).*cang(:,2).*sang(:,3), ...
            cang(:,1).*sang(:,2).*sang(:,3) - sang(:,1).*sang(:,2).*cang(:,3), ...
            cang(:,1).*sang(:,2).*cang(:,3) + sang(:,1).*sang(:,2).*sang(:,3), ...
            sang(:,1).*cang(:,2).*cang(:,3) + cang(:,1).*cang(:,2).*sang(:,3)];
    case 'zxy'
        q = [ cang(:,1).*cang(:,2).*cang(:,3) - sang(:,1).*sang(:,2).*sang(:,3), ...
            cang(:,1).*sang(:,2).*cang(:,3) - sang(:,1).*cang(:,2).*sang(:,3), ...
            cang(:,1).*cang(:,2).*sang(:,3) + sang(:,1).*sang(:,2).*cang(:,3), ...
            cang(:,1).*sang(:,2).*sang(:,3) + sang(:,1).*cang(:,2).*cang(:,3)];
    case 'zxz'
        q = [ cang(:,1).*cang(:,2).*cang(:,3) - sang(:,1).*cang(:,2).*sang(:,3), ...
            cang(:,1).*sang(:,2).*cang(:,3) + sang(:,1).*sang(:,2).*sang(:,3), ...
            sang(:,1).*sang(:,2).*cang(:,3) - cang(:,1).*sang(:,2).*sang(:,3), ...
            cang(:,1).*cang(:,2).*sang(:,3) + sang(:,1).*cang(:,2).*cang(:,3)];
    case 'yxz'
        q = [ cang(:,1).*cang(:,2).*cang(:,3) + sang(:,1).*sang(:,2).*sang(:,3), ...
            cang(:,1).*sang(:,2).*cang(:,3) + sang(:,1).*cang(:,2).*sang(:,3), ...
            sang(:,1).*cang(:,2).*cang(:,3) - cang(:,1).*sang(:,2).*sang(:,3), ...
            cang(:,1).*cang(:,2).*sang(:,3) - sang(:,1).*sang(:,2).*cang(:,3)];
    case 'yxy'
        q = [ cang(:,1).*cang(:,2).*cang(:,3) - sang(:,1).*cang(:,2).*sang(:,3), ...
            cang(:,1).*sang(:,2).*cang(:,3) + sang(:,1).*sang(:,2).*sang(:,3), ...
            sang(:,1).*cang(:,2).*cang(:,3) + cang(:,1).*cang(:,2).*sang(:,3), ...
            cang(:,1).*sang(:,2).*sang(:,3) - sang(:,1).*sang(:,2).*cang(:,3)];
    case 'yzx'
        q = [ cang(:,1).*cang(:,2).*cang(:,3) - sang(:,1).*sang(:,2).*sang(:,3), ...
            cang(:,1).*cang(:,2).*sang(:,3) + sang(:,1).*sang(:,2).*cang(:,3), ...
            cang(:,1).*sang(:,2).*sang(:,3) + sang(:,1).*cang(:,2).*cang(:,3), ...
            cang(:,1).*sang(:,2).*cang(:,3) - sang(:,1).*cang(:,2).*sang(:,3)];
    case 'yzy'
        q = [ cang(:,1).*cang(:,2).*cang(:,3) - sang(:,1).*cang(:,2).*sang(:,3), ...
            sang(:,1).*sang(:,2).*cang(:,3) - cang(:,1).*sang(:,2).*sang(:,3), ...
            cang(:,1).*cang(:,2).*sang(:,3) + sang(:,1).*cang(:,2).*cang(:,3), ...
            cang(:,1).*sang(:,2).*cang(:,3) + sang(:,1).*sang(:,2).*sang(:,3)];
    case 'xyz'
        q = [ cang(:,1).*cang(:,2).*cang(:,3) - sang(:,1).*sang(:,2).*sang(:,3), ...
            cang(:,1).*sang(:,2).*sang(:,3) + sang(:,1).*cang(:,2).*cang(:,3), ...
            cang(:,1).*sang(:,2).*cang(:,3) - sang(:,1).*cang(:,2).*sang(:,3), ...
            cang(:,1).*cang(:,2).*sang(:,3) + sang(:,1).*sang(:,2).*cang(:,3)];
    case 'xyx'
        q = [ cang(:,1).*cang(:,2).*cang(:,3) - sang(:,1).*cang(:,2).*sang(:,3), ...
            cang(:,1).*cang(:,2).*sang(:,3) + sang(:,1).*cang(:,2).*cang(:,3), ...
            cang(:,1).*sang(:,2).*cang(:,3) + sang(:,1).*sang(:,2).*sang(:,3), ...
            sang(:,1).*sang(:,2).*cang(:,3) - cang(:,1).*sang(:,2).*sang(:,3)];
    case 'xzy'
        q = [ cang(:,1).*cang(:,2).*cang(:,3) + sang(:,1).*sang(:,2).*sang(:,3), ...
            sang(:,1).*cang(:,2).*cang(:,3) - cang(:,1).*sang(:,2).*sang(:,3), ...
            cang(:,1).*cang(:,2).*sang(:,3) - sang(:,1).*sang(:,2).*cang(:,3), ...
            cang(:,1).*sang(:,2).*cang(:,3) + sang(:,1).*cang(:,2).*sang(:,3)];
    case 'xzx'
        q = [ cang(:,1).*cang(:,2).*cang(:,3) - sang(:,1).*cang(:,2).*sang(:,3), ...
            cang(:,1).*cang(:,2).*sang(:,3) + sang(:,1).*cang(:,2).*cang(:,3), ...
            cang(:,1).*sang(:,2).*sang(:,3) - sang(:,1).*sang(:,2).*cang(:,3), ...
            cang(:,1).*sang(:,2).*cang(:,3) + sang(:,1).*sang(:,2).*sang(:,3)];
    otherwise
        error(message('aero:angle2quat:unknownRotation', type));
end
