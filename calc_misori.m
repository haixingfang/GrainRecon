
function [ang,ax,ang_mtex]=calc_misori(ori1,ori2,orientationType)

cs = crystalSymmetry('cubic');
if nargin<3
    sprintf('Assum input as euler angles')
    U1=euler2u(ori1(1)*pi/180,ori1(2)*pi/180,ori1(3)*pi/180);
    U2=euler2u(ori2(1)*pi/180,ori2(2)*pi/180,ori2(3)*pi/180);
    [ang,ax]=misori2(U1,U2);
    rot1 = rotation('Euler',ori1(:,1)*degree,ori1(:,2)*degree,ori1(:,3)*degree);
    rot2 = rotation('Euler',ori2(:,1)*degree,ori2(:,2)*degree,ori2(:,3)*degree);
    o1=orientation(rot1,cs);
    o2=orientation(rot2,cs);
    ang_mtex=angle(o1,o2)/degree;
else
    switch orientationType
        case 'euler'
            U1=euler2u(ori1(1)*pi/180,ori1(2)*pi/180,ori1(3)*pi/180);
            U2=euler2u(ori2(1)*pi/180,ori2(2)*pi/180,ori2(3)*pi/180);
            [ang,ax]=misori2(U1,U2);
            rot1 = rotation('Euler',ori1(:,1)*degree,ori1(:,2)*degree,ori1(:,3)*degree);
            rot2 = rotation('Euler',ori2(:,1)*degree,ori2(:,2)*degree,ori2(:,3)*degree);
            o1=orientation(rot1,cs);
            o2=orientation(rot2,cs);
            ang_mtex=angle(o1,o2)/degree;
        case 'matrix'
            [ang,ax]=misori2(ori1,ori2);
            rot1 = rotation('matrix',ori1);
            rot2 = rotation('matrix',ori2);
            o1=orientation(rot1,cs);
            o2=orientation(rot2,cs);
            ang_mtex=angle(o1,o2)/degree;
        case 'rodrigues'
            x1=rod2euler(ori1);
            x2=rod2euler(ori2);
            U1=euler2u(x1(1)*pi/180,x1(2)*pi/180,x1(3)*pi/180);
            U2=euler2u(x2(1)*pi/180,x2(2)*pi/180,x2(3)*pi/180);
            [ang,ax]=misori2(U1,U2);
            ang_mtex=ang;
        case 'quaternion'
            U1=quaternion2U(ori1);
            U2=quaternion2U(ori2);
            [ang,ax]=misori2(U1,U2);
            ang_mtex=ang;
        otherwise
            sprintf('Assum input as euler angles')
            U1=euler2u(ori1(1)*pi/180,ori1(2)*pi/180,ori1(3)*pi/180);
            U2=euler2u(ori2(1)*pi/180,ori2(2)*pi/180,ori2(3)*pi/180);
            [ang,ax]=misori2(U1,U2);
            rot1 = rotation('Euler',ori1(:,1)*degree,ori1(:,2)*degree,ori1(:,3)*degree);
            rot2 = rotation('Euler',ori2(:,1)*degree,ori2(:,2)*degree,ori2(:,3)*degree);
            o1=orientation(rot1,cs);
            o2=orientation(rot2,cs);
            ang_mtex=angle(o1,o2)/degree;
    end
end
end


