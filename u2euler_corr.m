% Euler angles (phi1, PHI, phi2) from U matrix
% The formalism follows the ID11-3DXRD specs
% Note that there are two solutions
% (phi1, PHI, phi2) AND (phi1 + pi, -PHI, phi2 + pi)
% We pick the one with phi1 in the range [-pi/2 pi/2]
%
% Fails if U(3,2) or U(2,3) = 0 e.g. then U(3,3) = ~1
% If U(3,3) ~ 1 ph1 = ph2 = atan(U(2,1)/U(1,1))/2
% In this case there is only one solution.

% corrected on April 3, 2019

function phi_vector = u2euler_corr(U)

ph(1)=acos(U(3,3));
if ph(1) < 0.0001
    ph2(1) = atan(U(2,1)/U(1,1))/2.0;
    ph1(1) = ph2(1);
    phi1_deg = ph1(1)*180/pi;
    phi_deg = ph(1)*180/pi;
    phi2_deg = ph2(1)*180/pi;
else
    % there are 2 solutions to each of the equations
    ph1(1)=atan2(-U(1,3),U(2,3));
    ph2(1)=atan2(U(3,1),U(3,2));
    if ph1(1)<0
        ph1(1)=ph1(1)+2*pi;
        ph1(2)=ph1(1)-pi;
    else
        ph1(2)=ph1(1)+pi;
    end
    if ph2(1)<0
        ph2(1)=ph2(1)+2*pi;
        ph2(2)=ph2(1)-pi;
    else
        ph2(2)=ph2(1)+pi;
    end
    
    U1=euler2u(ph1(1),ph(1),ph2(1));
    U2=euler2u(ph1(2),ph(1),ph2(1));
    U3=euler2u(ph1(1),ph(1),ph2(2));
    U4=euler2u(ph1(2),ph(1),ph2(2));
    delta=[sum(sum(abs(U1-U))) sum(sum(abs(U2-U))) sum(sum(abs(U3-U))) sum(sum(abs(U4-U)))];
    Index=find(delta==min(delta));
    if length(Index)>1
        Index=Index(1);
    end
    if Index==1
        phi1_deg = ph1(1)*180/pi;
        phi_deg = ph(1)*180/pi;
        phi2_deg = ph2(1)*180/pi;
    end
    if Index==2
        phi1_deg = ph1(2)*180/pi;
        phi_deg = ph(1)*180/pi;
        phi2_deg = ph2(1)*180/pi;
    end
    if Index==3
        phi1_deg = ph1(1)*180/pi;
        phi_deg = ph(1)*180/pi;
        phi2_deg = ph2(2)*180/pi;
    end
    if Index==4
        phi1_deg = ph1(2)*180/pi;
        phi_deg = ph(1)*180/pi;
        phi2_deg = ph2(2)*180/pi;
    end
end
phi_vector = [phi1_deg phi_deg phi2_deg];



