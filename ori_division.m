% for cubic symmetry
% July 9, 2021
function OR=ori_division(nDiv)
    lims=[sqrt(2)-1 sqrt(2)-1 sqrt(2)-1]; %side of the truncated cube in Rodrigues space for a cubic geometry
    [lbr0,ubr0,rs0]=generate_constrains(nDiv);
    lbr=[lbr0(:,1) lbr0(:,2) lbr0(:,3)].*lims;
    ubr=[ubr0(:,1) ubr0(:,2) ubr0(:,3)].*lims;%ubr0*lims;
    rod=[rs0(:,1) rs0(:,2) rs0(:,3)].*lims;%rs0*lims;
    q=rod2quat(rod);
    q_b=q(:,2);
    q_c=q(:,3);
    q_d=q(:,4);
    [euler_angles(:,1), euler_angles(:,2), euler_angles(:,3)] = rod2angle(rod,'ZXZ');
    OR.len=length(q_b(:,1));
    OR.euler_angles=euler_angles;
    OR.rod=rod;
    OR.q=q;
    OR.q_b=q_b;
    OR.q_c=q_c;
    OR.q_d=q_d;
end

function [ lbr0,ubr0,rs0 ] = generate_constrains(nDiv)
    %Divides the rodrigues space in nDiv divisions
    CubesMatrix=ones(nDiv^3,3);
    for i=1:nDiv^3
        CubesMatrix(i,1)=ceil(i/nDiv^2);
        CubesMatrix(i,2)=ceil(i/nDiv)-(floor((i-1)/nDiv^2)*nDiv);
        CubesMatrix(i,3)=i-(CubesMatrix(i,2)-1)*nDiv-(CubesMatrix(i,1)-1)*nDiv^2;
    end
    rs0=((((2*CubesMatrix)-1)./(nDiv))-1);
    lbr0=rs0-(4/25);
    ubr0=rs0+(4/25);
    % scartter3(rs0(1,1,:),rs0(1,2,:),rs0(1,3,:))
    % set(gca,'XLim',[-1 1],'YLim',[-1 1],'ZLim',[-1 1])
    lbr0(lbr0<-1)=-1.01;
    ubr0(ubr0>1)=1.01;
end
