% X1, X2 are lists of euler angles, N*3 matrices
% x contain three angles (degree)
function ErrMean=Su_fitting(x,X1,X2)

Su_angle=x;
Su=euler2u(Su_angle(1)*pi/180,Su_angle(2)*pi/180,Su_angle(3)*pi/180);
Err=[];
for j=1:length(X2(:,1))
    E1=X1(j,:);
    E2=X2(j,:);
    U1=euler2u(E1(1)*pi/180,E1(2)*pi/180,E1(3)*pi/180);
    U1=Su*U1;
    U2=euler2u(E2(1)*pi/180,E2(2)*pi/180,E2(3)*pi/180);
    [ang,ax]=misori2(U1,U2);
    Err=[Err;ang]; % [deg]
end
ErrMean=mean(Err); % minimize the misorientation [deg]