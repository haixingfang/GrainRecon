function [ang, ax] = misori2(ori1, ori2, opt)

ori3 = ori1';
orit = ori2*ori3;
ang = 180; %Initial misorientation angle

i = 0;
while 1
    i = i+1;
    ori = eqmf(i)*orit;
    temp = (trace(ori)-1)/2;
    if temp <= -0.9999999999
        temp = 180;
    elseif temp >= 0.9999999999
        ang = 0;
        ax = [1,0,0];
        return
    else
        temp = acosd(temp);
        if temp <= 45 %In theory, correct!
            ang = temp;
            ax1 = ori(2,3)-ori(3,2);
            ax2 = ori(3,1)-ori(1,3);
            ax3 = ori(1,2)-ori(2,1);
            ax = [ax1 ax2 ax3];
            nax = norm(ax);
            if abs(nax) >= 0.000001
                ax = ax./nax;
            end
            if nargin==3
                if opt==1
                    ax=integer_vec(ax);
                end
            end
            return
        end
    end
    if temp < ang
        ang = temp;
        Index = i;
    end

    if i == 24 %Exit loop
        break
    end
end

if Index >0 && Index <= 24
    ori = eqmf(Index)*orit;
    ax1 = ori(2,3)-ori(3,2);
    ax2 = ori(3,1)-ori(1,3);
    ax3 = ori(1,2)-ori(2,1);
    ax = [ax1 ax2 ax3];
    nax = norm(ax);
    if abs(nax) >= 0.000001
        ax = ax./nax;
    end
    if nargin==3
        if opt==1
            ax=integer_vec(ax);
        end
    end
else
    ang = 180;
    ax = [1 0 0];
end

%equivalent orientation matrix
function eqm = eqmf(i)
switch i
    case 1
        eqm=[1 0 0;0 1 0;0 0 1];
    case 2
        eqm=[1 0 0;0 -1 0;0 0 -1];
    case 3
        eqm=[1 0 0;0 0 -1;0 1 0];
    case 4
        eqm=[1 0 0;0 0 1; 0 -1 0];
    case 5
        eqm=[-1 0 0;0 1 0;0 0 -1];
    case 6
        eqm=[-1 0 0;0 -1 0;0 0 1];
    case 7
        eqm=[-1 0 0;0 0 -1;0 -1 0];
    case 8
        eqm=[-1 0 0;0 0 1;0 1 0];
    case 9
        eqm=[0 1 0;-1 0 0;0 0 1];
    case 10
        eqm=[0 1 0;0 0 -1;-1 0 0];
    case 11
        eqm=[0 1 0;1 0 0;0 0 -1];
    case 12
        eqm=[0 1 0;0 0 1;1 0 0];
    case 13
        eqm=[0 -1 0;1 0 0;0 0 1];
    case 14
        eqm=[0 -1 0;0 0 -1;1 0 0];
    case 15
        eqm=[0 -1 0;-1 0 0;0 0 -1];
    case 16
        eqm=[0 -1 0;0 0 1;-1 0 0];
    case 17
        eqm=[0 0 1;0 1 0;-1 0 0];
    case 18
        eqm=[0 0 1;1 0 0;0 1 0];
    case 19
        eqm=[0 0 1;0 -1 0;1 0 0];
    case 20
        eqm=[0 0 1;-1 0 0;0 -1 0];
    case 21
        eqm=[0 0 -1;0 1 0;1 0 0];
    case 22
        eqm=[0 0 -1;-1 0 0;0 1 0];
    case 23
        eqm=[0 0 -1;0 -1 0;-1 0 0];
    case 24
        eqm=[0 0 -1;1 0 0;0 -1 0];
end