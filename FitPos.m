% fit for the position
% Dec 14, 2021
function center_shift=FitPos(DS_ref,DS_rec,PairIndex)
    x0=[0 0 0]; % [pixel]
    ErrMean=@(x)FitPos_fun(x,DS_ref,DS_rec,PairIndex);
    LB=x0-20;
    UB=x0+20;
    opts=optimset('Display','iter','Algorithm','trust-region-reflective', ...
        'MaxFunEvals',2000,'MaxIter',500,'TolFun',1e-10,'PlotFcns','optimplotfval');
    [x,fval1,exitflag,output] = fminsearchbnd(ErrMean,x0,LB,UB,opts);
    center_shift=x*DS_rec.VoxSize(1); % [mm]
end

function ErrMean=FitPos_fun(x,DS_ref,DS_rec,PairIndex)
    Err=[];
    for k=1:length(PairIndex(:,1))
        i=PairIndex(k,1);
        j=PairIndex(k,2);

        COM_ref=(DS_ref.Coord(i,:)-DS_ref.Dimension/2)+DS_ref.Center'./DS_ref.VoxSize(1); % [pixel]
        COM_rec=(DS_rec.Coord(j,:)-DS_rec.Dimension/2)+DS_rec.Center'./DS_rec.VoxSize(1); % [pixel]
        COM_rec=COM_rec+x; % [pixel]
        dis=sqrt(sum((COM_ref-COM_rec).^2)); % [pixel]
        Err=[Err;dis]; % [pixel]
    end
    ErrMean=mean(Err);
end