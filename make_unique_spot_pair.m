% make the spots pairing as unique, i.e. one-to-one for simu and exp spots
% Oct 29, 2021
function SpotsPair_unique=make_unique_spot_pair(SpotsPair)

SpotsPair_unique=[];
SpotsPair_cmp=SpotsPair;
stop_flag=0;
while stop_flag~=1
    if length(SpotsPair_cmp(:,1))>=2
        SpotsPair_rest=SpotsPair_cmp(2:end,:);
        ind=find(SpotsPair_rest(:,12)==SpotsPair_cmp(1,12) & SpotsPair_rest(:,15)==SpotsPair_cmp(1,15));
        if isempty(ind)
            SpotsPair_unique=[SpotsPair_unique;SpotsPair_cmp(1,:)];
            SpotsPair_cmp=SpotsPair_cmp(2:end,:);
        else
            if SpotsPair_rest(ind,20)<SpotsPair_cmp(1,20)
                SpotsPair_unique=[SpotsPair_unique;SpotsPair_rest(ind,:)];
            else
                SpotsPair_unique=[SpotsPair_unique;SpotsPair_cmp(1,:)];
            end
            SpotsPair_cmp=setdiff(SpotsPair_rest,SpotsPair_rest(ind,:),'rows');
        end
    else
        if ~isempty(SpotsPair_cmp)
            SpotsPair_unique=[SpotsPair_unique;SpotsPair_cmp(1,:)];
        end
        stop_flag=1;
    end
end










