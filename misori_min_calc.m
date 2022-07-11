function misori_min=misori_min_calc(ori)
misori_min=zeros(length(ori),3);
tic
for i=1:5%length(ori)
    for j=1:length(ori)
        if j~=i
            misori(j)=angle(ori(i),ori(j))./degree;
        end
    end
    [misori_minimum,ind]=min(nonzeros(misori));
    misori_min(i,:)=[misori_minimum i ind];
    i
end
toc