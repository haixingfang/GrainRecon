% get rolling median projection
function out=Rolling_Median_rev(im,Nrolling)
% update on Aug 5, 2020
% correct on Oct 20, 2020
if Nrolling<length(im)/2
    for i=1:length(im)
%         Obj(:,:,i)=double(im{i}(:,:));
        Obj(:,:,i)=double(im{i}(:,:));
    end
    for i=1:length(im)
        RollStart=i-fix((Nrolling-1)/2);
        RollEnd=Nrolling+RollStart-1;
        if RollStart<=0
            RollIndex=[length(im)-abs(RollStart):length(im) 1:RollEnd];
        elseif RollEnd>length(im)
            RollIndex=[RollStart:length(im) 1:RollEnd-length(im)];
        else
            RollIndex=RollStart:RollEnd;
        end

        medianObj(:,:,i)=median(Obj(:,:,RollIndex),3);
        out{i}=dip_image(medianObj(:,:,i),'uint16');
    end
else
    sprintf('Error: please set the number of rolling images smaller than %d !', Nrolling)
end


