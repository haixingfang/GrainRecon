function [CsI Swank]=ReadCsI()
fileID=fopen('CsI.txt','r');
CsI=[];
while(~feof(fileID))
    textdata=str2num(fgetl(fileID));
    if isempty(textdata)
        continue;
    else
        CsI=[CsI; textdata];
    end
end
fclose(fileID);

fileID=fopen('Swank_factor.txt','r');
Swank=[];
while(~feof(fileID))
    textdata=str2num(fgetl(fileID));
    if isempty(textdata)
        continue;
    else
        Swank=[Swank; textdata];
    end
end
fclose(fileID);



