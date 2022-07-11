%Initialize DIPLIB 
%Update this file to match your installation
%
% Henning Osholm S�rensen, Ris� National Laboratory, June 23, 2006.

%{
os = computer

if strcmp(os,'PCWIN') == 1
    diproot ='C:\Program Files\DIPimage 2.2\'
    addpath(diproot,[diproot,'\dipimage'],[diproot,'\dipimage\diplib'],'-begin')
    dipstart
    %dip_initialise
else
    diproot ='/usr/local/dip/toolbox'
    addpath([diproot,'/dipimage'],[diproot,'/diplib'],'-begin')
    dip_initialise
end
%}

%%%%%%%%%%% on Feb 28, 2019 modified by Haixing Fang
% function load_diplib
%     addpath('C:\Program Files\dipimage_2.9_win64\dip\common\dipimage');
%     dip_initialise;
%     dipsetpref('imagefilepath','C:\Program Files\dipimage_2.9_win64\dip\images')
% end
%{
function load_diplib
	addpath(strcat('D:\Documents\Matlab\LabDCT_simap','\dipimage_2.9_win64\dip\common\dipimage'));
    dip_initialise;
    dipsetpref('imagefilepath',strcat('\dipimage_2.9_win64\dip\images'));
end
%}
function load_diplib
	if isunix
		addpath('./dip/Linuxa64/lib/');
		addpath('./dip/common/dipimage/');
		addpath('./dip/common/dipimage/demos');
		dip_initialise;
	elseif ispc
		addpath('./dipimage_2.9_win64/dip/common/dipimage');
        dip_initialise;
        dipsetpref('imagefilepath','./dipimage_2.9_win64/dip/images');
	else
		disp('Platform not supported')
	end
end
