%WRITEAVI   Writes 3D image into an AVI file
%  Writes a 3D grey-value image into an avi movie file, where the 3rd dimension
%  (z-direction) is the time component.
%
% SYNOPSIS:
%  writeavi(image_in, filename, fps, comp)
%
% PARAMETERS:
%  image_in: image to write to file.
%  filename: string with name of file, optionally with path and extension.
%  fps:      frames per second to write
%  comp:     compression, one of:
%            'Motion JPEG': Compressed using Motion JPEG codec
%            'Uncompressed': Uncompressed
%            For MATLAB versions older than R2010b, the options are different:
%            'none': no compression, only option for non-Windows machines
%            'Indeo3','Indeo5','Cinepak','MSVC','RLE': various classical compression schemes
%
% DEFAULTS:
%  fps = 15
%  comp = 'Cinepak' under Windows, 'None' elsewhere, for MATLAB prior to release R2010b.
%         'Motion JPEG' for newer versions of MATLAB.
%
% SEE ALSO:
%  writedisplayavi

% (C) Copyright 1999-2014               Pattern Recognition Group
%     All rights reserved               Faculty of Applied Physics
%                                       Delft University of Technology
%                                       Lorentzweg 1
%                                       2628 CJ Delft
%                                       The Netherlands
%
% Bernd Rieger, May 2001.
% 11 Nov 2002:   Directly writes image data (ignoring display mode).
%                The old (display) version is now available as writedisplayavi. (MvG)
%                Now supports colour images.
% 10 March 2004: Changed os determination so that it works. (BR)
% 11 April 2005: Catch not closing avi when error occurs. (BR)
% 31 Oct 2014:   Rewrote now that AVIFILE no longer exists. (CL)

function out = writeavi(varargin)

if matlabver_ge([7,11])
   compressionopts={'Motion JPEG','Uncompressed'};
   compdef = compressionopts{1};
else
   if isunix
      compressionopts={'None'};
      compdef = compressionopts{1};
   else
      compressionopts={'None','Indeo3','Indeo5','Cinepak','MSVC','RLE'};
      compdef = compressionopts{4};
   end
end
d = struct('menu','File I/O',...
           'display','Write AVI',...
           'inparams',struct('name',       {'in',            'filename',                 'fps',               'compression'       },...
                             'description',{'Image to write','Name of the file to write','Frames per second', 'Compression method'},...
                             'type',       {'image',         'outfile',                  'array',             'option'            },...
                             'dim_check',  {0,               0,                          0,                   0                   },...
                             'range_check',{[],              '*.*',                      'R+',                compressionopts     },...
                             'required',   {1,               1,                          0,                   0                   },...
                             'default',    {'a',             '',                         15,                  compdef             }...
                            )...
          );
if nargin == 1
   s = varargin{1};
   if ischar(s) & strcmp(s,'DIP_GetParamList')
      out = d;
      return
   end
end
if matlabver_ge([7,11]) && nargin>=4
   if strcmp(varargin{4},'None')
      varargin{4} = 'Uncompressed';
   end
end
try
   [in,filename,fps,compression] = getparams(d,varargin{:});
catch
   if ~isempty(paramerror)
      error(paramerror)
   else
      error(firsterr)
   end
end

% Testing and converting input image.
sz = imsize(in);
if length(sz) ~= 3
   error('Input image must be 3D.');
end
nslices = sz(3);
col = false;
if iscolor(in)
   in = colorspace(in,'RGB');
   col = true;
end
in = dip_array(in,'uint8');

if matlabver_ge([7,11])

   % Using the VideoWriter object introduced in MATLAB R2010b (v7.11)
   switch compression
   case compressionopts{1}
      compression = 'Motion JPEG AVI';
   case compressionopts{2}
      if col
         compression = 'Uncompressed AVI';
      else
         compression = 'Grayscale AVI';
      end
   end
   videoObject = VideoWriter(filename,compression);
   videoObject.FrameRate = fps;
   videoObject.open;
   if col
      for ii=1:nslices
         writeVideo(videoObject,squeeze(in(:,:,ii,:)));
      end
   else
      for ii=1:nslices
         writeVideo(videoObject,in(:,:,ii));
      end
   end
   videoObject.close;

else

   % The old code
   aviobj = avifile(filename,'fps',fps,'compression',compression);
   if ~col
      try %this can go wrong depending on the compression and image content
         aviobj.Colormap = gray(256);
      catch
         close(aviobj);
         error(lasterr);
      end
   end
   if col
      for ii=1:nslices
         aviobj = addframe(aviobj, squeeze(in(:,:,ii,:)));
      end
   else
      for ii=1:nslices
         aviobj = addframe(aviobj, squeeze(in(:,:,ii)));
      end
   end
   aviobj = close(aviobj);

end
