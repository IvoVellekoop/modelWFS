%% Script to shear and combine TPM 2D intensity images to construct a 3D volume (output - TPM_3D)
% This script is used to shear the TPM images according to the conversion matrices.
% Edited by Abhilash Thendiyammal 2018

%% add path to the raw TPM image
data_dirname = '';                      % add data directory here  
dirname = [data_dirname,'raw_images\']; % pathway to TPM image of diffuser surface
filename = 'PDMS_diffuser_surface_1X512.tif';

%% Settings
im_width = 500;     % specified crop size of TPM image (in pixels) 
zoom = 1;           % zoom factor used in scan image
numPixels= 512;     % number of pixels 

%% Extract info from TPM image file from ScanImage 
info = imfinfo(filename);
Nslices = size(info,1);                           % Number of 2D image slices

%% Use conversion matrix to shear the TPM image using imwarp function in matlab
% find conversion matrix 
resX = 512/(numPixels*zoom);                      % Resolution in X direction
resY = 512/(numPixels*zoom);                      % Resolution in Y direction
tform =affine2d([-1.529*resX 0 0; -0.005*resX -1.499*resY 0; 0 0 1]); % Resolution in X and Y changes after converting images using conversion matrix

% shear images using conversion matrix and combine them
TPM_3D = zeros(im_width,im_width,Nslices);        % corrected image slices (in um)
for i = 1:Nslices
    TPMimage = imread(filename,i, 'Info', info);
    TPMimage= imwarp(TPMimage,tform);
    TPMimage(TPMimage<0)=0;                       % Remove negative values(background)
    TPMimage = TPMimage(floor(end/2)+(-im_width/2:im_width/2-1), ...
                        floor(end/2)+(-im_width/2:im_width/2-1),:); % crop TPM image to specified image size
    TPM_3D(:,:,i) = TPMimage;
end