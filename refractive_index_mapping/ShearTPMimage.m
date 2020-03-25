% Script to shear TPM 2D images and construct 3D volume
% This script is used to shear the TPM images according to the conversion matrices.
% Image is 180 degree rotated and scaling in X and Y axis are different.

addpath('Z:Org\TNW\BMPI\Projects\WAVEFRONTSHAPING\data\TPM\3rd gen\200310_OoC_collagenimaging');
filename = '512x70s_5um_OoC_00001.tif'; %Change this file name to shear the TPM images acquired using scan image

%% Settings
im_width = 750;         % specified crop size of TPM image (in um) 

% scan image settings
zoom = 1;
numPixels= 256;

%% Make 3D volime image from TPM images
info = imfinfo(filename);
Nslices = size(info,1);                           % Number of 2D images

%% Use conversion matrix to shear the TPM image using imwarp function in matlab
% find conversion matrix 
resX = 512/(numPixels*zoom);                      % Resolution in X direction
resY = 512/(numPixels*zoom);                      % Resolution in Y direction
tform =affine2d([-1.529*resX 0 0; -0.005*resX -1.499*resY 0; 0 0 1]); % Resolution in X and Y changes after converting images using conversion matrix

% shear imaging using conversion matrix and combine image stacks
TPM_3D = zeros(im_width,im_width,Nslices);  % corrected image slices (in um)
for i = 1:Nslices
    TPMimage = imread(filename,i, 'Info', info);
    TPMimage= imwarp(TPMimage,tform);
    TPMimage(TPMimage<0)=0;                                % Remove negative values(background)
    TPMimage = TPMimage(floor(end/2)+(-im_width/2:im_width/2-1), ...
                        floor(end/2)+(-im_width/2:im_width/2-1),:); % crop TPM image to specified image size
    TPM_3D(:,:,i) = TPMimage;
end 