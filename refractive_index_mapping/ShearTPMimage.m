% Script to shear TPM 2D images and construct 3D volume
% This script is used to shear the TPM images according to the conversion matrices.
% Image is 180 degree rotated and scaling in X and Y axis are different.

addpath('\\ad.utwente.nl\TNW\BMPI\Projects\WAVEFRONTSHAPING\data\TPM\3rd gen\191122_WFScomparison_vs_depth_PDMSdiffuser');
filename = 'file_00001.tif'; %Change this file name to shear the TPM images acquired using scan image

%% Make 3D volime image from TPM images
info = imfinfo(filename);
Nslices = size(info,1);                           % Number of 2D images

% scan image settings
zoom = 2;
numPixels= 256;

%% Use conversion matrix to shear the TPM image using imwarp function in matlab
% find conversion matrix 
resX = 512/(numPixels*zoom);                      % Resolution in X direction
resY = 512/(numPixels*zoom);                      % Resolution in Y direction
tform =affine2d([-1.529*resX 0 0; -0.005*resX -1.499*resY 0; 0 0 1]); % Resolution in X and Y changes after converting images using conversion matrix

% shear imaging using conversion matrix and combine image stacks
TPM_3D = [];
for i = 1:Nslices
    TPMImage = imread(filename,i, 'Info', info);
    TPMImage= imwarp(TPMImage,tform);
    TPMImage(TPMImage<0)=0;                                % Remove negative values(background)
    TPM_3D(:,:,i) = TPMImage;
end 

% crop stretched parts of image stack to get a square image slices
im_width = floor(min(size(TPM_3D,1),size(TPM_3D,2))/10)*10; % nearest number divisible by 10
TPM_3D = TPM_3D(floor(end/2)+(-im_width/2:im_width/2-1),floor(end/2)+(-im_width/2:im_width/2-1),:);