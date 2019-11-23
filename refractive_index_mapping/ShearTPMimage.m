% Script to shear TPM 2D images and construct 3D volume
% This script is used to shear the TPM images according to the conversion matrices.
% Image is 180 degree rotated and scaling in X and Y axis are different.

addpath('\\ad.utwente.nl\TNW\BMPI\Projects\WAVEFRONTSHAPING\data\TPM\3rd gen\191122_WFScomparison_vs_depth_PDMSdiffuser');
filename = 'file_00001.tif'; %Change this file name to shear the TPM images acquired using scan image

%% Make 3D volime image from TPM images
info = imfinfo(filename);
N = size(info,1);                                           % Number of 2D images
zoom = 2;
numPixels= 256;
resX = 512/(numPixels*zoom)                                 % Resolution in X direction
resY = 512/(numPixels*zoom)                                 % Resolution in Y direction

%% Use conversion matrix to shear the TPM image using imwarp function in matlab
tform =affine2d([-1.529*resX 0 0; -0.005*resX -1.499*resY 0; 0 0 1]); % Resolution in X and Y changes after converting images using conversion matrix

for i = 1:N
    TPMImage = imread(filename,i, 'Info', info);
    TPMImage= imwarp(TPMImage,tform);
    TPMImage(TPMImage<0)=0;                                           % Remove negative values(background)
    TPM_3D(:,:,i) = TPMImage;
%     imagesc(TPMImage);
%     colormap('hot');colorbar
%     drawnow
end 

%% Choose a desired side length of the TPM images
dnom=150;                                                              % How deep to focus
Reduced_FieldSize=round(dnom*tand(37)*2);                              % Field size (Distance between focus and interface*1.33(WATER)/1.51(GP))+Distance between SLM and interface)*tand(37)).
Reduced_FieldSize= Reduced_FieldSize-mod(Reduced_FieldSize,2);         % Round the value to nearest even number

TPM_width = round(Reduced_FieldSize/10)*10+10;                         % Round the value to nearest number divisible by 10
TPM_sub3D =TPM_3D(round(end/2)-TPM_width/2:round(end/2)+TPM_width/2-1, round(end/2)-TPM_width/2:round(end/2)+TPM_width/2-1, :);

%% create the Stack_3D

for i = 1:N
    currentImage = TPM_sub3D(:,:,i);
    currentImage=squeeze(mean(reshape(mean(reshape(currentImage,5,[])),size(currentImage,1)/5,5,[]),2)); %%Average and squeeze the image so that uniform intensity can be generated.
    currentImage(currentImage<0)=0;           %Remove negative values(background)
    Stack_3D(:,:,i) = currentImage;
%     imagesc(currentImage);
%     colormap('hot');colorbar
%     drawnow
end 

clearvars -except Stack_3D N hSI hSICtl 