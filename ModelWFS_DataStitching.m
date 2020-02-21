%% Code to stitch the two-photon images before and after correction at
% different depths. The correction is done at 13 different depths from 80
% to 320 um.  

addpath('C:\git\utilities');
addpath('P:\TNW\BMPI\Users\Abhilash Thendiyammal\Research@UT\Data\191223_WFScomparison_vs_depth_PDMSdiffuser\');

%% Load power and gain values used during the experiment
load info_power;                                           % Power and gain values during the experiment
d_nom=[80:20:320];                                         % Imaging depths in um
power_factor=(power_set./power_set(1)).^2;                 % Normalize the power values to the first value. Power factor by which intensity increase is sqaure of this ratio.
gain_factor=(gain_set./gain_set(1));                       % Normalize gain values applied to first value (Intensity is linear to gain)
gain=power_factor.*gain_factor;                            % We define a gain which is the comination of power factor and gain factor
filename=['d',num2str(d_nom(1),'%.3d'),'um_00003','.tif']; % 00001,00002 and 00003 corresponds to Reference, Feedback_based and Modelbased TPM Data
info = imfinfo(filename);

%% Parameters for converting TPM frames to correct dimensions in um
zoom = 30;                                                 % zoom factor from TPM scan image aquisition  
numPixels= 256;                                            % Pixels in TPM frame
resX = 512/(numPixels*zoom)                                % base x resolution of the image
resY = 512/(numPixels*zoom)                                % base y resolution
tform =affine2d([-1.529 0 0; -0.005 -1.499 0; 0 0 1]);     % Coversion matrix to convert the TPM frame to original size in um


%% Grab first image to find the dimension of the images
TPMImage0 = imread(filename,1, 'Info', info);              
TPMImage0= imwarp(TPMImage0,tform);
x_data= -(size(TPMImage0,1)/2)*resX:resX:(size(TPMImage0,1)/2)*resX;
y_data= -(size(TPMImage0,2)/2)*resY:resY:(size(TPMImage0,2)/2)*resY;
                                       
startFrame = 1;                                           
endFrame =size(info,1);

%% Make 3D images with the TPM stacks with proper normalization of the power and gain values 
% Do this once and save the stiched files into a folder

for k=1:numel(d_nom)
filename=['d',num2str(d_nom(k),'%.3d'),'um_00003','.tif'];
info = imfinfo(filename);
        for i = startFrame:endFrame
            TPMImage = imread(filename,i, 'Info', info);
            TPMImage= imwarp(TPMImage,tform);
%             TPMImage=TPMImage-min(min(TPMImage));
%             TPMImage(TPMImage<0)=0;                             % Remove negative values(background)
%             TPM_3D0(:,:,i) = TPMImage;                           % Raw (not normalized) frames
            TPM_3D0(:,:,i) = TPMImage./gain(k);               % Noramalized frames 
%             TPM_3D0(:,:,i) = TPMImage./gain_factor(k);          % Normalized backgrounds                                                                    
%             figure, imagesc(TPMImage);  
%             colormap('hot');colorbar
%             drawnow
        end
        TPM_3D=flip(TPM_3D0,3);                                 % Reverese the sample order to be in correct order with sample configuration
        
v = genvarname(['TPM_3D',num2str(k,'%d')])
eval([ v ' = TPM_3D;'])
end

% TPM3Dref=cat(3,TPM_3D1,TPM_3D2,TPM_3D3,TPM_3D4,TPM_3D5,TPM_3D6,TPM_3D7,TPM_3D8,TPM_3D9,TPM_3D10,TPM_3D11,TPM_3D12,TPM_3D13);
% TPM3Dfeedback=cat(3,TPM_3D1,TPM_3D2,TPM_3D3,TPM_3D4,TPM_3D5,TPM_3D6,TPM_3D7,TPM_3D8,TPM_3D9,TPM_3D10,TPM_3D11,TPM_3D12,TPM_3D13);
TPM3Dmodel=cat(3,TPM_3D1,TPM_3D2,TPM_3D3,TPM_3D4,TPM_3D5,TPM_3D6,TPM_3D7,TPM_3D8,TPM_3D9,TPM_3D10,TPM_3D11,TPM_3D12,TPM_3D13);

