%% Code to stitch the two-photon images before and after correction at
% different depths. The correction is done at 13 different depths from 80
% to 320 um. Edited by: Abhilash Thendiyammal

%% set data directory
data_dirname = '';                    % add data directory here
dirname = [data_dirname,'raw_images/'];       % combined directory name

%% Load power and gain values used during the experiment
load([dirname,'info.mat']);                 % Power and gain values during the experiment
d_nom=[80:20:320];                          % Imaging depths in um
power_factor=(power_set./power_set(1)).^2;  % Normalize the power values to the first value. Power factor by which intensity increase is sqaure of this ratio.
gain_factor=(gain_set./gain_set(1));        % Normalize gain values applied to first value (Intensity is linear to gain)
gain=power_factor.*gain_factor;             % We define a gain which is the comination of power factor and gain factor
datasets={'00001','00002','00003'};         % 00001,00002 and 00003 corresponds to Reference, Feedback_based and Modelbased TPM Data

%% Parameters for converting TPM frames to correct dimensions in um
zoom = 30;                                                 % zoom factor from TPM scan image aquisition  
numPixels= 256;                                            % Pixels in TPM frame
resX = 512/(numPixels*zoom);                               % base x resolution of the image
resY = 512/(numPixels*zoom);                               % base y resolution
tform =affine2d([-1.529 0 0; -0.005 -1.499 0; 0 0 1]);     % Coversion matrix to convert the TPM frame to original size in um

%% Make 3D images with the TPM stacks with proper normalization of the power and gain values 
% Do this once and save the stiched files into a folder
for i_set = 1:numel(datasets)  % loop through the reference, feedback-based and model-based image datasets
    dataset = datasets{i_set};
    for k=1:numel(d_nom)        % loop through all 3D image substacks
        filename=['d',num2str(d_nom(k),'%.3d'),'um_',dataset,'.tif'];
        info = imfinfo([dirname,filename]);
        Nframes =size(info,1);
        for i_frame = 1:Nframes % go through all individual image slices
            TPMImage = imread([dirname,filename],i_frame, 'Info', info);
            TPMImage= imwarp(TPMImage,tform);
            TPM_3D0(:,:,i_frame) = TPMImage./gain(k);                 % Normalized frames
        end
        TPM_3D=flip(TPM_3D0,3);                                 % Reverese the sample order to be in correct order with sample configuration
        
        % create new variable for each substack
        v = genvarname(['TPM_3D',num2str(k,'%d')]);
        eval([ v ' = TPM_3D;']);
    end
    
    % Generate stitched files for 3 different set of data. Enable one and comment rest during stitching.
    if strcmp(dataset,'00001')
        TPM3Dref=cat(3,TPM_3D1,TPM_3D2,TPM_3D3,TPM_3D4,TPM_3D5,TPM_3D6,TPM_3D7,TPM_3D8,TPM_3D9,TPM_3D10,TPM_3D11,TPM_3D12,TPM_3D13);
    elseif strcmp(dataset,'00002')
        TPM3Dfeedback=cat(3,TPM_3D1,TPM_3D2,TPM_3D3,TPM_3D4,TPM_3D5,TPM_3D6,TPM_3D7,TPM_3D8,TPM_3D9,TPM_3D10,TPM_3D11,TPM_3D12,TPM_3D13);
    elseif strcmp(dataset,'00003')
        TPM3Dmodel=cat(3,TPM_3D1,TPM_3D2,TPM_3D3,TPM_3D4,TPM_3D5,TPM_3D6,TPM_3D7,TPM_3D8,TPM_3D9,TPM_3D10,TPM_3D11,TPM_3D12,TPM_3D13);
    end
end


