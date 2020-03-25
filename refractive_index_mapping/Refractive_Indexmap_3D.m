%%% Script to convert the TPM 3D  image to a refractive index
%%% distribution. This code fits the the intensity along the each z-columns
%%% to sigmoid functions and find out the point of inflection. 
%%% Refractive index values at the both sides of the point of inflection
%%% are assigned from apriori knowledge of values. Here, 1.33 and 1.41
%%% for water-diluted fluorescein and PDMS medium respectively.
%%% Edited by Abhilash Thendiyammal 2019

%%
clf
% close all
% clear all
addpath('..\..\..\utilities');

n_slices = size(TPM_3D,3);

%% Fitting parameters
d_sample = 10;                  % downsample factor used to speed up fitting of surface
u_sample = 3;                   % upsample factor from tpm data to simulation grid
pdms_slices = 1:10;             % assumption: first 10 layers always consist of PDMS
fluo_slices = n_slices+(-9:0);  % assumption: last 10 layers always consist of fluorescene
n_pdms = 1.41;                  % refractive index PDMS
n_water = 1.33;                 % refractive index water

%% In order to find the interface, fit each Z column to a sigmoid function
% Downsample TPM 3D image to speed up fitting procedure
stack_3D = imresize(TPM_3D,1/d_sample);     % downsample in x and y
n_slices = size(stack_3D,3);                % number of slices in data set
z = (1:n_slices)';                          % depth coordinates TPM data

% preallocation
z_interface = zeros(size(stack_3D,1),size(stack_3D,2)); % fitted depth cordinates for PDMS interface                      

% find interface for each intensity line along depth by fitting
starttime = now;
for x_i=1:size(stack_3D,1)                   % loop through X-coordinates of each Image
    for y_i=1:size(stack_3D,2)               % loop through Y-coordinates of each image
        % select intensity along depth for given coordinates
        Idepth=stack_3D(x_i,y_i,:);                    
        
        % normalize data from 0 to 1
        Ibg = mean(Idepth(pdms_slices));
        Imax = mean(Idepth(fluo_slices));
        Idepth_norm = (Idepth(:)-Ibg)/(Imax-Ibg);

        % perform sigmoid fit and find point of inflection 
        fit=sigmoid_fit(z,Idepth_norm);% Fit a sigmoid function to eack z-column.
        z_interface(x_i,y_i) = fit.b;
%         figure(1);  clf; plot(fit); hold on; plot(z,Idepth_norm(:),'*b'); title(fit.b); pause(0.2); ylim([-0.2,1.2]);
    end    
    eta(x_i, size(stack_3D,1), starttime, 'console', 'Fitting stuff...', 0);
end

% upscale fitted surface to simulation grid size (only x and y)
z_interface_usample = imresize(z_interface,d_sample*u_sample);    % upsample fitted surface x and y

% determine center depth of corrugation
PDMS_thickness=round(median(z_interface(:)));

%% generate refractive index map using fitted surface
% convert depth axes to simulation grid
dz = (1/u_sample);
z_sim = reshape((dz:dz:n_slices),1,1,n_slices*u_sample);

% generate 3D refractive index map known values for water and PDMS
tic;
disp('generating 3D refractive index map...');
n_sample = n_water*ones(size(z_interface_usample,1),size(z_interface_usample,2),numel(z_sim));
n_sample(z_sim <= z_interface_usample) = n_pdms;
disp('done');
toc;

%% plot test results
% plot fitted surface
figure(1); surf(z_interface); shading interp; % check surface fit

% plot example 2d slice of refractive index and compare to TPM data
figure(2);
subplot(1,2,1); imagesc(squeeze(TPM_3D(end/2,:,:))');  axis off; set(gca,'YDir','normal');
subplot(1,2,2); imagesc(squeeze(n_sample(end/2,:,:))');axis off; set(gca,'YDir','normal');
