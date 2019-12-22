%%% Simulates the wave propagation through PDMS diffuser using beam propagation model.
%%% By manually changing the depth (variable-> d_nom) of the desired focus, we can simlate the
%%% correction wavefront : Edited by Abhilash Thendiyammal 2019

% clc;
close all; clear all
addpath('C:\git\bpm'); % add bpm path

%% simulation starting parameters
f_depth = 100;                  % focus depth (in um) % Note: update in between measurements
opt.pixel_size = 1/3;           % grid pixel size (in um)
n_pdms =1.41;                   % PDMS Refractive index 
n_water = 1.33;                 % water refractive index
opt.lambda = 0.804;             % wavelength in vacuum (in um)
focus_angle = 37;               % focusing angle of microscope objective (in water)

% refractive index map (path way)
dirname = '\\ad.utwente.nl\TNW\BMPI\Projects\WAVEFRONTSHAPING\data\TPM\3rd gen\191123_WFScomparison_vs_depth_PDMSdiffuser\';
filename = 'RefractiveIndex1X512.mat';

%% Import 3D refractive index medium 
load([dirname,filename]);

%% crop refractive index map to illuminated part only
crop_size = round(f_depth*tand(focus_angle)/opt.pixel_size)*2; % reduced grid size
n_sample = n_sample(end/2+(-crop_size/2+1:crop_size/2),end/2+(-crop_size/2+1:crop_size/2),:);

% Define simulation axes
Nx = size(n_sample);
x = (-size(n_sample,1)/2:size(n_sample,1)/2-1)*opt.pixel_size;
y = (-size(n_sample,2)/2:size(n_sample,2)/2-1)'*opt.pixel_size;
z = reshape((1:size(n_sample,3))*opt.pixel_size,1,1,size(n_sample,3));

% define thickness of simulated corrugated surface
d_layer = size(n_sample,3)*opt.pixel_size;

%% Step 1: Create a diverging beam as a starting wavefront at bottom of simulated surface
d_start = (f_depth*n_pdms/n_water)-d_layer;
r = sqrt(x.^2 + y.^2 + d_start.^2);             % distance to focus
D = crop_size;                                  % diameter of the gaussian window (um)
E0 = exp(2.0i * pi * r / opt.lambda );          % diverging wavefront

% create source field for BPM simulation
E0 = Field(E0, opt.pixel_size, opt.lambda, 'um');
E0 = aperture(E0, 'gaussian', D/2);             % Gaussian intensity profile

%% Step 2: Propagate through the sample layer          
[E2, E3D2] = propagate(E0, n_sample, d_layer);   % propagate through sample layer
figure(1); imagesc(E3D2(end/2, :, :));          % cross section in x-z plane
figure(); imagesc(E3D2(:,:,end));               % cross section in focal plane (xy)

%% Step 3: Phase conjugate the wavefront and focus at a distance f_m1
conjE2 = conj(E2); % Conjugate field 
[E3, E3D3] = propagate(conjE2,n_water*ones(Nx(1),Nx(2),10), f_depth); %propagate distance fm1

figure(2); imagesc(E3D3(end/2, :, :));                          % cross section in x-z plane
figure(3); imagesc(E3D3(:,:,end));                              % cross section in focal plane (xy)

E3=fftshift(fft2(ifftshift(E3)));                               % FFT of the distorted focus
Sample = E3.data;                                            

clear E3D3;

%% Define pupil plane based on NA of the objective
NA= 0.8;                                                        % NA of microscope objective
kS= 1/opt.pixel_size;                                           % Spatial sampling frequency, inverse microns
dk = kS /Reduced_npts;                                          % Spacing between discrete frequency coordinates, inverse microns
% kNA=NA/opt.lambda;                                            % radius of the pupil, inverse microns
kNA=1/1.33;                                                     % for an easy way to map the simulated wavefront on to the SLM we use unity NA for clipping the correction wavefront
NAradius=round(kNA/dk);                                         % in pixels

%% Correction wavefront
CorrectionWF=Sample(end/2-NAradius:end/2+NAradius,end/2-NAradius:end/2+NAradius);
x=linspace(-1,1,2*NAradius+1);
y=linspace(-1,1,2*NAradius+1)';
masked_C=mask(x,y, true);
CorrectionWF_SLM=masked_C(:,:,1).*CorrectionWF;                 % correction wavefront masked with unity NA
SLM_Correction_Pattern=exp(-1.0i*angle(CorrectionWF_SLM));     

%% Mapping to SLM
SLMCorrection=(angle(SLM_Correction_Pattern))*(255/2/pi);       % correction wavefront in gray values
SLMCorrection=flip(SLMCorrection,1);
figure(); imagesc(SLMCorrection);                       