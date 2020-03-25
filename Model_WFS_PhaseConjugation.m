%%% Simulates the wave propagation through PDMS diffuser using beam propagation model.
%%% By manually changing the depth (variable-> f_depth) of the desired focus, we can simulate the
%%% correction wavefront : Edited by Abhilash Thendiyammal 2019

% clc;
% close all; clear all
addpath('C:\git\bpm');          % add bpm path

%% simulation starting parameters
f_depth = 300;                  % focus depth (in um) % Note: update in between measurements
opt.pixel_size = 1/3;           % grid pixel size (in um)
n_pdms =1.41;                   % PDMS Refractive index 
n_water = 1.33;                 % water refractive index
opt.lambda = 0.804/n_water;     % wavelength in water (in um)
focus_angle = 37;               % focusing angle of microscope objective (in water)

% refractive index map (path way)
% dirname = 'Datapath';         % add data path
% filename = 'SampleProperties.mat';

%% Import 3D refractive index medium 
% load([dirname,filename]);

%% crop refractive index map to illuminated part only
crop_size = round(f_depth*tand(focus_angle)/opt.pixel_size)*2; % reduced grid size
n_sample_sim = n_sample(end/2+(-crop_size/2+1:crop_size/2),end/2+(-crop_size/2+1:crop_size/2),:)./n_water;
Nx = size(n_sample_sim);

% define thickness of simulated corrugated surface
d_layer = size(n_sample_sim,3)*opt.pixel_size + PDMS_thickness*(n_pdms/n_water-1);

%% Step 1: Create a diverging beam as a starting wavefront at bottom of simulated surface
d_start = f_depth-d_layer;                      % starting distance from focus to botom simulation volume

% create source field for BPM simulation
E1= ones(crop_size);                            % E: input field, plane wave here
E1 = Field(E1, opt.pixel_size, opt.lambda, 'um');
E1= lens(E1,-1*d_start);                        % diverging wavefront
E1 = aperture(E1, 'gaussian', crop_size/2);     % Gaussian intensity profile

%% Step 2: Propagate through the sample layer          
E2 = propagate(E1, n_sample_sim, d_layer);      % propagate through sample layer
% figure(1); imagesc(E2_3D(end/2, :, :));       % cross section in x-z plane
% figure(); imagesc(E2_3D(:,:,end));            % cross section in focal plane (xy)

%% Step 3: Phase conjugate the wavefront and focus at a distance f_m1
E3 = propagate(conj(E2),ones(Nx(1),Nx(2),10), f_depth);           % propagate back to a distance f_depth
% figure(2); imagesc(E3_3D(end/2, :, :));                         % cross section in x-z plane
% figure(3); imagesc(E3_3D(:,:,end));                             % cross section in focal plane (xy)

%% Step 4: Fourier transform to obtain wavefront at the SLM plane
E4 = fftshift(fft2(ifftshift(E3)));                               % FFT to propagate to pupil plane
E4 = E4.data;                                            

%% Define pupil plane based on NA of the objective
kS= 1/opt.pixel_size;                                             % Spatial sampling frequency, inverse microns
dk = kS /crop_size;                                               % Spacing between discrete frequency coordinates, inverse microns
kNA=1;                                                            % for an easy way to map the simulated wavefront on to the SLM we use unity kNA for clipping the correction wavefront
NAradius=round(kNA/dk);                                           % in pixels

% create circular SLM mask mimicking back aperture
x = linspace(-1,1,2*NAradius+1);
y = linspace(-1,1,2*NAradius+1)';
mask2D = x.^2 + y.^2 <= kNA;

%% Calculate correction wavefront
CorrectionWF=E4(end/2-NAradius:end/2+NAradius,end/2-NAradius:end/2+NAradius);
CorrectionWF_SLM= mask2D.*CorrectionWF;                           % correction wavefront masked with unity kNA

%% Mapping to SLM
SLMCorrection=angle(CorrectionWF_SLM)*(256/(2*pi));               % correction wavefront in gray values
SLMCorrection=flip(SLMCorrection,1);                              % pattern flipped because of 4f-system between objective and SLM
figure(); imagesc(SLMCorrection);  
%% save data
% dirname = '\\ad.utwente.nl\TNW\BMPI\Projects\WAVEFRONTSHAPING\data\TPM\3rd gen\200227_ModelWFS_daniel_merle\';
% filename = ['d',num2str(f_depth,'%.3d'),'um_model1.mat'];
% save([dirname,filename],'SLMCorrection','f_depth','CorrectionWF_SLM');