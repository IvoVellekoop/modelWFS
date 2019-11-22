%%% Simulates the wave propagation of a converging beam with NA of the TPM objective lens through
%%% two media with a corrogated interface using beam propagation model.
%%% By manuallychanging the depth (variable-> dnom) of the desired focus, we can simlate the
%%% correction wavefront
%%% Edited by Abhilash Thendiyammal 2019
% clc;
% close all; clear all
addpath('\\ad.utwente.nl\TNW\BMPI\Users\Abhilash Thendiyammal\Research@UT\Matlab programs\git\bpm');
addpath('\\ad.utwente.nl\TNW\BMPI\Users\Abhilash Thendiyammal\Research@UT\Matlab programs\git\wavesim');
addpath('\\ad.utwente.nl\TNW\BMPI\Users\Abhilash Thendiyammal\Research@UT\Matlab programs\git\AO');

%% simulation options
opt.lambda = 0.804;                                             % wavelength in vacuum (in um)
opt.pixel_size = 0.3333;                                        % grid pixel size (in um)

%% Define beam size for beam prpagation simulation

FieldSize = 800;                                                % field size in um (arbitrary) 
k0=2*pi/opt.lambda;                                             % wavenumber of the beam
npts = round((FieldSize/opt.pixel_size));                       % number of data points in x(or y) dimension, field is a square
Xe=linspace(-FieldSize/2,FieldSize/2,npts);                     % space scale in 1D for input field
[X,Y]=meshgrid(Xe,-Xe);
E0=exp(1i*(k0*0.*X+k0*0.*Y));                                   % E: input field, plane wave here

dnom=300;                                                        % How deep to focus
Reduced_FieldSize=round(dnom*tand(37)*2);                       % Field size (Distance between focus and interface*1.33(WATER)/1.51(GP))+Distance between SLM and interface)*tand(37)).
Reduced_FieldSize= Reduced_FieldSize-mod(Reduced_FieldSize,2);  % Round the value to nearest even number
Reduced_npts=round((Reduced_FieldSize/opt.pixel_size));         % Reduced number of data points
E0_reduced=E0(size(E0,1)/2-(Reduced_npts/2-1):size(E0,1)/2+Reduced_npts/2,size(E0,1)/2-(Reduced_npts/2-1):size(E0,1)/2+Reduced_npts/2);

clear FieldSize;
clear E0;
clear npts;

%% Import 3D refractive index medium 
% Sample interface from TPM (Import refractive index generated from TPM)

load('\\ad.utwente.nl\TNW\BMPI\Projects\WAVEFRONTSHAPING\data\TPM\3rd gen\191106_POP_ModelBasedWFS\340um\Sample_RI_stretched.mat');

Xshift=0;                                                      % Shift center if required
Yshift=0;                                                      % Shift center if required

Super_Stack_3D_reduced=Super_Stack_3D((size(Super_Stack_3D,1)/2-Xshift)-Reduced_FieldSize/2:(size(Super_Stack_3D,1)/2-Xshift)+Reduced_FieldSize/2-1,(size(Super_Stack_3D,2)/2+Yshift)-Reduced_FieldSize/2:(size(Super_Stack_3D,2)/2+Yshift)+Reduced_FieldSize/2-1,:);
Nz=size(Super_Stack_3D_reduced,3);

clear Super_Stack_3D;
clear Sub_Stack_3D Xshift Yshift

%% Interpolate the TPM data to satisfy sampling condition for simulation

IntX=linspace(1,Reduced_FieldSize,Reduced_npts);IntY=linspace(1,Reduced_FieldSize,Reduced_npts);IntZ=linspace(1,Nz,Nz); % Interpolation coordinates
[IntX,IntY,IntZ]=ndgrid(IntX,IntY,IntZ);                          
n_sample=interpn(Super_Stack_3D_reduced,IntX,IntY,IntZ);          % Refractive index of the sample
clear IntX; clear IntY; clear IntZ;

% n_sample = Interpolated_Stack_3D(:,:,1:Nz);                     % Refractive index of the sample

% clear Interpolated_Stack_3D;
clear Super_Stack_3D_reduced

%% Simulate Spherical Abberation through Corrogated Interface

f_air=dnom/1.33;                                                % focal length in air (um)
f_m1 = f_air*1.33;                                              % focal length in medium 1 (um)
f_m2 = f_air*1.41;                                              % focal length in medium 2 (um)
D = Reduced_FieldSize;                                          % diameter of the gaussian window (um)

E = Field(E0_reduced, opt.pixel_size, opt.lambda, 'um');        % Electric field in Filed format (refer to BPM folder in git repository for more information)
E = aperture(E, 'gaussian', D/2);                               % Gaussian field
figure(); imagesc(angle(E));

E=lens(E,f_air);                                                % lens function              
figure(2); imagesc(angle(E));                                   % display phase of the incident light (converging wave)

for j=0:size(n_sample,3)-1
    n_sample_rev(:,:,j+1)=n_sample(:,:,size(n_sample,3)-j);     % Reverese the sample order to be in correct order with sample configuration
end

clear n_sample;

%% Calculate the sample wavefront
%simulate through 3D refractive index media for 60 um

tic();
LayerThickness=(60-PDMS_thickness)+PDMS_thickness*1.06;                                      % TPM layer thickness                        
[E1, E3D1] = propagate(E, n_sample_rev, LayerThickness); 
figure(1); imagesc(E3D1(end/2, :, :));                          %cross section in x-z plane
figure(); imagesc(E3D1(:,:,end));                               %cross section in focal plane (xy)

clear E3D1;
clear n_sample_rev;

[E2, E3D2] = propagate(E1, ones(Reduced_npts,Reduced_npts,30)*1.41, f_m2-LayerThickness); %propagate to the focus through medium (n=1.41) in 32 steps

figure(2); imagesc(E3D2(end/2, :, :));                          %cross section in x-z plane
figure(3); imagesc(E3D2(:,:,end));                              %cross section in focal plane (xy)

clear E3D2;

E3=fftshift(fft2(ifftshift(E2)));                               % FFT of the distorted focus
Sample = E3.data;                                               
toc();

%% Calculate the Reference wavefront
[E2, E3D2] = propagate(E, ones(Reduced_npts,Reduced_npts,30)*1.33, f_m1); %propagate to the focus through medium (n=1.41) in 32 steps

figure(2); imagesc(E3D2(end/2, :, :));                          %cross section in x-z plane
figure(3); imagesc(E3D2(:,:,end));                              %cross section in focal plane (xy)

E3=fftshift(fft2(ifftshift(E2)));                               % FFT of the reference focus
Reference = E3.data;                                            

clear E3D1;
clear E3D2;
clear E1 E2 E3;
%% Define pupil plane based on NA of the objective
NA= 0.8;                                                        % NA of microscope objective
kS= 1/opt.pixel_size;                                           % Spatial sampling frequency, inverse microns
dk = kS /Reduced_npts;                                          % Spacing between discrete frequency coordinates, inverse microns
% kNA=NA/opt.lambda;                                            % radius of the pupil, inverse microns
kNA=1;                                                          % for an easy way to map the simulated wavefront on to the SLM we use unity NA for clipping the correction wavefront
NAradius=round(kNA/dk);                                         % in pixels

%% Correction wavefront
CorrectionWF=Sample(end/2-NAradius:end/2+NAradius,end/2-NAradius:end/2+NAradius)./Reference(end/2-NAradius:end/2+NAradius,end/2-NAradius:end/2+NAradius);
x=linspace(-1,1,2*NAradius+1);
y=linspace(-1,1,2*NAradius+1)';
masked_C=mask(x,y, true);
CorrectionWF_SLM=masked_C(:,:,1).*CorrectionWF;                 % correction wavefront masked with unity NA
SLM_Correction_Pattern=exp(-1.0i*angle(CorrectionWF_SLM));     

%% Mapping to SLM
SLMCorrection=(angle(SLM_Correction_Pattern))*(255/2/pi);       % correction wavefront in gray values
SLMCorrection=flip(SLMCorrection,1);
figure(); imagesc(SLMCorrection);                       

%% save the important files
dirname = 'P:\TNW\BMPI\Projects\WAVEFRONTSHAPING\data\TPM\3rd gen\191122_WFScomparison_vs_depth_PDMSdiffuser\';
filename = ['d',num2str(d_nom,'%.3d'),'um_model.mat'];
save([dirname,filename],'SLM_Correction_Pattern','SLMCorrection','dnom');