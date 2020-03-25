%% Script for finding conversion matrices. The image 
%  acquired from TPM images are in pixels %  This has to be converted into
%  micrometers for using the data in other codes like ShearTPMimage.m,
%  ModelWFS_DataStitching etc.

clc
clear all 
close all

%% Add relevant paths
addpath('C:\git\tpm\setup');
addpath('C:\git\PCT\Experiments');

%% setup SLM and experimetal parameters
active_devices.pmt_gain=true;
active_devices.sample_stage=true;
setup

%% Define Conversion matrices 
M=zeros(2); 
G=zeros(2);

%% Apply shift to sample with Zaber stage and acquire images
d =50;                                                       % displacement in um
zoom=10;
z2.moveDistance(d);
frame_refX = double(hSI.hDisplay.lastFrame{1});             % Record a reference frame

z2.moveDistance(d);                                         % Move sample in X direction 
frame_SampleXshifted = double(hSI.hDisplay.lastFrame{1});   % Record a shifted frame

z1.moveDistance(d);                                         % Move sample back to initial position (-X direction)
frame_refY = double(hSI.hDisplay.lastFrame{1});             % Record a reference frame

z1.moveDistance(d);                                        % Move sample in Y direction 
frame_SampleYshifted = double(hSI.hDisplay.lastFrame{1});   % Record a shifted frame

%% Calculate the correlation between the images to find X and Y shift

CCX1= xcorr2(frame_SampleXshifted,frame_refX);               % Cross-correlation of shifted frame with reference
[max1, imax1] = max(abs(CCX1(:)));
[ypeak1, xpeak1] = ind2sub(size(CCX1),imax1);
corr_offset1 = [(xpeak1-size(frame_refX,1)) (size(frame_refX,2))-ypeak1 ];

CCY1= xcorr2(frame_SampleYshifted,frame_refY);               % Cross-correlation of shifted frame with reference
[max2, imax2] = max(abs(CCY1(:)));
[ypeak2, xpeak2] = ind2sub(size(CCY1),imax2(1));
corr_offset2 = [(xpeak2-size(frame_refY,1)) (size(frame_refY,2))-ypeak2];

M(:,1)= corr_offset1;
M(:,2)= corr_offset2;

%% gradient parameters
bg_patch_id = 1;                                            % Background Gradient Patch ID
bg_pix = 1152;                                              % Number of SLM pixels (set it as slm width- set for -0.5 to 0.5)
bgtype = 'blaze';                                           % Type of grating ('sine', 'blaze', 'square')
pmin = 0;                                                   % Minimum phase
pmax = 255;                                                 % Maximum phase
ppp = 16;                                                   % pixels per period(2pi)

bg = bg_grating(bgtype, ppp, pmin, pmax, bg_pix);

%% Apply gradient to SLM and acquire TPM frames

slm.setRect(bg_patch_id, [sopt.cx sopt.cy 1 1]);             %keep the aspect ratio always 1:1
slm.setData(bg_patch_id, 0); slm.update;

frame_refX2 = double(hSI.hDisplay.lastFrame{1});             % Record a reference frame

slm.setData(bg_patch_id, bg); slm.update; 
frame_Xshifted = double(hSI.hDisplay.lastFrame{1});          % Record X shifted frame

slm.setData(bg_patch_id, 0); slm.update; 
frame_refY2 = double(hSI.hDisplay.lastFrame{1});             % Record a reference frame

slm.setData(bg_patch_id, bg'); slm.update; 
frame_Yshifted = double(hSI.hDisplay.lastFrame{1});          % Record Y shifted frame

%% Calculate the correlation between the images to find X and Y shift

CCX2= xcorr2(frame_Xshifted,frame_refX2);                     % Cross-sorrelation of shifted frame with reference
[max1, imax1] = max(abs(CCX2(:)));
[ypeak1, xpeak1] = ind2sub(size(CCX2),imax1(1));
corr_offset1 = [(xpeak1-size(frame_refX2,1)) (size(frame_refX2,2))-ypeak1];

CCY2= xcorr2(frame_Yshifted,frame_refY2);                     % Cross-sorrelation of shifted frame with reference
[max2, imax2] = max(abs(CCY2(:)));
[ypeak2, xpeak2] = ind2sub(size(CCY2),imax2(1));
corr_offset2 = [(xpeak2-size(frame_refY2,1)) (size(frame_refY2,2))-ypeak2];

G(:,1)= corr_offset2;
G(:,2)= corr_offset1;

%% Conversion matrix from gradient to angle.

mx=49.7; my=-50.06;                                           % Equivalent to 50um shift in X and Y direction respectively
m=[mx my]; 
M0=M./m;                                                      % Conversion matrix for TPM image calibration

G0=G.*8/pi;                                                   % Conversion matrix for SLM pixels to image shift

%% Calculate the coordinates for mapping the Simulated (kx,ky) to SLM pixels. Four vertices of a square 
% transforms to four vertices of a quadrilateral structure. 

A=[6.284 6.284];                                   % 2*pi (equivalent to kNA value of unity)
Q1 = (1/1152).* G0*(M0)^-1*A';
B=[6.284 -6.284];                                  % 2*pi (equivalent to kNA value of unity)
Q2 = (1/1152).* G0*(M0)^-1*B';
C=[-6.284 -6.284];                                 % 2*pi (equivalent to kNA value of unity)
Q3 = (1/1152).* G0*(M0)^-1*C';
D=[-6.284 6.284];                                  % 2*pi (equivalent to kNA value of unity)
Q4 = (1/1152).* G0*(M0)^-1*D';







