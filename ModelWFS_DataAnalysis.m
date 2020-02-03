%% Load saved stitched files (TPM frames with and without correction)
clear all
close all

load('P:\TNW\BMPI\Users\Abhilash Thendiyammal\Research@UT\Data\191223_WFScomparison_vs_depth_PDMSdiffuser\TPM3D_StichedFiles.mat')

%% plot parameters
Ns=10;                             % size of square considered when determining average intensity
Ithresh = 0;                     % intensity threshold used for bead segmentation
Nframes=150;                       % Number of frames selected for Max Intensity projection            

%% Parameters for converting TPM frames to correct dimensions in um
zoom = 30;                                                 % zoom factor from TPM scan image aquisition  
numPixels= 256;                                            % Pixels in TPM frame
resX = 512/(numPixels*zoom);                               % base x resolution of the image
resY = 512/(numPixels*zoom);                               % base y resolution
x_data= -(size(TPM3Dref,1)/2)*resX:resX:(size(TPM3Dref,1)/2)*resX;
y_data= -(size(TPM3Dref,2)/2)*resY:resY:(size(TPM3Dref,2)/2)*resY;

%% Depth parameters during TPM scanning(Objective piezo stage)
z=1:1:533;                                                                  % total z axis pixels
resZ=0.5;                                                                   % z-resolution
FStart=1;                                                                   % choose the frame to start with
dnom=40:20.5:306.5;                                                         % This distance is from the mean position of the PDMS surface. 
z_startframe=FStart*resZ;                                                   % z of starting frame in um
z_data=dnom(1)+(z_startframe:resZ:size(z,2)*resZ);                          % Relevant z_data for plotting and imaging
n_pdms=1.41;                                                                % PDMS refractive index
n_water=1.33;                                                               % Water refractive index
z_data=z_data*n_pdms/n_water;                                               % Original depth inside PDMS

%% Calculate Maximum Intensity projection
MaxIntensity_TPM3Dref=max(TPM3Dref(end/2-Nframes:end/2+Nframes-1,:,:),[],1);
MaxIntensity_TPM3Dfeedback=max(TPM3Dfeedback(end/2-Nframes:end/2+Nframes-1,:,:),[],1);
MaxIntensity_TPM3Dmodel=max(TPM3Dmodel(end/2-Nframes:end/2+Nframes-1,:,:),[],1);

%% Calculate the average bead intensity from each frame slice
%preallocation
Intensity_ref = zeros(size(TPM3Dref,3),1);
Intensity_feedback = zeros(size(TPM3Dfeedback,3),1);
Intensity_model = zeros(size(TPM3Dmodel,3),1);

% create moving average filter for detecting highest signal in image
f = 1/(Ns^2).*ones(Ns,Ns);

for fn=1:size(TPM3Dref,3)
% consider intensity distribution of single slice
frame_ref = TPM3Dref(:,:,fn);
frame_feedback = TPM3Dfeedback(:,:,fn);
frame_model = TPM3Dmodel(:,:,fn);

% apply moving average filter to git rid of noise
F_ref = conv2(frame_ref,f,'same');
F_feedback = conv2(frame_feedback,f,'same');
F_model = conv2(frame_model,f,'same');

% take maximum signal in low-pass filtered image
Intensity_ref(fn) = max(F_ref(:));
Intensity_feedback(fn) = max(F_feedback(:));
Intensity_model(fn) = max(F_model(:));
end

%% Average Noise level
% Noise=squeeze(TPM3Dref(round(end/2),round(end/2),1));
% Noise=mean(mean((TPM3Dref(1:20,1:20,1))));                                  % mean of the intensity across a small noisy region of the first frame
Noise = mean2(TPM3Dref(:,:,end));
%% Plot the 2D cross-section of intensities (Maximum Intensity projection)
A=squeeze(MaxIntensity_TPM3Dref(1,:,FStart:end));
B=squeeze(MaxIntensity_TPM3Dfeedback(1,:,FStart:end));
C=squeeze(MaxIntensity_TPM3Dmodel(1,:,FStart:end));

Imax_2 =max([A(:);B(:);C(:)]);
figure(1);colormap(hot); 
um = sprintf('(\x0B5m)');
subplot(1,3,1); imagesc(x_data,z_data, A',[0,Imax_2]);  ylabel(['Depth ' um]);  xlabel(['x ' um]); title('a');colorbar;axis on; set(gca,'FontSize',16);
subplot(1,3,2); imagesc(x_data,z_data,B',[0,Imax_2]);  ylabel(['z ' um]);  xlabel(['x ' um]); title('b');colorbar;axis off;set(gca,'FontSize',16);
subplot(1,3,3); imagesc(x_data,z_data,C',[0,Imax_2]); ylabel(['z ' um]);  xlabel(['x ' um]); title('c');colorbar ;axis off; set(gca,'FontSize',16);


%% Semilog plots of Intensity before and after correction
figure(2);
semilogy(z_data,Intensity_ref(:),'bd','MarkerSize',10,'MarkerEdgeColor','b','LineWidth',1.5)
hold on, semilogy(z_data,Intensity_feedback(:),'gs','MarkerSize',10,'MarkerEdgeColor','g','LineWidth',1.5)
hold on,semilogy(z_data,Intensity_model(:),'ro','MarkerSize',10,'MarkerEdgeColor','r','LineWidth',1.5)
D1=ones(1,size(z_data,2)).*double(Noise);
hold on,semilogy(z_data,D1(:),'black:','MarkerSize',20,'MarkerEdgeColor','black','LineWidth',3)
legend('No correction','Feedback based WFS','Model based WFS','average noise level');set(gca,'FontSize',20);ylabel(['Intensity (counts)']);  xlabel(['Depth ' um]);
set(gca,'box','on');set(legend,'box','off')
xlim([dnom(1)+round(FStart/2) 325]);ylim([3e1 1e4]);