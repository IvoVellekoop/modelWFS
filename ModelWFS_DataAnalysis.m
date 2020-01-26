%% Load saved stitched files (TPM frames with and without correction)
clear all
close all

load('P:\TNW\BMPI\Users\Abhilash Thendiyammal\Research@UT\Data\191223_WFScomparison_vs_depth_PDMSdiffuser\TPM3D_StichedFiles.mat')

%% Parameters for converting TPM frames to correct dimensions in um
zoom = 30;                                                 % zoom factor from TPM scan image aquisition  
numPixels= 256;                                            % Pixels in TPM frame
resX = 512/(numPixels*zoom)                                % base x resolution of the image
resY = 512/(numPixels*zoom)                                % base y resolution
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

Nframes=150;                                                                % Number of frames selected for Max Intensity projection            
MaxIntensity_TPM3Dref=max(TPM3Dref(end/2-Nframes:end/2+Nframes-1,:,:),[],1);
MaxIntensity_TPM3Dfeedback=max(TPM3Dfeedback(end/2-Nframes:end/2+Nframes-1,:,:),[],1);
MaxIntensity_TPM3Dmodel=max(TPM3Dmodel(end/2-Nframes:end/2+Nframes-1,:,:),[],1);

%% Calculate the average intensity counts from each frames around a window size of 11pixels X 11pixels

frame_number=1:1:size(TPM3Dref,3);                      % total frames

for fn=1:numel(frame_number)
    
frame_ref = TPM3Dref(:,:,frame_number(fn));
frame_feedback = TPM3Dfeedback(:,:,frame_number(fn));
frame_model = TPM3Dmodel(:,:,frame_number(fn));

% plot the frames (test)
% Imax = max([frame_ref(:);frame_feedback(:);frame_model(:)]);
% figure();colormap(hot); 
% subplot(1,3,1); imagesc(x_data,y_data,frame_ref,[0,Imax]); xlabel('x (um)');  ylabel('y (um)'); colorbar; axis image; title('TPMimage'); set(gca,'FontSize',16);
% subplot(1,3,2); imagesc(x_data,y_data,frame_feedback,[0,Imax]);  xlabel('x (um)');  ylabel('y (um)'); colorbar; axis image; title('TPMimage with feedbackWFS'); set(gca,'FontSize',16);
% subplot(1,3,3); imagesc(x_data,y_data,frame_model,[0,Imax]);  xlabel('x (um)');  ylabel('y (um)');  colorbar;axis image; title('TPMimage with  modelWFS'); set(gca,'FontSize',16);


% Choose a square region of the image.
window=150;
frame_ref_window=frame_ref(round(end/2)-window:round(end/2)+window-1,round(end/2)-window:round(end/2)+window-1);
frame_feedback_window=frame_feedback(round(end/2)-window:round(end/2)+window-1,round(end/2)-window:round(end/2)+window-1);
frame_model_window=frame_model(round(end/2)-window:round(end/2)+window-1,round(end/2)-window:round(end/2)+window-1);

% Find the maximum intensity at each frame. Add lower and Upper bounds for the data. These bounds are to make sure
% that the data is not close to the boundary where we can select a small
% square region for averaging intensity

lower=20; upper=350;
% reference data
[y_focus_ref_window,x_focus_ref_window] = find(frame_ref_window == max(frame_ref_window(:)),1);
[y_focus_ref,x_focus_ref] = find(frame_ref == max(frame_ref_window(:)),1);
y_focus_ref=y_focus_ref(y_focus_ref>lower);
y_focus_ref=y_focus_ref(y_focus_ref<upper);
x_focus_ref=x_focus_ref(x_focus_ref>lower);
x_focus_ref=x_focus_ref(x_focus_ref<upper);

%repeat for feedback-based
[y_focus_feedback_window,x_focus_feedback_window] = find(frame_feedback_window == max(frame_feedback_window(:)),1);
[y_focus_feedback,x_focus_feedback] = find(frame_feedback == max(frame_feedback_window(:)),1); 
% y_focus_feedback=y_focus_feedback(y_focus_feedback<upper);
% y_focus_feedback=y_focus_feedback(y_focus_feedback>lower);
% x_focus_feedback=x_focus_feedback(x_focus_feedback<upper);
% x_focus_feedback=x_focus_feedback(x_focus_feedback>lower);

%repeat for model-based
[y_focus_model_window,x_focus_model_window] = find(frame_model_window == max(frame_model_window(:)),1);
[y_focus_model,x_focus_model] = find(frame_model== max(frame_model_window(:)),1);
% y_focus_model=y_focus_model(y_focus_model<upper);
% y_focus_model=y_focus_model(y_focus_model>lower);
% x_focus_model=x_focus_model(x_focus_model<upper);
% x_focus_model=x_focus_model(x_focus_model>lower);

% Zoom the data around a maximum intensity point for each frames
Ns=5;                            %  Number of pixels to select square
frame_ref_zoom = frame_ref(y_focus_ref-Ns:y_focus_ref+Ns,x_focus_ref-Ns:x_focus_ref+Ns);
frame_feedback_zoom=frame_feedback(y_focus_feedback-Ns:y_focus_feedback+Ns,x_focus_feedback-Ns:x_focus_feedback+Ns);
frame_model_zoom = frame_model(y_focus_model-Ns:y_focus_model+Ns,x_focus_model-Ns:x_focus_model+Ns);

% figure, imagesc(frame_ref_zoom); ylabel('y (um)'); xlabel('z (um)');
% figure, imagesc(frame_feedback_zoom); ylabel('y (um)'); xlabel('z (um)');
% figure, imagesc(frame_model_zoom); ylabel('y (um)'); xlabel('z (um)');


Intensity_ref(1,fn)=mean(mean(frame_ref_zoom));
Intensity_feedback(1,fn)=mean(mean(frame_feedback_zoom));
Intensity_model(1,fn)=mean(mean(frame_model_zoom));
end

%% Average Noise level
% Noise=squeeze(TPM3Dref(round(end/2),round(end/2),1));
Noise=mean(mean((TPM3Dref(1:20,1:20,1))));                                  % mean of the intensity across a small noisy region of the first frame

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
xlim([dnom(1)+round(FStart/2) 325]);ylim([5e1 9e3]);