%%% Script to convert the TPM 3D  image to a refractive index
%%% distribution. This code fits the the intensity along the each z-columns
%%% to sigmoid function and find out the point of inflection. We assign the
%%% refractive index values at the both sides of the poit of inflection
%%% with the apriori knowledge of refractive indices. Here, 1.33 and 1.41
%%% for flurescene and PDMS meadia respectively.
%%% Edited by Abhilash Thendiyammal 2018

%%
clf
% close all
% clear all
addpath('D:\git\utilities');

%% In order to find the interface, fit each Z column to a sigmoid function
% Stack_3D is acquired from ShearTPMimage program. 

Sub_Stack_3D=Stack_3D(:,:,1:N);                   % Make a substack so that we can test the fitiing with less number of images.
xPoints=1:N;                                      % fit N2 number of data points in each Z column of Sub_Stack_3D
starttime = now;

for l=1:size(Sub_Stack_3D,1)                      % loop through X-coordinates of each Image
    for m=1:size(Sub_Stack_3D,2)                  % loop through Y-coordinates of each image
        a=Sub_Stack_3D(l,m,:);                    % select eack z-columns
        yPoints=a(:)';
        m1=min(min(a));                           % parameter for fitting
        m2=max(max(a));                           % parameter for fitting
        [maxvalue,indx]=max(max(a));
        maxIndex{l,m}=[l m indx];
        maxIndexZ(l, m)=maxIndex{l, m}(3);
        stde=std(a(maxIndexZ(l,m):N));        
        fit=createFit(xPoints,yPoints,m1,m2,stde);% Fit a sigmoid function to eack z-column.
        c=coeffvalues(fit); 
        fn=c(1)+((c(2)-c(1))./(1+exp(-c(3).*xPoints(:))*c(4))); % Fitted Sigmoid function
        dif=diff(fn);
        [maxv,idx]=max(dif);                      % Find the indices corresponding to the peak in first derivative of sigmoid function. 
        out{l,m}=[l m idx];                       % Indices corresponding to the point of inflection.
    end
    
    eta(l, size(Sub_Stack_3D,1), starttime, 'console', 'Fitting stuff...', 0);
end

%% Converting cell format to a array format for position of point of inflection
tic
Size_vec=size(out);
maximum=zeros(Size_vec);
for row=1:Size_vec(1)
    for column=1:Size_vec(2)
        maximum(row, column)=out{row, column}(3);  
    end
end
figure,surf(maximum);
PDMS_thickness=round(median(median(maximum)));
toc
%% Interpolate the Zmax corresponding to the point of inflection. 
zmultiplier=3;
IntX=linspace(1,size(Sub_Stack_3D,1),5*size(Sub_Stack_3D,1));IntY=linspace(1,size(Sub_Stack_3D,1),5*size(Sub_Stack_3D,1)); IntZ=linspace(1,60,60*zmultiplier);
[IntX,IntY]=ndgrid(IntX,IntY);

maximum_filtered = medfilt2(maximum);

% % Isolate pixels with the median filtered image
% maximum_mask = ones(size(maximum));
% maximum_mask(maximum > maximum_filtered*2) = 0;
% maximum_corrected = ~maximum_mask .* maximum_filtered + maximum_mask .* maximum;

MAX = round(zmultiplier*interpn(maximum_filtered,IntX,IntY, 'cubic'));

Super_Stack_3D=zeros(5*size(Sub_Stack_3D,1),5*size(Sub_Stack_3D,1),60*zmultiplier);
nInt=5;
Size_vec2=size(out)*nInt;
MAX(MAX<0)=[1];
for row2=1:Size_vec2(1)
    for column2=1:Size_vec2(2)
        Super_Stack_3D(row2,column2,1:(MAX(row2, column2)))=1.41;
        Super_Stack_3D(row2,column2,(MAX(row2, column2))+1:N*zmultiplier)=1.33;  
    end
end
        
%% Plot the refractive index and the TPM intensity

for kk=20:22
figure, imagesc(squeeze(Stack_3D(:,kk,:)))
figure, imagesc(squeeze(Super_Stack_3D(:,5*kk,:)))
end
Super_Stack_3D=Super_Stack_3D(:,:,1:180);

clearvars -except Stack_3D Super_Stack_3D MAX hSICtl hSI PDMS_thickness