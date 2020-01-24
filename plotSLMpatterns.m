
addpath('../../hardware/matlab');
%% General settings
% data location and figure saving location
dirname = 'P:\TNW\BMPI\Projects\WAVEFRONTSHAPING\data\TPM\3rd gen\191223_WFScomparison_vs_depth_PDMSdiffuser\';
fig_folder = [dirname,'patterns\'];

% displayed data
f_depth_set = 80:20:320; % all data

%% plot settings
diameter = 0.95;    % diameter of the wavefront (0.95 originally, but can be adjusted for more dark space)

% colormap used to represent phase of the images
load('../../utilities/Colormaps/cyclic_colormap2.mat');

%% Create SLM object to display wavefronts on
slm = SLM(3); slm.setData(0); slm.update
test_pattern = double( flipud(slm.getPixels'))*2*pi;
Nx = size(test_pattern);

%% open figure and set size
fig = figure(1); clf; 
set(fig,'Position',[984 61 919 876]);

%% generate figures of Hadamard optimization patterns
for f = 1:numel(f_depth_set)
    fig; clf;
    f_depth = f_depth_set(f);

    % load and generate pattern
    load([dirname,'d',num2str(f_depth,'%.3d'),'um_feedback.mat']);
    slm.setBlockGeometry(1, diameter, sopt.N_diameter, 0, 0);
    slm.setData(ideal_wavefront); slm.update;
    
    % retrieve image of pattern and convert phase to rads
    optimized_pattern = double( flipud(slm.getPixels'))*2*pi/256;
    optimized_pattern = optimized_pattern(:,end/2+(-Nx(1)/2+1:Nx(1)/2));
    
    % display figure
    imagesc(optimized_pattern,[0,2*pi]); colormap(cm);
    axis image; axis off;
    
    % save figures
    filename = ['d',num2str(f_depth,'%.3d'),'um_feedback_pattern'];
    saveas(fig,[fig_folder,filename,'.png']);
    saveas(fig,[fig_folder,filename,'.fig']);
    saveas(fig,[fig_folder,filename,'.pdf']);
end

%% generate figures of model WFS optimization patterns
for f = 1:numel(f_depth_set)
    fig; clf;
    f_depth = f_depth_set(f);

    % load and generate pattern
    load([dirname,'d',num2str(f_depth,'%.3d'),'um_model.mat']);
    slm.setQuadGeometry(1, [-0.4712 0.4886; 0.4698 0.4886; 0.4712 -0.4886; -0.4698 -0.4886]);
    slm.setData(SLMCorrection); slm.update;
    
    % retrieve image of pattern and convert phase to rads
    optimized_pattern = double( flipud(slm.getPixels'))*2*pi/256;
    optimized_pattern = optimized_pattern(:,end/2+(-Nx(1)/2+1:Nx(1)/2));
    
    % display figure
    imagesc(optimized_pattern,[0,2*pi]); colormap(cm);
    axis image; axis off;
    
    % save figures
    filename = ['d',num2str(f_depth,'%.3d'),'um_model_pattern'];
    saveas(fig,[fig_folder,filename,'.png']);
    saveas(fig,[fig_folder,filename,'.fig']);
    saveas(fig,[fig_folder,filename,'.pdf']);
end