%%% script used to compare feedback-based wavefront shaping with
%%% model-based wavefront shaping.
%%% First run feedback_optimization to obtain ideal_wavefront and model-based WFS scripts
%%% to obtain SLMCorrection

%% General information
f_depth = 300;                       % focus depth (in um) % Note: update in between measurements
data_dirname = '';                   % add data directory here  
wf_dirname = 'wavefronts/';          % directory with ideal wavefront data 
dirname = [data_dirname,wf_dirname]; % combined directory name

%% set flat wavefront
slm.setRect(1,[0,0,0,0]);
slm.setData(0); slm.update; 

%% Set wavefront obtained by feedback-based WFS
load([dirname,'d',num2str(f_depth,'%.3d'),'um_feedback.mat']);
sopt.N = slm.setBlockGeometry(1, sopt.diameter, sopt.N_diameter, sopt.cx, sopt.cy);
slm.setData(ideal_wavefront); slm.update;

%% Set wavefront obtained by model-based WFS
load([dirname,'d',num2str(f_depth,'%.3d'),'um_model.mat']);
slm.setQuadGeometry(1, [-0.4712 0.4886; 0.4698 0.4886; 0.4712 -0.4886; -0.4698 -0.4886]+[sopt.cx sopt.cy]);
slm.setData(SLMCorrection); slm.update;