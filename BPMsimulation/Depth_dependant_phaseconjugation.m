%% add relevant paths 

addpath('\\ad.utwente.nl\TNW\BMPI\Users\Abhilash Thendiyammal\Research@UT\Matlab programs\git\bpm');
addpath('\\ad.utwente.nl\TNW\BMPI\Users\Abhilash Thendiyammal\Research@UT\Matlab programs\git\wavesim');
addpath('\\ad.utwente.nl\TNW\BMPI\Users\Abhilash Thendiyammal\Research@UT\Matlab programs\git\AO');
addpath('D:\git\tpm\setup');
 
%% setup SLM 
setup;                                  
slm.setQuadGeometry(2, [-0.4712 0.4886; 0.4698 0.4886; 0.4712 -0.4886; -0.4698 -0.4886]+[sopt.cx-0.038 sopt.cy+0.01]);

%% depth dependant phase conjugation
for depth_i=100
    
    % move the sample stage up for 500 um
    Reduced_FieldSize=round(depth_i*tand(37)*2);                    % Field size (Distance between focus and interface*1.33(WATER)/1.51(GP))+Distance between SLM and interface)*tand(37)).
    Reduced_FieldSize= Reduced_FieldSize-mod(Reduced_FieldSize,2);  % Round the value to nearest even number
    NAradius=Reduced_FieldSize                                      % for a unity NA (please note that units are different)
    
    %% Correction wavefront
    CorrectionWF=Sample(end/2-NAradius:end/2+NAradius,end/2-NAradius:end/2+NAradius)./Reference(end/2-NAradius:end/2+NAradius,end/2-NAradius:end/2+NAradius);
    x=linspace(-1,1,2*NAradius+1);
    y=linspace(-1,1,2*NAradius+1)';
    masked_C=mask(x,y, true);
    CorrectionWF_SLM=masked_C(:,:,1).*CorrectionWF;                 % correction wavefront masked with unity NA
    SLM_Correction_Pattern=exp(-1.0i*angle(CorrectionWF_SLM));     
    SLM_angle=angle(SLM_Correction_Pattern);
    figure(); imagesc(SLM_angle);

%% Mapping to SLM

    SLMCorrection=(angle(SLM_Correction_Pattern))*(255/2/pi);       % correction wavefront in gray values
    figure(); imagesc(SLMCorrection);
    figure(); imagesc(flip(SLMCorrection,1));                       

    slm.setQuadGeometry(2, [-0.4712 0.4886; 0.4698 0.4886; 0.4712 -0.4886; -0.4698 -0.4886]+[sopt.cx-0.038 sopt.cy+0.01]); slm.update
    slm.setData(2,flip(SLMCorrection,1)); slm.update;
%     frame_ideal = double(hSI.hDisplay.lastFrame{1});
    slm.setData(2,0); slm.update;
%     frame_ref = double(hSI.hDisplay.lastFrame{1});
    
end

% TPM_3D  = frame_ideal;
% TPM_3D1 = frame_ref;

%% Find source of feedback
