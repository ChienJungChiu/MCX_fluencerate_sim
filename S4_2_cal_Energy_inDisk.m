%{
Reconstruct the fluence rate map from the slimmed file, and only calculate the fluence rate in GM disk area

Benjamin Kao
Last update: 2021/03/15
%}

clc;clear;close all;

%% param
subject_name_arr={'ZJ','WW','YF','YH','WH','KB','SJ','BT','SC'};
model_dir='models'; % the folder containing the voxel model of the subjects
fluence_dir='sim_2E8_literature_sCone1'; % the simulation result should be in fluence_dir / subject_name / fluence_subDir
fluence_subDir='litOP_1';
num_wl=2; % the number of wavelength in a folder
source_r=22.5; % mm

%% main
sbj_wl_flu_arr=[]; % store the fluence rate sum in GM for each subject and each wl
for sbj=1:length(subject_name_arr)
    %% load the voxel model
    fprintf('Processing %s\n',subject_name_arr{sbj});
    
    %% load the probe pos and dir
    load(fullfile(model_dir,[subject_name_arr{sbj} '_inDiskGM.mat']));
    
    for wl=1:num_wl
        %% load the slimmed fluence rate and reconstruct
        compressed_flu=load(fullfile(fluence_dir,subject_name_arr{sbj},fluence_subDir,['compressed_fluence_' num2str(wl) '.mat']));
        recon_flu=S3_2_fun_reCon_compressedFluence(compressed_flu);
        in_GM_flu=sum(recon_flu(superficial_inRange_noWM_GM>0));
        sbj_wl_flu_arr(sbj,wl)=in_GM_flu;
    end
end

disp('Done!');