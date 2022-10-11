%{
reconstruct the fluence rate map from the slimmed file, then compare to the original output file

Benjamin Kao
Last update: 2020/12/02
%}

clc;clear;close all;

%% param
subject_name_arr={'ZJ','WW','YF','YH','WH','KB','SJ','BT','SC'};
model_dir='models_test'; % the folder containing the voxel model of the subjects
fluence_dir='sim_2E8_literature_sCone1'; % the simulation result should be in fluence_dir / subject_name / fluence_subDir
fluence_subDir='litOP_1';
num_wl=2; % the number of wavelength in a folder

%% main
for sbj=1:length(subject_name_arr)
    %% load the voxel model
    fprintf('Processing %s\n',subject_name_arr{sbj});
    model=load(fullfile(model_dir,['headModel' subject_name_arr{sbj} '_EEG.mat']));
    
    for wl=1:num_wl
        %% load the slimmed fluence rate and reconstruct
        compressed_flu=load(fullfile(fluence_dir,subject_name_arr{sbj},fluence_subDir,['compressed_fluence_' num2str(wl) '.mat']));
        recon_flu=S3_2_fun_reCon_compressedFluence(compressed_flu);

        %% load the original fluence rate file
        orig_flu=load(fullfile(fluence_dir,subject_name_arr{sbj},fluence_subDir,['average_fluence_' num2str(wl) '.mat']));
        orig_flu_mask=ones(size(model.vol));
        orig_flu_mask(compressed_flu.to_save_voxel_index)=0;
        orig_flu_inLayer=orig_flu.average_fluence_rate;
        orig_flu_inLayer(orig_flu_mask==1)=0;

        fprintf('\tThere are %d voxels different from original simulation result.\n',length(find(orig_flu_inLayer~=recon_flu)));

        %% plot the recon and original
        figure('Units','pixels','position',[0 0 1920 1080]);
        ti=tiledlayout('flow','TileSpacing','compact','Padding','compact');
        nexttile()
        imagesc(log10(recon_flu(:,:,90)));
        title('re-constructed');
        ccaxis=caxis();
        colormap jet;
        colorbar;
        nexttile()
        imagesc(log10(orig_flu.average_fluence_rate(:,:,90)));
        title('original');
        caxis(ccaxis);
        colormap jet;
        colorbar;
        drawnow;
    end
end

disp('Done!');