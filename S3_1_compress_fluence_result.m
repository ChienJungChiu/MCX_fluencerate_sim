%{
Only store the fluence rate for the GM/WM layer, and larger than certain threshold, in order to compress the simulation result

Benjamin Kao
Last update: 2020/12/02
%}

clc;clear;close all;

%% param
subject_name_arr={'ZJ','WW','YF','YH','WH','KB','SJ','BT','SC'}; % the name of the subjects
model_dir='models_test'; % the folder containing the voxel model of the subjects
fluence_dir='sim_2E8_literature_sCone1'; % the simulation result should be in fluence_dir / subject_name / fluence_subDir
fluence_subDir='litOP_1';
to_save_layer=[4 5]; % the index for layers to save
num_wl=2; % the number of wavelength in a folder
energy_threshold=1E-11; % only store the voxel who's fluence rate/energy larger than this threshold.
num_layers=6; % the number of layer, not include type 0
do_delete_orig_file=0; % =1 to delete the original fluence rate file

%% main
for sbj=1:length(subject_name_arr)
    %% load the voxel model
    fprintf('Processing %s\n',subject_name_arr{sbj});
    model=load(fullfile(model_dir,['headModel' subject_name_arr{sbj} '_EEG.mat']));
    orig_vol_size=size(model.vol); % the size of the original voxel model
    sbj_to_save_voxel_index=zeros(size(model.vol)); % which voxel is the to_save_layer in all voxels
    for i=1:length(to_save_layer)
        sbj_to_save_voxel_index=sbj_to_save_voxel_index | model.vol==to_save_layer(i);
    end
    sbj_to_save_voxel_index=find(sbj_to_save_voxel_index>0);
    fprintf('\tThe number of voxel in layers to save / all voxel = %.2f%%\n',100*length(sbj_to_save_voxel_index)/length(model.vol(:)));
    
    for wl=1:num_wl
        fprintf('\tWavelength %d\n',wl);
        %% load the fluence rate simulation result
        sim_flu=load(fullfile(fluence_dir,subject_name_arr{sbj},fluence_subDir,['average_fluence_' num2str(wl) '.mat']));
        assert(sum(size(model.vol)~=size(sim_flu.average_fluence_rate))==0,'Error! The size of the model and the fluence rate are different!');

        %% calculate all energy in these layers
        layer_total_flu=zeros(1,num_layers);
        for L=1:num_layers
            layer_total_flu(L)=sum(sim_flu.average_fluence_rate(model.vol==L));
        end

        %% find the voxel have larger energy than threshold
        to_save_layer_flu=sim_flu.average_fluence_rate(sbj_to_save_voxel_index);
        larger_than_threshold_index=find(to_save_layer_flu>energy_threshold); % which voxel in the to_save_layer is larger than threshold
        fprintf('\tThe number of voxel larger than threshold / all voxel in layers to save = %.2f%%\n',100*length(larger_than_threshold_index)/length(sbj_to_save_voxel_index));
        fprintf('\tThe number of voxel larger than threshold / all voxel = %.2f%%\n',100*length(larger_than_threshold_index)/length(model.vol(:)));

        to_save_voxel_index=int32(sbj_to_save_voxel_index(larger_than_threshold_index)); % which voxel is larger than threshold in all voxels
        voxel_flu_arr=sim_flu.average_fluence_rate(to_save_voxel_index);

        %% calculate the total energy loss
        orig_total_flu=sum(layer_total_flu(to_save_layer));
        new_total_flu=sum(voxel_flu_arr);
        fprintf('\tTotal energy in layers to save drop from %.2e to %.2e, that is %.2f%%\n',orig_total_flu,new_total_flu,(new_total_flu/orig_total_flu-1)*100);

        %% save
        save(fullfile(fluence_dir,subject_name_arr{sbj},fluence_subDir,['compressed_fluence_' num2str(wl) '.mat']),'to_save_layer','energy_threshold','layer_total_flu','to_save_voxel_index','voxel_flu_arr','orig_vol_size')

        orig_file_info=dir(fullfile(fluence_dir,subject_name_arr{sbj},fluence_subDir,['average_fluence_' num2str(wl) '.mat']));
        new_file_info=dir(fullfile(fluence_dir,subject_name_arr{sbj},fluence_subDir,['compressed_fluence_' num2str(wl) '.mat']));
        fprintf('\tSlim file from %.2f MB to %.2f MB, that is %.2f%% of the original size.\n',orig_file_info.bytes/1E6,new_file_info.bytes/1E6,new_file_info.bytes/orig_file_info.bytes*100);

        if do_delete_orig_file==1 && exist(fullfile(fluence_dir,subject_name_arr{sbj},fluence_subDir,['compressed_fluence_' num2str(wl) '.mat']),'file')
            delete(fullfile(fluence_dir,subject_name_arr{sbj},fluence_subDir,['average_fluence_' num2str(wl) '.mat']));
        end
    end
end

disp('Done!');