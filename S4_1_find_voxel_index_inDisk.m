%{
Find the energy in the GM covered by disk area
Constrain the GM in superficial area.

Benjamin Kao
Last update: 2021/03/15
%}

clc;clear;close all;

%% param
subject_name_arr={'ZJ','WW','YF','YH','WH','KB','SJ','BT','SC'};
model_dir='models'; % the folder containing the voxel model of the subjects
source_r=22.5; % mm
max_GM_depth=10; % mm, the max depth in to calculate in GM

%% main
for sbj=1:length(subject_name_arr)
    %% load the voxel model
    fprintf('Processing %s\n',subject_name_arr{sbj});
    model=load(fullfile(model_dir,['headModel' subject_name_arr{sbj} '_EEG.mat']));
    voxel_size=model.voxel_size;
    
    %% load the probe pos and dir
    p_pos=load(fullfile(model_dir,[subject_name_arr{sbj} '_disk1_probe_pos.txt']));
    p_dir=load(fullfile(model_dir,[subject_name_arr{sbj} '_disk1_probe_dir.txt']));
    source_pos=p_pos(1,:);
    source_dir=p_dir(1,:);
    
    orig_GM_voxel=model.vol==4;

    %% find the voxel in the disk area
    [xx,yy,zz]=ndgrid(1:size(orig_GM_voxel,1),1:size(orig_GM_voxel,2),1:size(orig_GM_voxel,3));

    dist2source_square=(xx-source_pos(1)).^2+(yy-source_pos(2)).^2+(zz-source_pos(3)).^2;

    point_pos_arr=[reshape(xx,[],1) reshape(yy,[],1) reshape(zz,[],1)];

    % the vector from source to point
    source2point=point_pos_arr-source_pos;

    % dot of (source to point) and (source direction)
    dot_SP_SD=sum(source2point.*source_dir,2);

    % find the reference point on the axis that closest to the points
    on_axis_ref_arr=source_pos+source_dir.*dot_SP_SD;

    % calculate the distance between point to axis
    dist2axis=sqrt(sum((point_pos_arr-on_axis_ref_arr).^2,2));
    dist2axis=reshape(dist2axis,size(orig_GM_voxel,1),size(orig_GM_voxel,2),size(orig_GM_voxel,3));

    source_voxel_index=dist2axis<source_r/voxel_size;
    GM_inDisk_index=source_voxel_index & orig_GM_voxel;
    
    %% find the superficial GM and not cover by WM
    
    % Find the boundary between GM and CSF
    boundary_points=find_GM_CSF_boundary(model.vol);
    
    % the boundary point in the disk
    boundary_point_inDisk=boundary_points & GM_inDisk_index;
    
    % make a mask, which is the area cover by the first layer of GM
    boundary_mask=make_shift_mask_image(boundary_point_inDisk,source_dir,200);
    
    % shift the boundary points by 1 voxel
    shifted_boundary_mask=shift_mask(boundary_mask,source_dir);
    
    
    % find the first layer of GM
    first_layer_GM=boundary_point_inDisk;
    first_layer_GM(shifted_boundary_mask>0)=0;
    
    % find the superficial GM voxel
    superficial_GM=make_shift_mask_image(first_layer_GM,source_dir,max_GM_depth/voxel_size);
    superficial_GM=superficial_GM & orig_GM_voxel;
    
    % find the area cover by WM
    WM_mask=make_shift_mask_image(model.vol==5,source_dir,200);
    
    % find the GM not cover by WM
    superficial_noWM_GM=superficial_GM;
    superficial_noWM_GM(WM_mask>0)=0;
    
    %% make sure not select GM too deep
    % find the distance of source to the closest selected GM
    selected_points=point_pos_arr(superficial_noWM_GM>0,:);
    dist_selectedPoints_to_source=sqrt(sum((selected_points-source_pos).^2,2));
    min_dist_to_source=min(dist_selectedPoints_to_source);
    
    % find the area close to source
    in_range_index=sqrt(dist2source_square)<=min_dist_to_source+max_GM_depth+source_r/voxel_size;
    
    superficial_inRange_noWM_GM=superficial_noWM_GM;
    superficial_inRange_noWM_GM(in_range_index==0)=0;
    
    
%     a=superficial_inRange_noWM_GM+orig_GM_voxel;
%     sliceViewer(a);

    save(fullfile(model_dir,[subject_name_arr{sbj} '_inDiskGM.mat']),'superficial_inRange_noWM_GM');    
end

disp('Done!');

%% functions

% find the boundary between GM and CSF
% boundary_points is GM tissue.
function boundary_points=find_GM_CSF_boundary(head_model)
boundary_index_arr=zeros(size(head_model));

temp_slice=circshift(head_model,1,1);
boundary_index_arr=boundary_index_arr | (head_model-temp_slice)==1;
temp_slice=circshift(head_model,-1,1);
boundary_index_arr=boundary_index_arr | (head_model-temp_slice)==1;
temp_slice=circshift(head_model,1,2);
boundary_index_arr=boundary_index_arr | (head_model-temp_slice)==1;
temp_slice=circshift(head_model,-1,2);
boundary_index_arr=boundary_index_arr | (head_model-temp_slice)==1;
temp_slice=circshift(head_model,1,3);
boundary_index_arr=boundary_index_arr | (head_model-temp_slice)==1;
temp_slice=circshift(head_model,-1,3);
boundary_index_arr=boundary_index_arr | (head_model-temp_slice)==1;
boundary_index_arr=boundary_index_arr & head_model==4;
boundary_points=boundary_index_arr;
end

% find the area covered by a moving mask
% orig_mask: original moving mask, a 0/1 image
% dir: moving direction
% len: moving length
function masked_image=make_shift_mask_image(orig_mask,dir,len)
del_l=0.5;

len=ceil(len);

% make a dilate kernel
kernel_image=zeros(2*len+1,2*len+1,2*len+1);
center_point=[len+1,len+1,len+1];
for l=0:del_l:len
    kernel_image(ceil(center_point(1)+l*dir(1)),ceil(center_point(2)+l*dir(2)),ceil(center_point(3)+l*dir(3)))=1;
    kernel_image(round(center_point(1)+l*dir(1)),round(center_point(2)+l*dir(2)),round(center_point(3)+l*dir(3)))=1;
    kernel_image(floor(center_point(1)+l*dir(1)),floor(center_point(2)+l*dir(2)),floor(center_point(3)+l*dir(3)))=1;
end

% do the dilation to have a masked image
se = strel(kernel_image);
masked_image=imdilate(orig_mask,se);

end

% find the new mask, which is a shifted original mask
% orig_mask: original mask, a 0/1 image
% dir: shift direction
function shifted_mask=shift_mask(orig_mask,dir)
if dir(1)>0
    shifted_mask=orig_mask([1 1:end-1],:,:);
    shifted_mask(1,:,:)=0;
else
    shifted_mask=orig_mask([2:end end],:,:);
    shifted_mask(end,:,:)=0;
end
if dir(2)>0
    shifted_mask=shifted_mask(:,[1 1:end-1],:);
    shifted_mask(:,1,:)=0;
else
    shifted_mask=shifted_mask(:,[2:end end],:);
    shifted_mask(:,end,:)=0;
end
if dir(3)>0
    shifted_mask=shifted_mask(:,:,[1 1:end-1]);
    shifted_mask(:,:,1)=0;
else
    shifted_mask=shifted_mask(:,:,[2:end end]);
    shifted_mask(:,:,end)=0;
end
end