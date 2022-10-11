%{
Find the position and direction of the probe before simulation.

Benjamin Kao
Last update: 2020/12/02
%}

clc;clear;close all;

model_folder='models_test'; % the folder containing the model
subject_name_arr={'ZJ','WW','YF','YH','WH','KB','SJ','BT','SC'}; % the name of the subjects

find_cone_mode=2; % =1 for using fix source distance to head; =2 for using radius and angle to find the source location
if find_cone_mode==1
    % for using fix source distance to head to find the source location, the location of source is only decided by `shift_distance`, and the `equivalent_r_arr` and `setting_angle_arr` are only for plot figure
    shift_distance=[21 11 17.8]; % in mm, the distance between source and head
    equivalent_r_arr=[22.5 20 15]; % in mm
    setting_angle_arr=[30 45 60]; % in degree
elseif find_cone_mode==2
    % for using radius and angle to find the source location, the location of source is find by adjust the distance of source to head to let the cone ( with certain angle )shine on the head is like a circle ( with certain radius )
    equivalent_r_arr=[22.5 22.5 22.5 20 20 20 15 15 15 10 10 10]; % in mm, the radius of the circle in which the light will shine on
    setting_angle_arr=[30 45 60 30 45 60 30 45 60 30 45 60]; % in degree, the angle of the cone want to compare with the circle
end

to_plot_figure=1; % plot the probe position and direction or not

%% init
if find_cone_mode==1
    assert(length(shift_distance)==length(equivalent_r_arr),'Error: source setting not match');
    assert(length(shift_distance)==length(setting_angle_arr),'Error: source setting not match');
    assert(length(equivalent_r_arr)==length(setting_angle_arr),'Error: source setting not match');
elseif find_cone_mode==2
    assert(length(equivalent_r_arr)==length(setting_angle_arr),'Error: source setting not match');
end

%% main
for sbj=1:length(subject_name_arr)
    subject_name=subject_name_arr{sbj};
    MRI_voxel_file=fullfile(model_folder,['headModel' subject_name '_EEG.mat']); % containing the MRI voxel model, also the EEG point and head surface mesh
    
    MRI_model=load(MRI_voxel_file);
    if isfield(MRI_model,'model_version')==0
        MRI_model.model_version=1;
    end
    
    %% find the position and direction of probes
    [orig_p_pos,p_dir]=fun_cal_probe_position_and_direction_2D(MRI_model.EEG,[0],[0],MRI_model.vol,MRI_model.headsurf,MRI_model.voxel_size,[2],MRI_model.model_version);
    
    save(fullfile(model_folder,[subject_name '_orig_probe_pos.txt']),'orig_p_pos','-ascii','-tabs');
    
    %% change location for cone source
    % using the setting angle and setting radius
    
    for di=1:length(equivalent_r_arr)
        p_pos=orig_p_pos;
        if find_cone_mode==1
            p_pos(1,:)=p_pos(1,:)-p_dir(1,:)*shift_distance(di)/MRI_model.voxel_size;
        elseif find_cone_mode==2
            changed_p_pos=fitting_for_cone_source(MRI_model.headsurf,p_pos(1,:),p_dir(1,:),equivalent_r_arr(di),setting_angle_arr(di),MRI_model.voxel_size);
            p_pos(1,:)=changed_p_pos;
        end
        save(fullfile(model_folder,[subject_name '_cone' num2str(di) '_probe_pos.txt']),'p_pos','-ascii','-tabs');
        save(fullfile(model_folder,[subject_name '_cone' num2str(di) '_probe_dir.txt']),'p_dir','-ascii','-tabs');
        
        %% plot the probes
        if to_plot_figure
            plot_probe_position(MRI_model.vol,MRI_model.headsurf,0,p_pos,p_dir,p_pos(1,:),p_dir(1,:),equivalent_r_arr(di),setting_angle_arr(di),model_folder,subject_name,MRI_model.model_version);
        end
        close all;
    end
end

disp('Done!');

%% functions

%% plot the 2D and 3D plot for the probe position to make sure the probe is on the subject's head surface
function plot_probe_position(vol,headsurf,num_SDS,p_pos,p_dir,src_pos,src_dir,set_r,set_angle,model_folder,subject_name,model_version)
%% param
do_plot_2D=0;
do_plot_3D=1;
do_plot_disk_source=1;
do_plot_cone_source=1;
cone_hight=70;

dir_len=10;

%% main

image_ratio=size(vol,1)/size(vol,2);
if image_ratio>1
    image_ratio=1/image_ratio;
end


if do_plot_2D
    for s=1:num_SDS+1
        figure('Units','pixels','position',[0 0 round(image_ratio*800) 800]);
        imagesc(vol(:,:,round(p_pos(s,3))));
        hold on;
        if model_version==1
            plot(p_pos(s,1),p_pos(s,2),'xr','LineWidth',2,'MarkerSize',24);
            quiver(p_pos(s,1)-p_dir(s,1)*dir_len,p_pos(s,2)-p_dir(s,2)*dir_len,p_dir(s,1),p_dir(s,2),dir_len,'k','LineWidth',3,'MarkerSize',24);
        elseif model_version>=2
            plot(p_pos(s,2),p_pos(s,1),'xr','LineWidth',2,'MarkerSize',24); % swap the x and y
            quiver(p_pos(s,2)-p_dir(s,2)*dir_len,p_pos(s,1)-p_dir(s,1)*dir_len,p_dir(s,2),p_dir(s,1),dir_len,'k','LineWidth',3,'MarkerSize',24); % swap the x and y
        end
        if s==1
            title(['source at layer ' num2str(round(p_pos(s,3)))]);
        else
            title(['SDS ' num2str(s-1) ' at layer ' num2str(round(p_pos(s,3)))]);
        end
        saveas(gcf,fullfile(model_folder,[subject_name '_SDS_' num2str(s-1) '.png']));
        close all;
    end
end

if do_plot_3D
    %% plot the head and the SDSs
    figure('Units','pixels','position',[0 0 1400 1080]);
    patch(headsurf,'FaceColor',[1,.75,.65],'EdgeColor','none','FaceAlpha','0.9');
    lighting gouraud
    lightangle(0,30);
    lightangle(120,30);
    lightangle(-120,30);
    axis('image');
    hold on;
    plot3(p_pos(:,1),p_pos(:,2),p_pos(:,3),'ro','MarkerSize',10,'LineWidth',2);
    for s=1:size(p_pos,1)
        if s==1
            quiver3(p_pos(s,1)-p_dir(s,1)*dir_len,p_pos(s,2)-p_dir(s,2)*dir_len,p_pos(s,3)-p_dir(s,3)*dir_len,p_dir(s,1)*dir_len,p_dir(s,2)*dir_len,p_dir(s,3)*dir_len,'k','LineWidth',3,'MarkerSize',24);
            text(p_pos(s,1)-p_dir(s,1)*dir_len,p_pos(s,2)-p_dir(s,2)*dir_len,p_pos(s,3)-p_dir(s,3)*dir_len,'source','FontSize',20,'Color','g');
        else
            quiver3(p_pos(s,1),p_pos(s,2),p_pos(s,3),p_dir(s,1)*dir_len,p_dir(s,2)*dir_len,p_dir(s,3)*dir_len,'k','LineWidth',3,'MarkerSize',24);
            text(p_pos(s,1)+p_dir(s,1)*dir_len,p_pos(s,2)+p_dir(s,2)*dir_len,p_pos(s,3)+p_dir(s,3)*dir_len,['SDS ' num2str(s-1)],'FontSize',20,'Color','g');
        end
    end
    
    if do_plot_disk_source
        %% generate cylinder for disk source
        [xx,yy,zz]=cylinder(set_r,100);
        zz=zz.*cone_hight;
        
        % find the angle to turn
        orig_cyl_dir=[0 0 1];
        phy=acos(src_dir(3)); % altitude angle
        orig_cyl_dir=[cos(phy) 0 sin(phy);0 1 0;-sin(phy) 0 cos(phy)]*orig_cyl_dir'; % rotate on x-z plan
        theta=acos(src_dir(1)./sqrt(sum(src_dir(1:2).^2))); % azimuth angle
        if model_version>=2
            theta=-theta;
        end
        orig_cyl_dir=[cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1]*orig_cyl_dir;  % rotate on x-y plan
        orig_cyl_dir=orig_cyl_dir'; % this should be the same as p_dir
        
        temp_cyl_arr=[reshape(xx,[],1) reshape(yy,[],1) reshape(zz,[],1)];
        temp_cyl_arr=[cos(phy) 0 sin(phy);0 1 0;-sin(phy) 0 cos(phy)]*temp_cyl_arr';
        temp_cyl_arr=[cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1]*temp_cyl_arr;
        temp_cyl_arr=temp_cyl_arr';
        xx=reshape(temp_cyl_arr(:,1),size(xx))+src_pos(1);
        yy=reshape(temp_cyl_arr(:,2),size(yy))+src_pos(2);
        zz=reshape(temp_cyl_arr(:,3),size(zz))+src_pos(3);
        
        surf(xx,yy,zz);
    end
    
    if do_plot_cone_source
        [xx,yy,zz]=cylinder([0 cone_hight*tan(set_angle/180*pi)],100);
        zz=zz*cone_hight;
        
        % find the angle to turn
        orig_cyl_dir=[0 0 1];
        phy=acos(src_dir(3)); % altitude angle
        orig_cyl_dir=[cos(phy) 0 sin(phy);0 1 0;-sin(phy) 0 cos(phy)]*orig_cyl_dir'; % rotate on x-z plan
        theta=acos(src_dir(1)./sqrt(sum(src_dir(1:2).^2))); % azimuth angle
        if model_version>=2
            theta=-theta;
        end
        orig_cyl_dir=[cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1]*orig_cyl_dir;  % rotate on x-y plan
        orig_cyl_dir=orig_cyl_dir'; % this should be the same as p_dir
        
        temp_cyl_arr=[reshape(xx,[],1) reshape(yy,[],1) reshape(zz,[],1)];
        temp_cyl_arr=[cos(phy) 0 sin(phy);0 1 0;-sin(phy) 0 cos(phy)]*temp_cyl_arr';
        temp_cyl_arr=[cos(theta) -sin(theta) 0; sin(theta) cos(theta) 0; 0 0 1]*temp_cyl_arr;
        temp_cyl_arr=temp_cyl_arr';
        xx=reshape(temp_cyl_arr(:,1),size(xx))+src_pos(1);
        yy=reshape(temp_cyl_arr(:,2),size(yy))+src_pos(2);
        zz=reshape(temp_cyl_arr(:,3),size(zz))+src_pos(3);
        
        surf(xx,yy,zz);
    end
    
    view(-src_dir);
    saveas(gcf,fullfile(model_folder,['headModel' subject_name '_fiber_plot.fig']));
end
end

% the old function to find the cone source position
function changed_p_pos=find_cone_source_to_shift(vol,orig_pos,orig_dir,set_r,set_angle,voxel_size)
%% param
r_multiply=1.5; % if the distance between found point and source is wierd, try to adjust this value.

%% main
set_r=set_r/voxel_size; % change the set_r from mm to voxel
%% find the layer 1 boundary
layer_1_boundary=zeros(size(vol));
for i=1:3
    for j=[-1 1]
        layer_1_boundary=layer_1_boundary | (vol-circshift(vol,j,i))==1;
    end
end
layer_1_boundary=layer_1_boundary & vol==1;
is_boundary=find(layer_1_boundary==1);

%% find the distance of the points to the source, choose the close points
[xx,yy,zz]=ndgrid(1:size(vol,1),1:size(vol,2),1:size(vol,3));
xx=xx(is_boundary);
yy=yy(is_boundary);
zz=zz(is_boundary);

vol_mid_point=size(vol)/2;
mid_to_points=[xx yy zz]-vol_mid_point;
mid_to_source=orig_pos([2 1 3])-vol_mid_point; % swap the x and y
is_half=find(sum(mid_to_points.*mid_to_source,2)>=0);

%% find the distance of the points to the source axis
xx=xx(is_half);
yy=yy(is_half);
zz=zz(is_half);

S2P=[xx yy zz]-orig_pos([2 1 3]); % swap the x and y
distSquare_2_source=sum(S2P.^2,2);
dot_SP_sDir=sum(S2P.*orig_dir([2 1 3]),2); % the dot of (source to point) and source direction
distSquare_2_axis=distSquare_2_source-dot_SP_sDir.^2;

%% find the point closest to the set r
temp=abs(distSquare_2_axis-set_r^2);
[~,min_index]=min(temp(:));
fprintf('The distance between find point (close to set_r) and source is %f voxel, and %f voxel to axis\n',sqrt(sum(([xx(min_index) yy(min_index) zz(min_index)]-orig_pos([2 1 3])).^2)),sqrt(distSquare_2_axis(min_index)));

%% calculate the pos that source should be
source_offset=set_r*cot(set_angle/180*pi);
changed_p_pos=orig_pos-(source_offset-dot_SP_sDir(min_index))*orig_dir;
fprintf('Move the position for cone source, %f voxel for plan, %f voxel for this head model\n',source_offset,sqrt(sum((changed_p_pos-orig_pos).^2)));

end

%{
Fitting to find the position of light source that the shined area of the cone source and the disk source are almost the same.

Inputs:
headsurf: the head surface mesh structure
orig_pos: the original source position on the head surface
orig_dir: the direction of the source
set_r: the radius of the shined area
set_angle: the angle of cone source
voxel_size: the side length of each voxel

set_angle: the angle of cone source
boundary_half_points: the coordinate of the points that is on the half head close to the source
is_shined_index: the index of points shined by disk source.

Outputs:
changed_p_pos: The new position the source should be.
%}
function changed_p_pos=fitting_for_cone_source(headsurf,orig_pos,orig_dir,set_r,set_angle,voxel_size)
%% main
set_r=set_r/voxel_size; % change the set_r from mm to voxel

%% cut the half head far from the source
headsurf=headsurf.vertices;

vol_mid_point=mean(headsurf,1);
mid_to_points=headsurf-vol_mid_point;
mid_to_source=orig_pos-vol_mid_point;
is_half=find(sum(mid_to_points.*mid_to_source,2)>=0);

%% find the distance of the points to the source axis
headsurf=headsurf(is_half,:);
S2P=headsurf-orig_pos;
dot_SP_sDir=sum(S2P.*orig_dir,2); % the dot of (source to point) and source direction
distSquare_2_source=sum(S2P.^2,2);
distSquare_2_axis=distSquare_2_source-dot_SP_sDir.^2;

%% find the points will be shined by disk source
is_shined=find(distSquare_2_axis<=set_r^2);

%% fitting to find the cone source
source_offset=set_r*cot(set_angle/180*pi); % the offset if is on plane
x0=source_offset-1;
x=fminsearch(@(x)find_dist_to_shift(x,orig_pos,orig_dir,set_angle,headsurf,is_shined),x0); % swap the x and y
changed_p_pos=orig_pos-x*orig_dir;
fprintf('Move the position for cone source, %f voxel for plan, %f voxel for this head model\n',source_offset,sqrt(sum((changed_p_pos-orig_pos).^2)));
end

%{
Find the shined voxel by cone source, and compare them to those shined by disk source

Inputs:
x: how long the cone source should move
orig_pos: the original source position on the head surface
orig_dir: the direction of the source
set_angle: the angle of cone source
boundary_half_points: the coordinate of the points that is on the half head close to the source
is_shined_index: the index of points shined by disk source.

Outputs:
error: The unsimilarity of the voxel shined by 2 sources
%}
function error=find_dist_to_shift(x,orig_pos,orig_dir,set_angle,boundary_half_points,is_shined_index)
new_pos=orig_pos-x*orig_dir;
newPos_to_points=boundary_half_points-new_pos;
dist_to_points=sqrt(sum(newPos_to_points.^2,2));
newPos_to_points_norm=newPos_to_points./dist_to_points; % normalize the length
cosine_NPP_dir=sum(newPos_to_points_norm.*orig_dir,2); % the cosine of (newPos to points) and the ( source direction ), because both the length are 1, the cosine = dot.
cone_shined_index=find(cosine_NPP_dir>=cos(set_angle/180*pi) & sqrt(sum((boundary_half_points-orig_pos).^2,2))<1.5*22.5); % find the voxel index in the setting angle
error=1-(find_array_similarity(is_shined_index,cone_shined_index));
% fprintf('x=%d, error= %10f\n',x,error);
end

% Find the simularity of 2 array
% if the 2 input array are identical, the output will be 1
function sim=find_array_similarity(a,b)
a=sort(a);
b=sort(b);
a_index=1;
b_index=1;
equal_count=0;
while a_index<=length(a) && b_index<=length(b)
    if a(a_index)==b(b_index)
        equal_count=equal_count+1;
        a_index=a_index+1;
        b_index=b_index+1;
    elseif a(a_index)<b(b_index)
        a_index=a_index+1;
    elseif a(a_index)>b(b_index)
        b_index=b_index+1;
    end
end
sim=equal_count*2/(length(a)+length(b));
end