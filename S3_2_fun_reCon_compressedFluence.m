%{
reconstruct the fluence rate map from the slimmed file

slim_flu: the loaded slimmed fluence rate file structure

recon_flu: the reconstructed fluence rate voxel map

Benjamin Kao
Last update: 2020/12/02
%}

function recon_flu=S3_2_fun_reCon_compressedFluence(slim_flu)

recon_flu=zeros(slim_flu.orig_vol_size);
recon_flu(slim_flu.to_save_voxel_index)=slim_flu.voxel_flu_arr;
    
end