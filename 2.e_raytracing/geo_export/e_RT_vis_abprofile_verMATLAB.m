% eLeaf: 3D model of dicot leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 0.1

%% visualize light absorption in 3D using MATLAB
%% also can use COMSOL to plot the light profile
%% - light absorption of each chloroplast comes ray tracing; 
%%      -- after use trace_recal to summerize all files from paralleled tasks
%%      -- in fout_srf, last two columns are curobj->I_absorb_o and curobj->I_absorb_i. For chloroplasts, curobj->I_absorb_i is enough 
%% - volume of chloroplasts from e_geo_AVR_cal_v0_1.m; 
%%      -- VAR_NAME=chlo_vol_profile
%%      -- MAT_NAME=chlo_vol_profile.mat or SAVE_e_geom4RT.mat

tmp_RTresults_path='/home/xiaoyi/eLeaf_dicot/eleaf_dicot_v0_1_PAL_r2c2/4.e_RT_recal_comparison/directRT_blue445nm_rep1/';
tmp_RTresults_file_name_rtsum='results_merged_rtsum_445nm_500x_rep1';
tmp_RTresults_file_name_absrf='results_merged_absrf_445nm_500x_rep1';

tmp_RTresults_rtsum=importdata([tmp_RTresults_path,tmp_RTresults_file_name_rtsum]);
tmp_RTresults_absrf=importdata([tmp_RTresults_path,tmp_RTresults_file_name_absrf])/tmp_RTresults_rtsum.data(5);

tmp_idx=find(tmp_RTresults_absrf(:,3)~=0);
tmp_ab_chlo=tmp_RTresults_absrf(tmp_idx,7);

load SAVE_e_geom4RT.mat chlo_vol_profile

tmp_ab_chlo_pervol=tmp_ab_chlo./chlo_vol_profile;
tmp_ab_chlo_pervol_log=log(tmp_ab_chlo_pervol);


%% plot 3D leaf and colormap chloroplasts light absorption
alpha=0.4%0.25;
linecolor='none';%[0.5,0.5,0.5];%'none';%[0.8,0.8,0.8];
color_wall=[0.3,0.3,0.3];
%color_chl=[0.2,0.8,0.2];
%color_vac=[0.2,0.2,0.8];

cmap=colormap('hot');%% return a 256*3 vector
cval4cmapmin=min(tmp_ab_chlo_pervol_log);
cval4cmapmax=max(tmp_ab_chlo_pervol_log);

load('SAVE_RT_geo_export.mat', 'GLB_tag_*')

count_chlo=0;
%%%% 2. upper epidermis
for loop_i=1:numel(GLB_tag_EPL_u_set)
    tmp_name=['../nonMS/EPL_u_',num2str(loop_i),'.ply'];
    [tri,pts]=ply_read(tmp_name,'tri');
    trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'Facecolor',color_wall,'facealpha',alpha,'EdgeColor',linecolor);hold on;
end

%%%% 3-5. PAL
for loop_i=1:4%numel(GLB_tag_PAL_MS)
    %%%% 3. PAL MS
    tmp_name=['../MS/PAL_MS_',num2str(loop_i),'.ply'];
    [tri,pts]=ply_read(tmp_name,'tri');
    %trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'Facecolor',color_wall,'facealpha',alpha,'EdgeColor',linecolor);hold on;
    
    %%%% 4. PAL VAC
    tmp_name=['../MS/PAL_VAC_',num2str(loop_i),'.ply'];
    [tri,pts]=ply_read(tmp_name,'tri');
    %trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'Facecolor',color_vac,'facealpha',alpha,'EdgeColor',linecolor);hold on;
    
    %%%% 5. PAL CHL
    tmptag_set_PAL_CHL=GLB_tag_PAL_CHL{loop_i};
    for loop_i_CHL=1:numel(tmptag_set_PAL_CHL)
        tmp_name=['../MS/PAL_MS_',num2str(loop_i),'_CHL_',num2str(loop_i_CHL),'.ply'];
        [tri,pts]=ply_read(tmp_name,'tri');
        count_chlo=count_chlo+1;
        tmp_color_idx=round((tmp_ab_chlo_pervol_log(count_chlo)-cval4cmapmin)/(cval4cmapmax-cval4cmapmin)*255+1);
        trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'Facecolor',cmap(tmp_color_idx,:),'facealpha',alpha,'EdgeColor',linecolor);hold on;
    end 
end

%%%% 6-8. SPO
for loop_i=1:numel(GLB_tag_SPO_MS)
    %%%% 6. SPO MS
    tmp_name=['../MS/SPO_MS_',num2str(loop_i),'.ply'];
    [tri,pts]=ply_read(tmp_name,'tri');
    %trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'Facecolor',color_wall,'facealpha',alpha,'EdgeColor',linecolor);hold on;
    
    %%%% 7. SPO VAC
    tmp_name=['../MS/SPO_VAC_',num2str(loop_i),'.ply'];
    [tri,pts]=ply_read(tmp_name,'tri');
    %trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'Facecolor',color_vac,'facealpha',alpha,'EdgeColor',linecolor);hold on;
    
    %%%% 8. SPO CHL
    tmptag_set_SPO_CHL=GLB_tag_SPO_CHL{loop_i};
    for loop_i_CHL=1:numel(tmptag_set_SPO_CHL)
        tmp_name=['../MS/SPO_MS_',num2str(loop_i),'_CHL_',num2str(loop_i_CHL),'.ply'];
        [tri,pts]=ply_read(tmp_name,'tri');
        count_chlo=count_chlo+1;
        tmp_color_idx=round((tmp_ab_chlo_pervol_log(count_chlo)-cval4cmapmin)/(cval4cmapmax-cval4cmapmin)*255+1);
        trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'Facecolor',cmap(tmp_color_idx,:),'facealpha',alpha,'EdgeColor',linecolor);hold on;
    end
end

%%%% 1. lower epidermis
for loop_i=1:numel(GLB_tag_EPL_l_set)
    tmp_name=['../nonMS/EPL_l_',num2str(loop_i),'.ply'];
    [tri,pts]=ply_read(tmp_name,'tri');
    trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'Facecolor',color_wall,'facealpha',alpha,'EdgeColor',linecolor);hold on;
end

axis equal;

hbar=colorbar;
%% map hbar.Limits to [cval4cmapmin, cval4cmapmax];
hbar_tick_cval=29:34;
hbar_tick_labels={'29','30','31','32','33','34'};
hbar_tick_true=(hbar_tick_cval-cval4cmapmin)./(cval4cmapmax-cval4cmapmin)*(hbar.Limits(2)-hbar.Limits(1))+hbar.Limits(1);
set(hbar,'Ticks',hbar_tick_true);
set(hbar,'TickLabels',hbar_tick_labels);
