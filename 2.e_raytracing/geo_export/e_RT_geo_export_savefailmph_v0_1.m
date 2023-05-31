% eLeaf: 3D model of dicot leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 0.1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update 2020-Jan
% excute after e_RT_geo_export_vis_geo
% - which will delete the rest geometry
% - save mph for manually check

clearvars

load ../../1.e_geom/TMP_CFG.mat

import com.comsol.model.*
import com.comsol.model.util.*
ModelUtil.remove('Model_3Dleaf_dicot4RT')
%mphtags -show
load('TMP_matlab.mat')
load('SAVE_RT_geo_export.mat','GLB_threshold','GLB_max_mesh_length_*');

if isnan(CFG_VARARGIN(1))==1
    %% input mode for this script
    prompt = 'Which MAT file to load the SAVE_FAIL_LIST_ext_*?:\n1=SAVE_RT_geo_export_aftervisualcheck_full.mat, \n2=TMP_SAVE_RT_geo_export_aftervisualcheck_full.mat, \n3=SAVE_RT_geo_export.mat, \n';
    GLB_savefailmph_mode = input(prompt);
else
    GLB_savefailmph_mode = CFG_VARARGIN(1);
end

if GLB_savefailmph_mode==1
    load('SAVE_RT_geo_export_aftervisualcheck_full.mat', 'SAVE*');%SAVE_FAIL_list_ext_???
elseif GLB_savefailmph_mode==2
    %% supplementary savefailmph model
    %% after visualization full geometry or ray tracing, find new obj with mesh error
    %% SUGGEST write nodeleterepair information also to this mat file.
    load('TMP_SAVE_RT_geo_export_aftervisualcheck_full.mat', 'SAVE*');%SAVE_FAIL_list_ext_???
else %GLB_savefailmph_mode==3
    % skipped manual check
    load('SAVE_RT_geo_export.mat', 'SAVE*');
    SAVE_FAIL_LIST_ext_EPL_u=SAVE_FAIL_LIST_EPL_u;
    SAVE_FAIL_LIST_ext_EPL_l=SAVE_FAIL_LIST_EPL_l;
    SAVE_FAIL_LIST_ext_PAL_MS=SAVE_FAIL_LIST_PAL_MS;
    SAVE_FAIL_LIST_ext_PAL_VAC=SAVE_FAIL_LIST_PAL_VAC;
    SAVE_FAIL_LIST_ext_PAL_CHL=SAVE_FAIL_LIST_PAL_CHL;
    SAVE_FAIL_LIST_ext_PAL_VACCHL=[];
    SAVE_FAIL_LIST_ext_SPO_MS=SAVE_FAIL_LIST_SPO_MS;
    SAVE_FAIL_LIST_ext_SPO_VACCHL=[];
    SAVE_FAIL_LIST_ext_SPO_VAC=SAVE_FAIL_LIST_SPO_VAC;%only used in 2nd batch
    SAVE_FAIL_LIST_ext_SPO_CHL=SAVE_FAIL_LIST_SPO_CHL;%only used in 2nd batch
end
if exist('SAVE_FAIL_LIST_ext_EPL_l','var')==0
    error('[ErrorMsg]: No input list for saving failed mph files.')
end

if isnan(CFG_VARARGIN(2))==1
    %% select include delete repair or not; used later
    prompt = 'Include Repair/Vitual Operation for saved failed mph files?:\n1=Yes (Default, suggested), 2=No\n';
    GLB_nodeleterepair_mode = input(prompt);
else
    GLB_nodeleterepair_mode = CFG_VARARGIN(2);
end

model=mphopen('TMP_pre_RT_geo_export.mph','Model_3Dleaf_dicot4RT');

count_csur=1;
%% 2. delete the resting objs; convert domain to surface; and create selection
tmptag_del=['del',num2str(count_del)];count_del=count_del+1;%% create del and csur operations
model.component('comp1').geom('geom1').feature.create(tmptag_del, 'Delete');
tmptag_csur=['csur', num2str(count_csur)];count_csur=count_csur+1;
model.component('comp1').geom('geom1').create(tmptag_csur, 'ConvertToSurface');

mphsave(model, 'TMP2_pre_RT_geo_export.mph')
GLB_count_fail=0;

% GLB_threshold=1e-8; %defined in e_RT_geo_export_v0_1.m
% SAVE_FAIL_LIST_EPL_l=[];%No.1
% SAVE_FAIL_LIST_EPL_u=[];
% SAVE_FAIL_LIST_PAL_MS=[];
% SAVE_FAIL_LIST_PAL_VAC=[];
% SAVE_FAIL_LIST_PAL_CHL=[];
% SAVE_FAIL_LIST_SPO_MS=[];
% SAVE_FAIL_LIST_SPO_VAC=[];
% SAVE_FAIL_LIST_SPO_CHL=[];%No.8
% SAVE_COUNT_FAIL=zeros(1,8);
SAVE_COUNT_FAIL_r2=zeros(1,8);

%%%% 2.1 lower epidermis
TMP_count_fail=0;
for TMP_loop_i=1:numel(SAVE_FAIL_LIST_ext_EPL_l)
    loop_i=SAVE_FAIL_LIST_ext_EPL_l(TMP_loop_i);
    GLB_count_fail=GLB_count_fail+1;
    
    if GLB_count_fail~=1
        tic
        ModelUtil.remove('Model_3Dleaf_dicot4RT');
        model=mphopen('TMP2_pre_RT_geo_export.mph','Model_3Dleaf_dicot4RT');
        toc
    end
    
    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').init;
    clear del_list
    del_list=[GLB_tag_EPL_l_set([1:loop_i-1,loop_i+1:numel(GLB_tag_EPL_l_set)]), ...
        GLB_tag_EPL_u_set, ...
        GLB_tag_PAL_MS, ...
        GLB_tag_PAL_VAC, ...
        cat(2,GLB_tag_PAL_CHL{:}), ...
        GLB_tag_SPO_MS, ...
        GLB_tag_SPO_VAC, ...
        cat(2,GLB_tag_SPO_CHL{:})];
    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(del_list);
    
    model.component('comp1').geom('geom1').feature('csur1').selection('input').set(GLB_tag_EPL_l_set(loop_i));
    
    if GLB_nodeleterepair_mode==1
        %%%%%% delete small faces and edges
        model.component('comp1').geom('geom1').defeaturing('SmallFaces').set('entsize', GLB_threshold);
        model.component('comp1').geom('geom1').defeaturing('SmallFaces').find;
        if (model.component('comp1').geom('geom1').defeaturing('SmallFaces').detail().size())~=0
            try
                model.component('comp1').geom('geom1').defeaturing('SmallFaces').deleteAll('dsf1');
            catch
                disp('[ErrorMsg]: cannot delete small faces');
            end
        end
        model.component('comp1').geom('geom1').defeaturing('ShortEdges').set('entsize', GLB_threshold);
        model.component('comp1').geom('geom1').defeaturing('ShortEdges').find;
        if (model.component('comp1').geom('geom1').defeaturing('ShortEdges').detail().size())~=0
            try
                model.component('comp1').geom('geom1').defeaturing('ShortEdges').deleteAll('dse1');
            catch
                disp('[ErrorMsg]: cannot delete short edges');
            end
        end
    end
    
    model.component('comp1').mesh.create('mesh1');
    model.component('comp1').mesh('mesh1').create('size1', 'Size');
    model.component('comp1').mesh('mesh1').feature('size1').set('hauto', 4);
    %model.component('comp1').mesh('mesh1').feature('size1').set('hmax', GLB_max_mesh_length);
    model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
    model.component('comp1').mesh('mesh1').feature('ftri1').selection.all;
%     try
%         model.component('comp1').mesh('mesh1').run('ftri1');
%         
%         [meshstats,meshdata] = mphmeshstats(model,'mesh1');
%         pts=meshdata.vertex;
%         pts=pts';
%         tri=meshdata.elem{2};
%         tri=(tri+1)';
%         %trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facealpha',0.5);
%         %axis equal;hold on
%         clear face_ms;
%         for i=1:size(tri,1)
%             face_ms(i)={[tri(i,1)-1,tri(i,2)-1,tri(i,3)-1]};
%         end
%         data.vertex.x=pts(:,1);
%         data.vertex.y=pts(:,2);
%         data.vertex.z=pts(:,3);
%         data.face.vertex_indices=face_ms;
%         tmp_filename=['EPL_l_',num2str(loop_i),'.ply'];
%         ply_write(data,tmp_filename,'ascii','double');
%     catch
    if GLB_nodeleterepair_mode==1
        tmp_filename=['Debug_EPL_l_',num2str(loop_i),'.mph'];
    else
        tmp_filename=['Debug_nodel_EPL_l_',num2str(loop_i),'.mph'];
    end
        mphsave(model, tmp_filename)
        TMP_count_fail=TMP_count_fail+1;
%     end
end
SAVE_COUNT_FAIL_r2(1)=TMP_count_fail;

%%%% 2.2 upper epidermis
TMP_count_fail=0;
for TMP_loop_i=1:numel(SAVE_FAIL_LIST_ext_EPL_u)
    loop_i=SAVE_FAIL_LIST_ext_EPL_u(TMP_loop_i);
    GLB_count_fail=GLB_count_fail+1;
    
    if GLB_count_fail~=1
        tic
        ModelUtil.remove('Model_3Dleaf_dicot4RT');
        model=mphopen('TMP2_pre_RT_geo_export.mph','Model_3Dleaf_dicot4RT');
        toc
    end
    
    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').init;
    clear del_list
    del_list=[GLB_tag_EPL_l_set, ...
        GLB_tag_EPL_u_set([1:loop_i-1,loop_i+1:numel(GLB_tag_EPL_u_set)]), ...
        GLB_tag_PAL_MS, ...
        GLB_tag_PAL_VAC, ...
        cat(2,GLB_tag_PAL_CHL{:}), ...
        GLB_tag_SPO_MS, ...
        GLB_tag_SPO_VAC, ...
        cat(2,GLB_tag_SPO_CHL{:})];
    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(del_list);
    
    model.component('comp1').geom('geom1').feature('csur1').selection('input').set(GLB_tag_EPL_u_set(loop_i));
    
    if GLB_nodeleterepair_mode==1
        %%%%%% delete small faces and edges
        model.component('comp1').geom('geom1').defeaturing('SmallFaces').set('entsize', GLB_threshold);
        model.component('comp1').geom('geom1').defeaturing('SmallFaces').find;
        if (model.component('comp1').geom('geom1').defeaturing('SmallFaces').detail().size())~=0
            try
                model.component('comp1').geom('geom1').defeaturing('SmallFaces').deleteAll('dsf1');
            catch
                disp('[ErrorMsg]: cannot delete small faces');
            end
        end
        model.component('comp1').geom('geom1').defeaturing('ShortEdges').set('entsize', GLB_threshold);
        model.component('comp1').geom('geom1').defeaturing('ShortEdges').find;
        if (model.component('comp1').geom('geom1').defeaturing('ShortEdges').detail().size())~=0
            try
                model.component('comp1').geom('geom1').defeaturing('ShortEdges').deleteAll('dse1');
            catch
                disp('[ErrorMsg]: cannot delete short edges');
            end
        end
    end
    
    model.component('comp1').mesh.create('mesh1');
    model.component('comp1').mesh('mesh1').create('size1', 'Size');
    model.component('comp1').mesh('mesh1').feature('size1').set('hauto', 4);
    %model.component('comp1').mesh('mesh1').feature('size1').set('hmax', GLB_max_mesh_length);
    model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
    model.component('comp1').mesh('mesh1').feature('ftri1').selection.all;
%     try
%         model.component('comp1').mesh('mesh1').run('ftri1');
%         
%         [meshstats,meshdata] = mphmeshstats(model,'mesh1');
%         pts=meshdata.vertex;
%         pts=pts';
%         tri=meshdata.elem{2};
%         tri=(tri+1)';
%         %trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facealpha',0.5);
%         %axis equal;hold on
%         clear face_ms;
%         for i=1:size(tri,1)
%             face_ms(i)={[tri(i,1)-1,tri(i,2)-1,tri(i,3)-1]};
%         end
%         data.vertex.x=pts(:,1);
%         data.vertex.y=pts(:,2);
%         data.vertex.z=pts(:,3);
%         data.face.vertex_indices=face_ms;
%         tmp_filename=['EPL_u_',num2str(loop_i),'.ply'];
%         ply_write(data,tmp_filename,'ascii','double');
%     catch
    if GLB_nodeleterepair_mode==1
        tmp_filename=['Debug_EPL_u_',num2str(loop_i),'.mph'];
    else
        tmp_filename=['Debug_nodel_EPL_u_',num2str(loop_i),'.mph'];
    end
        mphsave(model, tmp_filename)
        TMP_count_fail=TMP_count_fail+1;
%     end
end
SAVE_COUNT_FAIL_r2(2)=TMP_count_fail;

%%%% 2.3 PAL cell
TMP_count_fail=0;
for TMP_loop_i=1:numel(SAVE_FAIL_LIST_ext_PAL_MS)
    loop_i=SAVE_FAIL_LIST_ext_PAL_MS(TMP_loop_i);
    GLB_count_fail=GLB_count_fail+1;
    
    if GLB_count_fail~=1
        tic
        ModelUtil.remove('Model_3Dleaf_dicot4RT');
        model=mphopen('TMP2_pre_RT_geo_export.mph','Model_3Dleaf_dicot4RT');
        toc
    end
    
    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').init;
    clear del_list
    del_list=[GLB_tag_EPL_l_set, ...
        GLB_tag_EPL_u_set, ...
        GLB_tag_PAL_MS([1:loop_i-1,loop_i+1:numel(GLB_tag_PAL_MS)]), ...
        GLB_tag_PAL_VAC, ...
        cat(2,GLB_tag_PAL_CHL{:}), ...
        GLB_tag_SPO_MS, ...
        GLB_tag_SPO_VAC, ...
        cat(2,GLB_tag_SPO_CHL{:})];
    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(del_list);
    
    model.component('comp1').geom('geom1').feature('csur1').selection('input').set(GLB_tag_PAL_MS(loop_i));
    
    if GLB_nodeleterepair_mode==1
        %%%%%% delete small faces and edges
        model.component('comp1').geom('geom1').defeaturing('SmallFaces').set('entsize', GLB_threshold);
        model.component('comp1').geom('geom1').defeaturing('SmallFaces').find;
        if (model.component('comp1').geom('geom1').defeaturing('SmallFaces').detail().size())~=0
            try
                model.component('comp1').geom('geom1').defeaturing('SmallFaces').deleteAll('dsf1');
            catch
                disp('[ErrorMsg]: cannot delete small faces');
            end
        end
        model.component('comp1').geom('geom1').defeaturing('ShortEdges').set('entsize', GLB_threshold);
        model.component('comp1').geom('geom1').defeaturing('ShortEdges').find;
        if (model.component('comp1').geom('geom1').defeaturing('ShortEdges').detail().size())~=0
            try
                model.component('comp1').geom('geom1').defeaturing('ShortEdges').deleteAll('dse1');
            catch
                disp('[ErrorMsg]: cannot delete short edges');
            end
        end
    end
    
    model.component('comp1').mesh.create('mesh1');
    model.component('comp1').mesh('mesh1').create('size1', 'Size');
    model.component('comp1').mesh('mesh1').feature('size1').set('hauto', 4);
    model.component('comp1').mesh('mesh1').feature('size1').set('hmax', GLB_max_mesh_length_MSVAC);
    model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
    model.component('comp1').mesh('mesh1').feature('ftri1').selection.all;
%     try
%         model.component('comp1').mesh('mesh1').run('ftri1');
%         
%         [meshstats,meshdata] = mphmeshstats(model,'mesh1');
%         pts=meshdata.vertex;
%         pts=pts';
%         tri=meshdata.elem{2};
%         tri=(tri+1)';
%         %trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facealpha',0.5);
%         %axis equal;hold on
%         clear face_ms;
%         for i=1:size(tri,1)
%             face_ms(i)={[tri(i,1)-1,tri(i,2)-1,tri(i,3)-1]};
%         end
%         data.vertex.x=pts(:,1);
%         data.vertex.y=pts(:,2);
%         data.vertex.z=pts(:,3);
%         data.face.vertex_indices=face_ms;
%         tmp_filename=['PAL_MS_',num2str(loop_i),'.ply'];
%         ply_write(data,tmp_filename,'ascii','double');
%     catch
    if GLB_nodeleterepair_mode==1
        tmp_filename=['Debug_PAL_MS_',num2str(loop_i),'.mph'];
    else
        tmp_filename=['Debug_nodel_PAL_MS_',num2str(loop_i),'.mph'];
    end
        mphsave(model, tmp_filename)
        TMP_count_fail=TMP_count_fail+1;
%     end
end
SAVE_COUNT_FAIL_r2(3)=TMP_count_fail;

%%%% 2.4 PAL VAC
TMP_count_fail=0;
for TMP_loop_i=1:numel(SAVE_FAIL_LIST_ext_PAL_VAC)
    loop_i=SAVE_FAIL_LIST_PAL_ext_VAC(TMP_loop_i);
    GLB_count_fail=GLB_count_fail+1;
    
    if GLB_count_fail~=1
        tic
        ModelUtil.remove('Model_3Dleaf_dicot4RT');
        model=mphopen('TMP2_pre_RT_geo_export.mph','Model_3Dleaf_dicot4RT');
        toc
    end
    
    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').init;
    clear del_list
    del_list=[GLB_tag_EPL_l_set, ...
        GLB_tag_EPL_u_set, ...
        GLB_tag_PAL_MS, ...
        GLB_tag_PAL_VAC([1:loop_i-1,loop_i+1:numel(GLB_tag_PAL_VAC)]), ...
        cat(2,GLB_tag_PAL_CHL{:}), ...
        GLB_tag_SPO_MS, ...
        GLB_tag_SPO_VAC, ...
        cat(2,GLB_tag_SPO_CHL{:})];
    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(del_list);
    
    model.component('comp1').geom('geom1').feature('csur1').selection('input').set(GLB_tag_PAL_VAC(loop_i));
    
    if GLB_nodeleterepair_mode==1
        %%%%%% delete small faces and edges
        model.component('comp1').geom('geom1').defeaturing('SmallFaces').set('entsize', GLB_threshold);
        model.component('comp1').geom('geom1').defeaturing('SmallFaces').find;
        if (model.component('comp1').geom('geom1').defeaturing('SmallFaces').detail().size())~=0
            try
                model.component('comp1').geom('geom1').defeaturing('SmallFaces').deleteAll('dsf1');
            catch
                disp('[ErrorMsg]: cannot delete small faces');
            end
        end
        model.component('comp1').geom('geom1').defeaturing('ShortEdges').set('entsize', GLB_threshold);
        model.component('comp1').geom('geom1').defeaturing('ShortEdges').find;
        if (model.component('comp1').geom('geom1').defeaturing('ShortEdges').detail().size())~=0
            try
                model.component('comp1').geom('geom1').defeaturing('ShortEdges').deleteAll('dse1');
            catch
                disp('[ErrorMsg]: cannot delete short edges');
            end
        end
    end
    
    model.component('comp1').mesh.create('mesh1');
    model.component('comp1').mesh('mesh1').create('size1', 'Size');
    model.component('comp1').mesh('mesh1').feature('size1').set('hauto', 4);
    model.component('comp1').mesh('mesh1').feature('size1').set('hmax', GLB_max_mesh_length_MSVAC);
    model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
    model.component('comp1').mesh('mesh1').feature('ftri1').selection.all;
%     try
%         model.component('comp1').mesh('mesh1').run('ftri1');
%         
%         [meshstats,meshdata] = mphmeshstats(model,'mesh1');
%         pts=meshdata.vertex;
%         pts=pts';
%         tri=meshdata.elem{2};
%         tri=(tri+1)';
%         %trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facealpha',0.5);
%         %axis equal;hold on
%         clear face_ms;
%         for i=1:size(tri,1)
%             face_ms(i)={[tri(i,1)-1,tri(i,2)-1,tri(i,3)-1]};
%         end
%         data.vertex.x=pts(:,1);
%         data.vertex.y=pts(:,2);
%         data.vertex.z=pts(:,3);
%         data.face.vertex_indices=face_ms;
%         tmp_filename=['PAL_VAC_',num2str(loop_i),'.ply'];
%         ply_write(data,tmp_filename,'ascii','double');
%     catch
    if GLB_nodeleterepair_mode==1
        tmp_filename=['Debug_PAL_VAC_',num2str(loop_i),'.mph'];
    else
        tmp_filename=['Debug_nodel_PAL_VAC_',num2str(loop_i),'.mph'];
    end
        mphsave(model, tmp_filename)
        TMP_count_fail=TMP_count_fail+1;
%     end
end
SAVE_COUNT_FAIL_r2(4)=TMP_count_fail;

%%%% 2.5 PAL CHL
TMP_count_fail=0;
for TMP_loop_i=1:size(SAVE_FAIL_LIST_ext_PAL_CHL,1)
    loop_i=SAVE_FAIL_LIST_ext_PAL_CHL(TMP_loop_i,1);
    loop_i_CHL=SAVE_FAIL_LIST_ext_PAL_CHL(TMP_loop_i,2);
    GLB_count_fail=GLB_count_fail+1;
    
    if GLB_count_fail~=1
        tic
        ModelUtil.remove('Model_3Dleaf_dicot4RT');
        model=mphopen('TMP2_pre_RT_geo_export.mph','Model_3Dleaf_dicot4RT');
        toc
    end
    
    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').init;
    clear del_list
    GLB_tag_PAL_CHL_tmpdel=GLB_tag_PAL_CHL;
    tmp_list=GLB_tag_PAL_CHL{loop_i};
    %GLB_tag_PAL_CHL_tmpdel{loop_i}=tmp_list([1:loop_i_CHL-1,loop_i_CHL+1:GLB_count_PAL_chl(loop_i)]);
    GLB_tag_PAL_CHL_tmpdel{loop_i}=tmp_list([1:loop_i_CHL-1,loop_i_CHL+1:numel(tmp_list)]);
    del_list=[GLB_tag_EPL_l_set, ...
        GLB_tag_EPL_u_set, ...
        GLB_tag_PAL_MS, ...
        GLB_tag_PAL_VAC, ...
        cat(2,GLB_tag_PAL_CHL_tmpdel{:}), ...
        GLB_tag_SPO_MS, ...
        GLB_tag_SPO_VAC, ...
        cat(2,GLB_tag_SPO_CHL{:})];
    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(del_list);
    
    model.component('comp1').geom('geom1').feature('csur1').selection('input').set(tmp_list(loop_i_CHL));
    
    if GLB_nodeleterepair_mode==1
        %%%%%% delete small faces and edges
        model.component('comp1').geom('geom1').defeaturing('SmallFaces').set('entsize', GLB_threshold);
        model.component('comp1').geom('geom1').defeaturing('SmallFaces').find;
        if (model.component('comp1').geom('geom1').defeaturing('SmallFaces').detail().size())~=0
            try
                model.component('comp1').geom('geom1').defeaturing('SmallFaces').deleteAll('dsf1');
            catch
                disp('[ErrorMsg]: cannot delete small faces');
            end
        end
        model.component('comp1').geom('geom1').defeaturing('ShortEdges').set('entsize', GLB_threshold);
        model.component('comp1').geom('geom1').defeaturing('ShortEdges').find;
        if (model.component('comp1').geom('geom1').defeaturing('ShortEdges').detail().size())~=0
            try
                model.component('comp1').geom('geom1').defeaturing('ShortEdges').deleteAll('dse1');
            catch
                disp('[ErrorMsg]: cannot delete short edges');
            end
        end
    end
    
    model.component('comp1').mesh.create('mesh1');
    model.component('comp1').mesh('mesh1').create('size1', 'Size');
    model.component('comp1').mesh('mesh1').feature('size1').set('hauto', 4);
    model.component('comp1').mesh('mesh1').feature('size1').set('hmax', GLB_max_mesh_length_CHL);
    model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
    model.component('comp1').mesh('mesh1').feature('ftri1').selection.all;
    model.component('comp1').mesh('mesh1').feature('ftri1').set('zscale', 0.5);%not tested
    %%if zscale works, the following try-catch would save lots of time
    try
        %model.component('comp1').mesh('mesh1').feature('ftri1').set('xscale', 1.1);
        %model.component('comp1').mesh('mesh1').feature('ftri1').set('yscale', 1.1);
        model.component('comp1').mesh('mesh1').run('ftri1');
        
        [meshstats,meshdata] = mphmeshstats(model,'mesh1');
        pts=meshdata.vertex;
        pts=pts';
        tri=meshdata.elem{2};
        tri=(tri+1)';
        %trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facealpha',0.5);
        %axis equal;hold on
        clear face_ms;
        for i=1:size(tri,1)
            face_ms(i)={[tri(i,1)-1,tri(i,2)-1,tri(i,3)-1]};
        end
        data.vertex.x=pts(:,1);
        data.vertex.y=pts(:,2);
        data.vertex.z=pts(:,3);
        data.face.vertex_indices=face_ms;
        tmp_filename=['PAL_MS_',num2str(loop_i),'_CHL_',num2str(loop_i_CHL),'.ply'];
        ply_write(data,tmp_filename,'ascii','double');
    catch
    end
    if GLB_nodeleterepair_mode==1
        tmp_filename=['Debug_PAL_MS_',num2str(loop_i),'_CHL_',num2str(loop_i_CHL),'.mph'];
    else
        tmp_filename=['Debug_nodel_PAL_MS_',num2str(loop_i),'_CHL_',num2str(loop_i_CHL),'.mph'];
    end
        mphsave(model, tmp_filename)
        TMP_count_fail=TMP_count_fail+1;
%     end
end
SAVE_COUNT_FAIL_r2(5)=TMP_count_fail;

%%%% 2.6 SPO cell
TMP_count_fail=0;
for TMP_loop_i=1:numel(SAVE_FAIL_LIST_ext_SPO_MS)
    loop_i=SAVE_FAIL_LIST_ext_SPO_MS(TMP_loop_i);
    GLB_count_fail=GLB_count_fail+1;
    
    if GLB_count_fail~=1
        tic
        ModelUtil.remove('Model_3Dleaf_dicot4RT');
        model=mphopen('TMP2_pre_RT_geo_export.mph','Model_3Dleaf_dicot4RT');
        toc
    end
    
    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').init;
    clear del_list
    del_list=[GLB_tag_EPL_l_set, ...
        GLB_tag_EPL_u_set, ...
        GLB_tag_PAL_MS, ...
        GLB_tag_PAL_VAC, ...
        cat(2,GLB_tag_PAL_CHL{:}), ...
        GLB_tag_SPO_MS([1:loop_i-1,loop_i+1:numel(GLB_tag_SPO_MS)]), ...
        GLB_tag_SPO_VAC, ...
        cat(2,GLB_tag_SPO_CHL{:})];
    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(del_list);
    
    model.component('comp1').geom('geom1').feature('csur1').selection('input').set(GLB_tag_SPO_MS(loop_i));
    
    if GLB_nodeleterepair_mode==1
        %%%%%% delete small faces and edges
        model.component('comp1').geom('geom1').defeaturing('SmallFaces').set('entsize', GLB_threshold);
        model.component('comp1').geom('geom1').defeaturing('SmallFaces').find;
        if (model.component('comp1').geom('geom1').defeaturing('SmallFaces').detail().size())~=0
            try
                model.component('comp1').geom('geom1').defeaturing('SmallFaces').deleteAll('dsf1');
            catch
                disp('[ErrorMsg]: cannot delete small faces');
            end
        end
        model.component('comp1').geom('geom1').defeaturing('ShortEdges').set('entsize', GLB_threshold);
        model.component('comp1').geom('geom1').defeaturing('ShortEdges').find;
        if (model.component('comp1').geom('geom1').defeaturing('ShortEdges').detail().size())~=0
            try
                model.component('comp1').geom('geom1').defeaturing('ShortEdges').deleteAll('dse1');
            catch
                disp('[ErrorMsg]: cannot delete short edges');
            end
        end
    end
    
    model.component('comp1').mesh.create('mesh1');
    model.component('comp1').mesh('mesh1').create('size1', 'Size');
    model.component('comp1').mesh('mesh1').feature('size1').set('hauto', 4);
    model.component('comp1').mesh('mesh1').feature('size1').set('hmax', GLB_max_mesh_length_MSVAC);
    model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
    model.component('comp1').mesh('mesh1').feature('ftri1').selection.all;
%     try
%         model.component('comp1').mesh('mesh1').run('ftri1');
%         
%         [meshstats,meshdata] = mphmeshstats(model,'mesh1');
%         pts=meshdata.vertex;
%         pts=pts';
%         tri=meshdata.elem{2};
%         tri=(tri+1)';
%         %trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facealpha',0.5);
%         %axis equal;hold on
%         clear face_ms;
%         for i=1:size(tri,1)
%             face_ms(i)={[tri(i,1)-1,tri(i,2)-1,tri(i,3)-1]};
%         end
%         data.vertex.x=pts(:,1);
%         data.vertex.y=pts(:,2);
%         data.vertex.z=pts(:,3);
%         data.face.vertex_indices=face_ms;
%         tmp_filename=['SPO_MS_',num2str(loop_i),'.ply'];
%         ply_write(data,tmp_filename,'ascii','double');
%     catch
    if GLB_nodeleterepair_mode==1
        tmp_filename=['Debug_SPO_MS_',num2str(loop_i),'.mph'];
    else
        tmp_filename=['Debug_nodel_SPO_MS_',num2str(loop_i),'.mph'];
    end
        mphsave(model, tmp_filename)
        TMP_count_fail=TMP_count_fail+1;
%     end
end
SAVE_COUNT_FAIL_r2(6)=TMP_count_fail;

%%%% 2.7 SPO VAC
TMP_count_fail=0;
for TMP_loop_i=1:numel(SAVE_FAIL_LIST_ext_SPO_VAC)
    loop_i=SAVE_FAIL_LIST_ext_SPO_VAC(TMP_loop_i);
    GLB_count_fail=GLB_count_fail+1;
    
    if GLB_count_fail~=1
        tic
        ModelUtil.remove('Model_3Dleaf_dicot4RT');
        model=mphopen('TMP2_pre_RT_geo_export.mph','Model_3Dleaf_dicot4RT');
        toc
    end
    
    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').init;
    clear del_list
    del_list=[GLB_tag_EPL_l_set, ...
        GLB_tag_EPL_u_set, ...
        GLB_tag_PAL_MS, ...
        GLB_tag_PAL_VAC, ...
        cat(2,GLB_tag_PAL_CHL{:}), ...
        GLB_tag_SPO_MS, ...
        GLB_tag_SPO_VAC([1:loop_i-1,loop_i+1:numel(GLB_tag_SPO_VAC)]), ...
        cat(2,GLB_tag_SPO_CHL{:})];
    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(del_list);
    
    model.component('comp1').geom('geom1').feature('csur1').selection('input').set(GLB_tag_SPO_VAC(loop_i));
    
    if GLB_nodeleterepair_mode==1
        %%%%%% delete small faces and edges
        model.component('comp1').geom('geom1').defeaturing('SmallFaces').set('entsize', GLB_threshold);
        model.component('comp1').geom('geom1').defeaturing('SmallFaces').find;
        if (model.component('comp1').geom('geom1').defeaturing('SmallFaces').detail().size())~=0
            try
                model.component('comp1').geom('geom1').defeaturing('SmallFaces').deleteAll('dsf1');
            catch
                disp('[ErrorMsg]: cannot delete small faces');
            end
        end
        model.component('comp1').geom('geom1').defeaturing('ShortEdges').set('entsize', GLB_threshold);
        model.component('comp1').geom('geom1').defeaturing('ShortEdges').find;
        if (model.component('comp1').geom('geom1').defeaturing('ShortEdges').detail().size())~=0
            try
                model.component('comp1').geom('geom1').defeaturing('ShortEdges').deleteAll('dse1');
            catch
                disp('[ErrorMsg]: cannot delete short edges');
            end
        end
    end
    
    model.component('comp1').mesh.create('mesh1');
    model.component('comp1').mesh('mesh1').create('size1', 'Size');
    model.component('comp1').mesh('mesh1').feature('size1').set('hauto', 4);
    model.component('comp1').mesh('mesh1').feature('size1').set('hmax', GLB_max_mesh_length_MSVAC);
    model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
    model.component('comp1').mesh('mesh1').feature('ftri1').selection.all;
%     try
%         model.component('comp1').mesh('mesh1').run('ftri1');
%         
%         [meshstats,meshdata] = mphmeshstats(model,'mesh1');
%         pts=meshdata.vertex;
%         pts=pts';
%         tri=meshdata.elem{2};
%         tri=(tri+1)';
%         %trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facealpha',0.5);
%         %axis equal;hold on
%         clear face_ms;
%         for i=1:size(tri,1)
%             face_ms(i)={[tri(i,1)-1,tri(i,2)-1,tri(i,3)-1]};
%         end
%         data.vertex.x=pts(:,1);
%         data.vertex.y=pts(:,2);
%         data.vertex.z=pts(:,3);
%         data.face.vertex_indices=face_ms;
%         tmp_filename=['SPO_VAC_',num2str(loop_i),'.ply'];
%         ply_write(data,tmp_filename,'ascii','double');
%     catch
    if GLB_nodeleterepair_mode==1
        tmp_filename=['Debug_SPO_VAC_',num2str(loop_i),'.mph'];
    else
        tmp_filename=['Debug_nodel_SPO_VAC_',num2str(loop_i),'.mph'];
    end
        mphsave(model, tmp_filename)
        TMP_count_fail=TMP_count_fail+1;
%     end
end
SAVE_COUNT_FAIL_r2(7)=TMP_count_fail;

%%%% 2.8 SPO CHL
TMP_count_fail=0;
for TMP_loop_i=1:size(SAVE_FAIL_LIST_ext_SPO_CHL,1)
    loop_i=SAVE_FAIL_LIST_ext_SPO_CHL(TMP_loop_i,1);
    loop_i_CHL=SAVE_FAIL_LIST_ext_SPO_CHL(TMP_loop_i,2);
    GLB_count_fail=GLB_count_fail+1;
    
    if GLB_count_fail~=1
        tic
        ModelUtil.remove('Model_3Dleaf_dicot4RT');
        model=mphopen('TMP2_pre_RT_geo_export.mph','Model_3Dleaf_dicot4RT');
        toc
    end
    
    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').init;
    clear del_list
    GLB_tag_SPO_CHL_tmpdel=GLB_tag_SPO_CHL;
    tmp_list=GLB_tag_SPO_CHL{loop_i};
    %GLB_tag_SPO_CHL_tmpdel{loop_i}=tmp_list([1:loop_i_CHL-1,loop_i_CHL+1:GLB_count_SPO_chl(loop_i)]);
    GLB_tag_SPO_CHL_tmpdel{loop_i}=tmp_list([1:loop_i_CHL-1,loop_i_CHL+1:numel(tmp_list)]);
    del_list=[GLB_tag_EPL_l_set, ...
        GLB_tag_EPL_u_set, ...
        GLB_tag_PAL_MS, ...
        GLB_tag_PAL_VAC, ...
        cat(2,GLB_tag_PAL_CHL{:}), ...
        GLB_tag_SPO_MS, ...
        GLB_tag_SPO_VAC, ...
        cat(2,GLB_tag_SPO_CHL_tmpdel{:})];
    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(del_list);
    
    model.component('comp1').geom('geom1').feature('csur1').selection('input').set(tmp_list(loop_i_CHL));
    
    if GLB_nodeleterepair_mode==1
        %%%%%% delete small faces and edges
        model.component('comp1').geom('geom1').defeaturing('SmallFaces').set('entsize', GLB_threshold);
        model.component('comp1').geom('geom1').defeaturing('SmallFaces').find;
        if (model.component('comp1').geom('geom1').defeaturing('SmallFaces').detail().size())~=0
            try
                model.component('comp1').geom('geom1').defeaturing('SmallFaces').deleteAll('dsf1');
            catch
                disp('[ErrorMsg]: cannot delete small faces');
            end
        end
        model.component('comp1').geom('geom1').defeaturing('ShortEdges').set('entsize', GLB_threshold);
        model.component('comp1').geom('geom1').defeaturing('ShortEdges').find;
        if (model.component('comp1').geom('geom1').defeaturing('ShortEdges').detail().size())~=0
            try
                model.component('comp1').geom('geom1').defeaturing('ShortEdges').deleteAll('dse1');
            catch
                disp('[ErrorMsg]: cannot delete short edges');
            end
        end
    end
    
    model.component('comp1').mesh.create('mesh1');
    model.component('comp1').mesh('mesh1').create('size1', 'Size');
    model.component('comp1').mesh('mesh1').feature('size1').set('hauto', 4);
    model.component('comp1').mesh('mesh1').feature('size1').set('hmax', GLB_max_mesh_length_CHL);
    model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
    model.component('comp1').mesh('mesh1').feature('ftri1').selection.all;
    try
        %model.component('comp1').mesh('mesh1').feature('ftri1').set('xscale', 1.1);
        %model.component('comp1').mesh('mesh1').feature('ftri1').set('yscale', 1.1);
        model.component('comp1').mesh('mesh1').run('ftri1');
        
        [meshstats,meshdata] = mphmeshstats(model,'mesh1');
        pts=meshdata.vertex;
        pts=pts';
        tri=meshdata.elem{2};
        tri=(tri+1)';
        %trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facealpha',0.5);
        %axis equal;hold on
        clear face_ms;
        for i=1:size(tri,1)
            face_ms(i)={[tri(i,1)-1,tri(i,2)-1,tri(i,3)-1]};
        end
        data.vertex.x=pts(:,1);
        data.vertex.y=pts(:,2);
        data.vertex.z=pts(:,3);
        data.face.vertex_indices=face_ms;
        tmp_filename=['SPO_MS_',num2str(loop_i),'_CHL_',num2str(loop_i_CHL),'.ply'];
        ply_write(data,tmp_filename,'ascii','double');
    catch
    end
    if GLB_nodeleterepair_mode==1
        tmp_filename=['Debug_SPO_MS_',num2str(loop_i),'_CHL_',num2str(loop_i_CHL),'.mph'];
    else
        tmp_filename=['Debug_nodel_SPO_MS_',num2str(loop_i),'_CHL_',num2str(loop_i_CHL),'.mph'];
    end
        mphsave(model, tmp_filename)
        TMP_count_fail=TMP_count_fail+1;
%     end
end
SAVE_COUNT_FAIL_r2(8)=TMP_count_fail;

%% ERROR Msg for manually debug
if any(SAVE_COUNT_FAIL_r2)==1
    error('[ERROR Msg] some surfaces cannot be meshed, please check recent Debug_EPL/PAL/SPO_*.mph files manually.')
end