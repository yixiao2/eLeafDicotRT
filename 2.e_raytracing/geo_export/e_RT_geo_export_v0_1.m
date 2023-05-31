% eLeaf: 3D model of dicot leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 0.1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update 2021-May
% - constrain max mesh size
%
% update 2020-Jan
% - add try/catch to e_RT_geo_export_v0_1.m
% - save failed obj during mesh export; [Type, N_cell, N_chl(optional)]
% AUTO_FAIL_CHECK algorithm
% - save intermediate mat&mph results in e_RT_geo_export_v0_1.m
% - pass [failed information, mat, mph] to e_RT_geo_export_autocheckfail_v0_1.m
% - which will try delete the rest geometry and export mesh again
% - if failed again, save mph for manually check

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update 2020-Jan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% begin with the outputs at the middle of e_geo_main
% 1. cut obj with box
% 2. convert domain to surface; and create selection
% 3. create (partially) mesh and export

%% 0. load mat and mph
clear all;
import com.comsol.model.*
import com.comsol.model.util.*
load('tmpsave4RT_geo_export.mat')
model=mphopen('tmpmph4RT_geo_export.mph','Model_3Dleaf_dicot4RT');
disp('[Msg][RT geo export] MAT and MPH files loaded for triangle meshed surface export')

%GLB_count_PAL_chl=[];
%GLB_count_SPO_chl=[];
%GLB_count_nonMS=numel(GLB_tag_EPL_l_set)+numel(GLB_tag_EPL_u_set);

%% 1. cut obj with box
%%%% 1.1 PAL
TMP_tag_set_final_del={};tmp_count=1;
for loop_i=1:Num_pal%numel(GLB_tag_PAL_MS)
    if GLB_mat_contact_PAL_MSVAC_idx(loop_i)==1
        %%%%%% PAL wall
        tmptag_int=['int',num2str(count_int)];count_int=count_int+1;
        model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
        model.component('comp1').geom('geom1').feature(tmptag_int).set('keep', true);
        model.component('comp1').geom('geom1').feature(tmptag_int).set('intbnd', false);
        model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set([cellstr(GLB_tag_PAL_blk), cellstr(GLB_tag_PAL_MS{loop_i})]);
        TMP_tag_set_final_del{tmp_count}=GLB_tag_PAL_MS{loop_i};tmp_count=tmp_count+1;
        GLB_tag_PAL_MS{loop_i}=tmptag_int;
        %%%%%% PAL VAC
        tmptag_int=['int',num2str(count_int)];count_int=count_int+1;
        model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
        model.component('comp1').geom('geom1').feature(tmptag_int).set('keep', true);
        model.component('comp1').geom('geom1').feature(tmptag_int).set('intbnd', false);
        model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set([cellstr(GLB_tag_PAL_blk), cellstr(GLB_tag_PAL_VAC{loop_i})]);
        TMP_tag_set_final_del{tmp_count}=GLB_tag_PAL_VAC{loop_i};tmp_count=tmp_count+1;
        GLB_tag_PAL_VAC{loop_i}=tmptag_int;
        %%%%%% PAL CHL
        tmptag_set_PAL_CHL=GLB_tag_PAL_CHL{loop_i};
        tmp_mat_contact_PAL_CHL_idx=GLB_mat_contact_PAL_CHL_idx{loop_i};
        %GLB_count_PAL_chl(loop_i)=numel(tmptag_set_PAL_CHL);
        for loop_i_CHL=1:numel(tmptag_set_PAL_CHL)
            if tmp_mat_contact_PAL_CHL_idx(loop_i_CHL,loop_i_CHL)==1
                tmptag_int=['int',num2str(count_int)];count_int=count_int+1;
                model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
                model.component('comp1').geom('geom1').feature(tmptag_int).set('keep', true);
                model.component('comp1').geom('geom1').feature(tmptag_int).set('intbnd', false);
                model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set([cellstr(GLB_tag_PAL_blk), cellstr(tmptag_set_PAL_CHL{loop_i_CHL})]);
                TMP_tag_set_final_del{tmp_count}=tmptag_set_PAL_CHL{loop_i_CHL};tmp_count=tmp_count+1;
                tmptag_set_PAL_CHL{loop_i_CHL}=tmptag_int;
            end
        end
        GLB_tag_PAL_CHL{loop_i}=tmptag_set_PAL_CHL;%update tmp to GLB
        
    end
end
%%% delete TMP_tag_set_final_del in each loop
tmptag_del=['del',num2str(count_del)];count_del=count_del+1;
model.component('comp1').geom('geom1').create(tmptag_del, 'Delete');
model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').init;
model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(TMP_tag_set_final_del);

%%%% 1.2 SPO
TMP_tag_set_final_del={};tmp_count=1;
for loop_i=1:Num_spo%numel(GLB_tag_SPO_MS)
    if GLB_mat_contact_SPO_MSVAC_idx(loop_i)==1
        %%%%%% SPO wall
        tmptag_int=['int',num2str(count_int)];count_int=count_int+1;
        model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
        model.component('comp1').geom('geom1').feature(tmptag_int).set('keep', true);
        model.component('comp1').geom('geom1').feature(tmptag_int).set('intbnd', false);
        model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set([cellstr(GLB_tag_SPO_blk), cellstr(GLB_tag_SPO_MS{loop_i})]);
        TMP_tag_set_final_del{tmp_count}=GLB_tag_SPO_MS{loop_i};tmp_count=tmp_count+1;
        GLB_tag_SPO_MS{loop_i}=tmptag_int;
        %%%%%% SPO VAC
        tmptag_int=['int',num2str(count_int)];count_int=count_int+1;
        model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
        model.component('comp1').geom('geom1').feature(tmptag_int).set('keep', true);
        model.component('comp1').geom('geom1').feature(tmptag_int).set('intbnd', false);
        model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set([cellstr(GLB_tag_SPO_blk), cellstr(GLB_tag_SPO_VAC{loop_i})]);
        TMP_tag_set_final_del{tmp_count}=GLB_tag_SPO_VAC{loop_i};tmp_count=tmp_count+1;
        GLB_tag_SPO_VAC{loop_i}=tmptag_int;
        %%%%%% SPO CHL
        tmptag_set_SPO_CHL=GLB_tag_SPO_CHL{loop_i};
        tmp_mat_contact_SPO_CHL_idx=GLB_mat_contact_SPO_CHL_idx{loop_i};
        %GLB_count_SPO_chl(loop_i)=numel(tmptag_set_SPO_CHL);
        for loop_i_CHL=1:numel(tmptag_set_SPO_CHL)
            if tmp_mat_contact_SPO_CHL_idx(loop_i_CHL,loop_i_CHL)==1
                tmptag_int=['int',num2str(count_int)];count_int=count_int+1;
                model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
                model.component('comp1').geom('geom1').feature(tmptag_int).set('keep', true);
                model.component('comp1').geom('geom1').feature(tmptag_int).set('intbnd', false);
                model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set([cellstr(GLB_tag_SPO_blk), cellstr(tmptag_set_SPO_CHL{loop_i_CHL})]);
                TMP_tag_set_final_del{tmp_count}=tmptag_set_SPO_CHL{loop_i_CHL};tmp_count=tmp_count+1;
                tmptag_set_SPO_CHL{loop_i_CHL}=tmptag_int;
            end
        end
        GLB_tag_SPO_CHL{loop_i}=tmptag_set_SPO_CHL;%update tmp to GLB
        
    end
end
%%% delete TMP_tag_set_final_del in each loop
tmptag_del=['del',num2str(count_del)];count_del=count_del+1;
model.component('comp1').geom('geom1').create(tmptag_del, 'Delete');
model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').init;
model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set([TMP_tag_set_final_del,GLB_tag_PAL_blk,GLB_tag_SPO_blk]);
model.geom('geom1').run(tmptag_del);
disp('[Msg][RT geo export] MPH objects trimmed. Start to convert objs to surfaces')

save TMP_matlab.mat -regexp '^(?!(model|ans)$).'
mphsave(model, 'TMP_pre_RT_geo_export.mph')

count_csur=1;
%% 2. convert domain to surface; and create selection
%%%% 2.1 lower epidermis
GLB_tag_EPL_l_set_csur={};
for loop_i=1:numel(GLB_tag_EPL_l_set)
    tmptag_csur=['csur', num2str(count_csur)];count_csur=count_csur+1;
    model.component('comp1').geom('geom1').create(tmptag_csur, 'ConvertToSurface');
    model.component('comp1').geom('geom1').feature(tmptag_csur).set('selresult', true);
    model.component('comp1').geom('geom1').feature(tmptag_int).set('selresultshow', 'bnd');
    model.component('comp1').geom('geom1').feature(tmptag_csur).selection('input').set(GLB_tag_EPL_l_set{loop_i});
    GLB_tag_EPL_l_set_csur{loop_i}=tmptag_csur;
end
disp('[Msg][RT geo export] ')

%%%% 2.2 upper epidermis
GLB_tag_EPL_u_set_csur={};
for loop_i=1:numel(GLB_tag_EPL_u_set)
    tmptag_csur=['csur', num2str(count_csur)];count_csur=count_csur+1;
    model.component('comp1').geom('geom1').create(tmptag_csur, 'ConvertToSurface');
    model.component('comp1').geom('geom1').feature(tmptag_csur).set('selresult', true);
    model.component('comp1').geom('geom1').feature(tmptag_int).set('selresultshow', 'bnd');
    model.component('comp1').geom('geom1').feature(tmptag_csur).selection('input').set(GLB_tag_EPL_u_set{loop_i});
    GLB_tag_EPL_u_set_csur{loop_i}=tmptag_csur;
end

%%%% 2.3-2.5 PAL
GLB_tag_PAL_MS_csur={};
GLB_tag_PAL_VAC_csur={};
GLB_tag_PAL_CHL_csur={};%TYPE=cell array of cell array
for loop_i=1:Num_pal%numel(GLB_tag_PAL_MS)
    %%%% 2.3 PAL cell
    tmptag_csur=['csur', num2str(count_csur)];count_csur=count_csur+1;
    model.component('comp1').geom('geom1').create(tmptag_csur, 'ConvertToSurface');
    model.component('comp1').geom('geom1').feature(tmptag_csur).set('selresult', true);
    model.component('comp1').geom('geom1').feature(tmptag_int).set('selresultshow', 'bnd');
    model.component('comp1').geom('geom1').feature(tmptag_csur).selection('input').set(GLB_tag_PAL_MS{loop_i});
    GLB_tag_PAL_MS_csur{loop_i}=tmptag_csur;
    %%%% 2.4 PAL VAC
    tmptag_csur=['csur', num2str(count_csur)];count_csur=count_csur+1;
    model.component('comp1').geom('geom1').create(tmptag_csur, 'ConvertToSurface');
    model.component('comp1').geom('geom1').feature(tmptag_csur).set('selresult', true);
    model.component('comp1').geom('geom1').feature(tmptag_int).set('selresultshow', 'bnd');
    model.component('comp1').geom('geom1').feature(tmptag_csur).selection('input').set(GLB_tag_PAL_VAC{loop_i});
    GLB_tag_PAL_VAC_csur{loop_i}=tmptag_csur;
    %%%% 2.5 PAL CHL
    tmptag_set_PAL_CHL_csur={};
    tmptag_set_PAL_CHL=GLB_tag_PAL_CHL{loop_i};
    for loop_i_CHL=1:numel(tmptag_set_PAL_CHL)
        tmptag_csur=['csur', num2str(count_csur)];count_csur=count_csur+1;
        model.component('comp1').geom('geom1').create(tmptag_csur, 'ConvertToSurface');
        model.component('comp1').geom('geom1').feature(tmptag_csur).set('selresult', true);
        model.component('comp1').geom('geom1').feature(tmptag_int).set('selresultshow', 'bnd');
        model.component('comp1').geom('geom1').feature(tmptag_csur).selection('input').set(tmptag_set_PAL_CHL{loop_i_CHL});
        tmptag_set_PAL_CHL_csur{loop_i_CHL}=tmptag_csur;
    end
    GLB_tag_PAL_CHL_csur{loop_i}=tmptag_set_PAL_CHL_csur;
end

%%%% 2.6-2.8 SPO
GLB_tag_SPO_MS_csur={};
GLB_tag_SPO_VAC_csur={};
GLB_tag_SPO_CHL_csur={};%TYPE=cell array of cell array
for loop_i=1:numel(GLB_tag_SPO_MS)
    %%%% 2.6 SPO cell
    tmptag_csur=['csur', num2str(count_csur)];count_csur=count_csur+1;
    model.component('comp1').geom('geom1').create(tmptag_csur, 'ConvertToSurface');
    model.component('comp1').geom('geom1').feature(tmptag_csur).set('selresult', true);
    model.component('comp1').geom('geom1').feature(tmptag_int).set('selresultshow', 'bnd');
    model.component('comp1').geom('geom1').feature(tmptag_csur).selection('input').set(GLB_tag_SPO_MS{loop_i});
    GLB_tag_SPO_MS_csur{loop_i}=tmptag_csur;
    %%%% 2.7 SPO VAC
    tmptag_csur=['csur', num2str(count_csur)];count_csur=count_csur+1;
    model.component('comp1').geom('geom1').create(tmptag_csur, 'ConvertToSurface');
    model.component('comp1').geom('geom1').feature(tmptag_csur).set('selresult', true);
    model.component('comp1').geom('geom1').feature(tmptag_int).set('selresultshow', 'bnd');
    model.component('comp1').geom('geom1').feature(tmptag_csur).selection('input').set(GLB_tag_SPO_VAC{loop_i});
    GLB_tag_SPO_VAC_csur{loop_i}=tmptag_csur;
    %%%% 2.8 SPO CHL
    tmptag_set_SPO_CHL_csur={};
    tmptag_set_SPO_CHL=GLB_tag_SPO_CHL{loop_i};
    for loop_i_CHL=1:numel(tmptag_set_SPO_CHL)
        tmptag_csur=['csur', num2str(count_csur)];count_csur=count_csur+1;
        model.component('comp1').geom('geom1').create(tmptag_csur, 'ConvertToSurface');
        model.component('comp1').geom('geom1').feature(tmptag_csur).set('selresult', true);
        model.component('comp1').geom('geom1').feature(tmptag_int).set('selresultshow', 'bnd');
        model.component('comp1').geom('geom1').feature(tmptag_csur).selection('input').set(tmptag_set_SPO_CHL{loop_i_CHL});
        tmptag_set_SPO_CHL_csur{loop_i_CHL}=tmptag_csur;
    end
    GLB_tag_SPO_CHL_csur{loop_i}=tmptag_set_SPO_CHL_csur;
end

%% 2.5 delete small faces before form assembly (from CAD module)
disp('[Msg][RT geo export] Start to delete small faces.')
GLB_threshold=1e-8;

model.component('comp1').geom('geom1').defeaturing('SmallFaces').set('entsize', GLB_threshold);
model.component('comp1').geom('geom1').defeaturing('SmallFaces').find;
model.component('comp1').geom('geom1').defeaturing('SmallFaces').deleteAll('dsf1');
model.component('comp1').geom('geom1').defeaturing('ShortEdges').set('entsize', GLB_threshold);
model.component('comp1').geom('geom1').defeaturing('ShortEdges').find;
model.component('comp1').geom('geom1').defeaturing('ShortEdges').deleteAll('dse1');

model.geom('geom1').feature('fin').name('Form Assembly');
model.geom('geom1').feature('fin').set('action', 'assembly');
model.component('comp1').geom('geom1').feature('fin').set('imprint', false);
model.component('comp1').geom('geom1').feature('fin').set('createpairs', false);
model.geom('geom1').run;
disp('[Msg][RT geo export] Form assembly. Start to export ply mesh files.')

%% 3. create mesh and export
TagsCreSel=model.selection.tags();
count_mesh=1;

%%%% initialize SAVE lists for AUTO_FAIL_CHECK
SAVE_FAIL_LIST_EPL_l=[];%No.1
SAVE_FAIL_LIST_EPL_u=[];
SAVE_FAIL_LIST_PAL_MS=[];
SAVE_FAIL_LIST_PAL_VAC=[];
SAVE_FAIL_LIST_PAL_CHL=[];
SAVE_FAIL_LIST_SPO_MS=[];
SAVE_FAIL_LIST_SPO_VAC=[];
SAVE_FAIL_LIST_SPO_CHL=[];%No.8
SAVE_COUNT_FAIL=zeros(1,8);

%%%% constrain max mesh size
GLB_max_mesh_length_CHL=min([PAL_radius_chl_sphere*2,PAL_diameter_chl_patch,SPO_radius_chl_sphere*2,SPO_diameter_chl_patch])/5;
GLB_max_mesh_length_MSVAC=min([PAL_radius_chl_sphere*2,PAL_diameter_chl_patch,SPO_radius_chl_sphere*2,SPO_diameter_chl_patch])/2;

%%%% 3.1 lower epidermis
TMP_count_fail=0;
for loop_i=1:numel(GLB_tag_EPL_l_set)
    tmptag_cresel=TagsCreSel(count_mesh);
    count_mesh=count_mesh+1;
    if (count_mesh-1)==1
        model.component('comp1').mesh.create('mesh1');
        model.component('comp1').mesh('mesh1').create('size1', 'Size');
        model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
    end
    model.component('comp1').mesh('mesh1').feature('size1').selection.geom('geom1', 2);
    model.component('comp1').mesh('mesh1').feature('size1').selection.named(tmptag_cresel);
    model.component('comp1').mesh('mesh1').feature('size1').set('hauto', 4);
    %model.component('comp1').mesh('mesh1').feature('size1').set('hmax', GLB_max_mesh_length);
    model.component('comp1').mesh('mesh1').feature('ftri1').selection.named(tmptag_cresel);
    model.component('comp1').mesh('mesh1').run('size1');
    try
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
        tmp_filename=['EPL_l_',num2str(loop_i),'.ply'];
        ply_write(data,tmp_filename,'ascii','double');
    catch
        %tmp_filename=['Debug_EPL_u_',num2str(loop_i),'.mph'];
        %mphsave(model, tmp_filename)
        TMP_count_fail=TMP_count_fail+1;
        SAVE_FAIL_LIST_EPL_l(TMP_count_fail)=loop_i;
    end
end
SAVE_COUNT_FAIL(1)=TMP_count_fail;
disp('[Msg][RT geo export] Successfully export all lower epidermis cells.')

%%%% 3.2 upper epidermis
TMP_count_fail=0;
for loop_i=1:numel(GLB_tag_EPL_u_set)
    tmptag_cresel=TagsCreSel(count_mesh);
    count_mesh=count_mesh+1;
%     if (count_mesh-1)==1
%         model.component('comp1').mesh.create('mesh1');
%         model.component('comp1').mesh('mesh1').create('size1', 'Size');
%         model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
%     end
    model.component('comp1').mesh('mesh1').feature('size1').selection.geom('geom1', 2);
    model.component('comp1').mesh('mesh1').feature('size1').selection.named(tmptag_cresel);
    model.component('comp1').mesh('mesh1').feature('size1').set('hauto', 4);
    %model.component('comp1').mesh('mesh1').feature('size1').set('hmax', GLB_max_mesh_length);
    model.component('comp1').mesh('mesh1').feature('ftri1').selection.named(tmptag_cresel);
    model.component('comp1').mesh('mesh1').run('size1');
    try
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
        tmp_filename=['EPL_u_',num2str(loop_i),'.ply'];
        ply_write(data,tmp_filename,'ascii','double');
    catch
        %tmp_filename=['Debug_EPL_u_',num2str(loop_i),'.mph'];
        %mphsave(model, tmp_filename)
        TMP_count_fail=TMP_count_fail+1;
        SAVE_FAIL_LIST_EPL_u(TMP_count_fail)=loop_i;
    end
end
SAVE_COUNT_FAIL(2)=TMP_count_fail;
disp('[Msg][RT geo export] Successfully export all upper epidermis cells.')

%%%% 3.3-3.5 PAL
TMP_count_fail_MS=0;
TMP_count_fail_VAC=0;
TMP_count_fail_CHL=0;
%GLB_PAL_ftri_zscale=1/(min(10,PAL_box_thick/(PAL_MS_radius*2)));
%GLB_PAL_CHL_ftri_zscale=0.5;%potentially decrease the fail rate, but not
%working here in COMSOL 5.5.
for loop_i=1:numel(GLB_tag_PAL_MS)
    %%%% 3.3 PAL MS
    tmptag_cresel=TagsCreSel(count_mesh);
    count_mesh=count_mesh+1;
    model.component('comp1').mesh('mesh1').feature('size1').selection.geom('geom1', 2);
    model.component('comp1').mesh('mesh1').feature('size1').selection.named(tmptag_cresel);
    model.component('comp1').mesh('mesh1').feature('size1').set('hauto', 4);
    model.component('comp1').mesh('mesh1').feature('size1').set('hmax', GLB_max_mesh_length_MSVAC);
    model.component('comp1').mesh('mesh1').feature('ftri1').selection.named(tmptag_cresel);
    %model.component('comp1').mesh('mesh1').feature('ftri1').set('zscale', GLB_PAL_ftri_zscale);
    model.component('comp1').mesh('mesh1').run('size1');
    try
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
        tmp_filename=['PAL_MS_',num2str(loop_i),'.ply'];
        ply_write(data,tmp_filename,'ascii','double');
    catch
        %tmp_filename=['Debug_PAL_MS_',num2str(loop_i),'.mph'];
        %mphsave(model, tmp_filename)
        TMP_count_fail_MS=TMP_count_fail_MS+1;
        SAVE_FAIL_LIST_PAL_MS(TMP_count_fail_MS)=loop_i;
    end
    
    %%%% 3.4 PAL VAC
    tmptag_cresel=TagsCreSel(count_mesh);
    count_mesh=count_mesh+1;
    model.component('comp1').mesh('mesh1').feature('size1').selection.geom('geom1', 2);
    model.component('comp1').mesh('mesh1').feature('size1').selection.named(tmptag_cresel);
    model.component('comp1').mesh('mesh1').feature('size1').set('hauto', 4);
    model.component('comp1').mesh('mesh1').feature('size1').set('hmax', GLB_max_mesh_length_MSVAC);
    model.component('comp1').mesh('mesh1').feature('ftri1').selection.named(tmptag_cresel);
    %model.component('comp1').mesh('mesh1').feature('ftri1').set('zscale', GLB_PAL_ftri_zscale);
    model.component('comp1').mesh('mesh1').run('size1');
    try
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
        tmp_filename=['PAL_VAC_',num2str(loop_i),'.ply'];
        ply_write(data,tmp_filename,'ascii','double');
    catch
        %tmp_filename=['Debug_PAL_VAC_',num2str(loop_i),'.mph'];
        %mphsave(model, tmp_filename)
        TMP_count_fail_VAC=TMP_count_fail_VAC+1;
        SAVE_FAIL_LIST_PAL_VAC(TMP_count_fail_VAC)=loop_i;
    end
    
    %%%% 3.5 PAL CHL
    tmptag_set_PAL_CHL=GLB_tag_PAL_CHL{loop_i};
    for loop_i_CHL=1:numel(tmptag_set_PAL_CHL)
        tmptag_cresel=TagsCreSel(count_mesh);
        count_mesh=count_mesh+1;
        model.component('comp1').mesh('mesh1').feature('size1').selection.geom('geom1', 2);
        model.component('comp1').mesh('mesh1').feature('size1').selection.named(tmptag_cresel);
        model.component('comp1').mesh('mesh1').feature('size1').set('hauto', 4);
        model.component('comp1').mesh('mesh1').feature('size1').set('hmax', GLB_max_mesh_length_CHL);
        model.component('comp1').mesh('mesh1').feature('ftri1').selection.named(tmptag_cresel);
        %model.component('comp1').mesh('mesh1').feature('ftri1').set('zscale', GLB_PAL_CHL_ftri_zscale);
        model.component('comp1').mesh('mesh1').run('size1');
        try
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
            %tmp_filename=['Debug_PAL_MS_',num2str(loop_i),'_CHL_',num2str(loop_i_CHL),'.mph'];
            %mphsave(model, tmp_filename)
            TMP_count_fail_CHL=TMP_count_fail_CHL+1;
            SAVE_FAIL_LIST_PAL_CHL(TMP_count_fail_CHL,1)=loop_i;
            SAVE_FAIL_LIST_PAL_CHL(TMP_count_fail_CHL,2)=loop_i_CHL;
        end
    end
end
SAVE_COUNT_FAIL(3)=TMP_count_fail_MS;
SAVE_COUNT_FAIL(4)=TMP_count_fail_VAC;
SAVE_COUNT_FAIL(5)=TMP_count_fail_CHL;
disp('[Msg][RT geo export] Successfully export all palisade cells/chloroplasts/vacuoles.')

%%%% 3.6-3.8 SPO
TMP_count_fail_MS=0;
TMP_count_fail_VAC=0;
TMP_count_fail_CHL=0;
for loop_i=1:numel(GLB_tag_SPO_MS)
    %%%% 3.3 SPO MS
    tmptag_cresel=TagsCreSel(count_mesh);
    count_mesh=count_mesh+1;
    model.component('comp1').mesh('mesh1').feature('size1').selection.geom('geom1', 2);
    model.component('comp1').mesh('mesh1').feature('size1').selection.named(tmptag_cresel);
    model.component('comp1').mesh('mesh1').feature('size1').set('hauto', 4);
    model.component('comp1').mesh('mesh1').feature('size1').set('hmax', GLB_max_mesh_length_MSVAC);
    model.component('comp1').mesh('mesh1').feature('ftri1').selection.named(tmptag_cresel);
    model.component('comp1').mesh('mesh1').run('size1');
    try
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
        tmp_filename=['SPO_MS_',num2str(loop_i),'.ply'];
        ply_write(data,tmp_filename,'ascii','double');
    catch
        %tmp_filename=['Debug_SPO_MS_',num2str(loop_i),'.mph'];
        %mphsave(model, tmp_filename)
        TMP_count_fail_MS=TMP_count_fail_MS+1;
        SAVE_FAIL_LIST_SPO_MS(TMP_count_fail_MS)=loop_i;
    end
    
    %%%% 3.4 SPO VAC
    tmptag_cresel=TagsCreSel(count_mesh);
    count_mesh=count_mesh+1;
    model.component('comp1').mesh('mesh1').feature('size1').selection.geom('geom1', 2);
    model.component('comp1').mesh('mesh1').feature('size1').selection.named(tmptag_cresel);
    model.component('comp1').mesh('mesh1').feature('size1').set('hauto', 4);
    model.component('comp1').mesh('mesh1').feature('size1').set('hmax', GLB_max_mesh_length_MSVAC);
    model.component('comp1').mesh('mesh1').feature('ftri1').selection.named(tmptag_cresel);
    model.component('comp1').mesh('mesh1').run('size1');
    try
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
        tmp_filename=['SPO_VAC_',num2str(loop_i),'.ply'];
        ply_write(data,tmp_filename,'ascii','double');
    catch
        %tmp_filename=['Debug_SPO_VAC_',num2str(loop_i),'.mph'];
        %mphsave(model, tmp_filename)
        TMP_count_fail_VAC=TMP_count_fail_VAC+1;
        SAVE_FAIL_LIST_SPO_VAC(TMP_count_fail_VAC)=loop_i;
    end
    
    %%%% 3.5 SPO CHL
    tmptag_set_SPO_CHL=GLB_tag_SPO_CHL{loop_i};
    for loop_i_CHL=1:numel(tmptag_set_SPO_CHL)
        tmptag_cresel=TagsCreSel(count_mesh);
        count_mesh=count_mesh+1;
        model.component('comp1').mesh('mesh1').feature('size1').selection.geom('geom1', 2);
        model.component('comp1').mesh('mesh1').feature('size1').selection.named(tmptag_cresel);
        model.component('comp1').mesh('mesh1').feature('size1').set('hauto', 4);
        model.component('comp1').mesh('mesh1').feature('size1').set('hmax', GLB_max_mesh_length_CHL);
        model.component('comp1').mesh('mesh1').feature('ftri1').selection.named(tmptag_cresel);
        model.component('comp1').mesh('mesh1').run('size1');
        try
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
            %tmp_filename=['Debug_SPO_MS_',num2str(loop_i),'_CHL_',num2str(loop_i_CHL),'.mph'];
            %mphsave(model, tmp_filename)
            TMP_count_fail_CHL=TMP_count_fail_CHL+1;
            SAVE_FAIL_LIST_SPO_CHL(TMP_count_fail_CHL,1)=loop_i;
            SAVE_FAIL_LIST_SPO_CHL(TMP_count_fail_CHL,2)=loop_i_CHL;
        end
    end
end
SAVE_COUNT_FAIL(6)=TMP_count_fail_MS;
SAVE_COUNT_FAIL(7)=TMP_count_fail_VAC;
SAVE_COUNT_FAIL(8)=TMP_count_fail_CHL;
disp('[Msg][RT geo export] Successfully export all spongy cells/chloroplasts/vacuoles.')

% load SAVE_e_geom2.mat chl_con_PAL_4RT chl_con_SPO_4RT
% fid=fopen('count_chl4RT','w');
% fprintf(fid,'%d %d\n', Num_pal, Num_spo);
% fprintf(fid,'%d %e\n',[GLB_count_PAL_chl;chl_con_PAL_4RT*ones(1,Num_pal)]);
% fprintf(fid,'%d %e\n',[GLB_count_SPO_chl;chl_con_SPO_4RT*ones(1,Num_spo)]);
% fclose(fid);
% 
save SAVE_RT_geo_export.mat -regexp '^(?!(model|ans)$).'
mphsave(model, 'SAVE_RT_geo_export.mph')
