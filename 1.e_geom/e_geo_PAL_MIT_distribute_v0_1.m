% eLeaf: 3D model of dicot leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 0.1

% PART ONE
% e_geo_PAL_MS_distribute_v0_1.m
% distribute PAL in PAL box to match porosity
%
% PART TWO
% e_geo_PAL_CH_distribute_v0_1.m
% generate centroid for chloroplast
%
%% PART THREE %%
% e_geo_PAL_MIT_distribute_v0_1.m
% generate centroid for mitochondria

clear all;
load parainput.mat

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model_PAL_MIT_distribute');

model.param.set('PAL_box_thick', [num2str(PAL_box_thick),'[m]']);
model.param.set('PAL_box_length', [num2str(PAL_box_length),'[m]']);
model.param.set('SPO_box_thick', [num2str(SPO_box_thick),'[m]']);
model.param.set('MODEL_ddis', [num2str(MODEL_ddis),'[m]']);
model.param.set('PAL_MS_radius', [num2str(PAL_MS_radius),'[m]']);
model.param.set('PAL_MS_lcap_flat', [num2str(PAL_MS_lcap_flat)]);
model.param.set('PAL_MS_ucap_flat', [num2str(PAL_MS_ucap_flat)]);
model.param.set('PAL_MS_height', 'PAL_box_thick-MODEL_ddis*2-(PAL_MS_radius)*(1-PAL_MS_ucap_flat)-(PAL_MS_radius)*(1-PAL_MS_lcap_flat)');

model.component.create('comp1', true);
model.component('comp1').geom.create('geom1', 3);

model.component('comp1').geom('geom1').create('blk1', 'Block');
model.component('comp1').geom('geom1').feature('blk1').set('pos', {'0' '0' 'SPO_box_thick'});
model.component('comp1').geom('geom1').feature('blk1').set('size', {'PAL_box_length' 'PAL_box_length' 'PAL_box_thick'});

%% PART THREE
model.param.set('PAL_dis_chl2wall', [num2str(PAL_dis_chl2wall),'[m]']);
model.param.set('PAL_thick_chl_layer', [num2str(PAL_thick_chl_layer),'[m]']);
%model.param.set('PAL_diameter_chl_patch', [num2str(PAL_diameter_chl_patch),'[m]']);
%model.param.set('PAL_radius_chl_sphere', [num2str(PAL_radius_chl_sphere),'[m]']);

model.param.set('PAL_dis_mit2chl_rt', num2str(PAL_dis_mit2chl_rt));%relative to PAL_MS_radius
model.param.set('PAL_radius_mit_rt', num2str(PAL_radius_mit_rt));
model.param.set('PAL_dis_vac2mit_rt', num2str(PAL_dis_vac2mit_rt));
model.param.set('PAL_diameter_mit_patch', [num2str(PAL_diameter_mit_patch),'[m]']);

load('tmp_PAL_MS_distribute.mat', 'tmp_rand_pos_mat', 'delta_cylinder_d_mat')

count_cyl=1;
count_elp=1;
count_uni=1;
count_csur=1; %ConvertToSurface
count_int=1;
count_del=1;
count_sph=1;
count_PAL=1;

for loop_i=1:N_PAL_rows
    for loop_j=1:N_PAL_cols
        tmptag_pal_1=['cyl', num2str(count_cyl)];
        model.component('comp1').geom('geom1').create(tmptag_pal_1, 'Cylinder');
        count_cyl=count_cyl+1;
        tmpstr_displx=['(PAL_MS_radius*2+MODEL_ddis)*',num2str(loop_j-1)];
        tmpstr_disply=['(PAL_MS_radius*2+MODEL_ddis)*',num2str(loop_i-1)];
        tmpstr_displz=['SPO_box_thick+MODEL_ddis+(PAL_MS_radius)*(1-PAL_MS_lcap_flat)+PAL_dis_chl2wall+PAL_thick_chl_layer+PAL_MS_radius*(PAL_dis_mit2chl_rt)+',num2str(tmp_rand_pos_mat(loop_i,loop_j))];
        eval(['model.component(''comp1'').geom(''geom1'').feature(tmptag_pal_1).set(''pos'', {''',tmpstr_displx,''' ''',tmpstr_disply,''' ''',tmpstr_displz,'''});'])
        model.component('comp1').geom('geom1').feature(tmptag_pal_1).set('r', 'PAL_MS_radius-PAL_dis_chl2wall-PAL_thick_chl_layer-PAL_MS_radius*(PAL_dis_mit2chl_rt)');
        eval(['model.component(''comp1'').geom(''geom1'').feature(tmptag_pal_1).set(''h'', ''PAL_MS_height-2*PAL_dis_chl2wall-2*PAL_thick_chl_layer-2*PAL_MS_radius*(PAL_dis_mit2chl_rt)-',num2str(delta_cylinder_d_mat(loop_i,loop_j)),''');'])
        
        tmptag_pal_2=['elp', num2str(count_elp)];
        model.component('comp1').geom('geom1').create(tmptag_pal_2, 'Ellipsoid');
        count_elp=count_elp+1;
        eval(['model.component(''comp1'').geom(''geom1'').feature(tmptag_pal_2).set(''pos'', {''',tmpstr_displx,''' ''',tmpstr_disply,''' ''',tmpstr_displz,'''});'])
        model.component('comp1').geom('geom1').feature(tmptag_pal_2).set('semiaxes', {'PAL_MS_radius-PAL_dis_chl2wall-PAL_thick_chl_layer-PAL_MS_radius*(PAL_dis_mit2chl_rt)' 'PAL_MS_radius-PAL_dis_chl2wall-PAL_thick_chl_layer-PAL_MS_radius*(PAL_dis_mit2chl_rt)' '(PAL_MS_radius-PAL_dis_chl2wall-PAL_thick_chl_layer-PAL_MS_radius*(PAL_dis_mit2chl_rt))*(1-PAL_MS_lcap_flat)'});
        
        tmptag_pal_3=['elp', num2str(count_elp)];
        count_elp=count_elp+1;
        tmpstr_displz=['SPO_box_thick+MODEL_ddis+(PAL_MS_radius)*(1-PAL_MS_lcap_flat)+PAL_MS_height+',num2str(tmp_rand_pos_mat(loop_i,loop_j)),'-PAL_dis_chl2wall-PAL_thick_chl_layer-PAL_MS_radius*(PAL_dis_mit2chl_rt)-',num2str(delta_cylinder_d_mat(loop_i,loop_j))];
        model.component('comp1').geom('geom1').create(tmptag_pal_3, 'Ellipsoid');
        eval(['model.component(''comp1'').geom(''geom1'').feature(tmptag_pal_3).set(''pos'', {''',tmpstr_displx,''' ''',tmpstr_disply,''' ''',tmpstr_displz,'''});'])
        model.component('comp1').geom('geom1').feature(tmptag_pal_3).set('semiaxes', {'PAL_MS_radius-PAL_dis_chl2wall-PAL_thick_chl_layer-PAL_MS_radius*(PAL_dis_mit2chl_rt)' 'PAL_MS_radius-PAL_dis_chl2wall-PAL_thick_chl_layer-PAL_MS_radius*(PAL_dis_mit2chl_rt)' '(PAL_MS_radius-PAL_dis_chl2wall-PAL_thick_chl_layer-PAL_MS_radius*(PAL_dis_mit2chl_rt))*(1-PAL_MS_ucap_flat)'});
        
        tmptag_uni=['uni', num2str(count_uni)];
        model.component('comp1').geom('geom1').create(tmptag_uni, 'Union');
        count_uni=count_uni+1;
        model.component('comp1').geom('geom1').feature(tmptag_uni).set('intbnd', false);
        tmptag_set={tmptag_pal_1,tmptag_pal_2,tmptag_pal_3};
        model.component('comp1').geom('geom1').feature(tmptag_uni).selection('input').set(tmptag_set);
        
        % convert obj to surface; integrate with box if it is on boudaries of box;
        % generate selection
        tmpnum_displx=(PAL_MS_radius*2+MODEL_ddis)*(loop_j-1);
        tmpnum_disply=(PAL_MS_radius*2+MODEL_ddis)*(loop_i-1);
        if (tmpnum_displx<PAL_MS_radius) ...,
                || (tmpnum_displx>(PAL_box_length-PAL_MS_radius)) ...,
                || (tmpnum_disply<PAL_MS_radius) ...,
                || (tmpnum_disply>(PAL_box_length-PAL_MS_radius))
            %%% if uni? intesect with box
            tmptag_csur=['csur', num2str(count_csur)];count_csur=count_csur+1;
            model.component('comp1').geom('geom1').create(tmptag_csur, 'ConvertToSurface');
            %model.component('comp1').geom('geom1').feature(tmptag_csur).set('selresult', true);
            model.component('comp1').geom('geom1').feature(tmptag_csur).selection('input').set(tmptag_uni);%[??] tmptag_uni is str; cellstr() is cell;
            %%% intersect csur with box
            tmptag_int=['int', num2str(count_int)];count_int=count_int+1;
            model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
            model.component('comp1').geom('geom1').feature(tmptag_int).set('selresult', true);
            model.component('comp1').geom('geom1').feature(tmptag_int).set('selresultshow', 'bnd');
            model.component('comp1').geom('geom1').feature(tmptag_int).set('keep', true);
            model.component('comp1').geom('geom1').feature(tmptag_int).set('intbnd', false);
            tmpcellstr=[cellstr('blk1'),cellstr(tmptag_csur)];
            model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set(tmpcellstr);%set({'blk1' 'csur1'})
            PAL_MIT_tag_set{count_PAL}=tmptag_int;
            count_PAL=count_PAL+1;
            %%% delete csur
            tmptag_del=['del', num2str(count_del)];count_del=count_del+1;
            model.component('comp1').geom('geom1').create(tmptag_del, 'Delete');
            model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').init;
            model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(tmptag_csur);
        else
            tmptag_csur=['csur', num2str(count_csur)];count_csur=count_csur+1;
            model.component('comp1').geom('geom1').create(tmptag_csur, 'ConvertToSurface');
            model.component('comp1').geom('geom1').feature(tmptag_csur).set('selresult', true);
            model.component('comp1').geom('geom1').feature(tmptag_csur).selection('input').set(tmptag_uni);%[??] tmptag_uni is str; cellstr() is cell;
            PAL_MIT_tag_set{count_PAL}=tmptag_csur;
            count_PAL=count_PAL+1;
        end
    end
end
tmptag_del=['del', num2str(count_del)];count_del=count_del+1;
model.component('comp1').geom('geom1').create(tmptag_del, 'Delete');
model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').init;
model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set({'blk1'});

model.geom('geom1').feature('fin').name('Form Assembly');
model.geom('geom1').feature('fin').set('action', 'assembly');
model.geom('geom1').feature('fin').set('createpairs', 'on');
model.geom('geom1').feature('fin').set('imprint', 'on');
%model.geom('geom1').feature('fin').set('repairtol', '1e-5');
%model.geom('geom1').feature(tmptag_IAS).set('repairtoltype', 'auto');
model.geom('geom1').run;

NBnd=model.geom('geom1').getNBoundaries();% number of boundaries
NCreSel=model.selection.size();% number of created selection, i.e. MS surface
TagsCreSel=model.selection.tags();% almost equal to PAL_tag_set

count_mesh=1;
count_mesh_size=1;
count_mesh_ftri=1;
PAL_MIT_center_trimesh_pts={};%store chloroplasts centroid for all PAL cells
PAL_MIT_center_trimesh_tri={};
for loop_i=1:numel(TagsCreSel)
    tmptag_mesh=['mesh',num2str(count_mesh)];count_mesh=count_mesh+1;
    model.component('comp1').mesh.create(tmptag_mesh);
    
    tmptag_cresel=TagsCreSel(loop_i);
    tmptag_size=['size',num2str(count_mesh_size)];count_mesh_size=count_mesh_size+1;
    model.component('comp1').mesh(tmptag_mesh).create(tmptag_size, 'Size');
    model.component('comp1').mesh(tmptag_mesh).feature(tmptag_size).selection.geom('geom1', 2);
    model.component('comp1').mesh(tmptag_mesh).feature(tmptag_size).selection.named(tmptag_cresel);
    model.component('comp1').mesh(tmptag_mesh).feature(tmptag_size).set('custom', true);
    model.component('comp1').mesh(tmptag_mesh).feature(tmptag_size).set('hminactive', true);
    model.component('comp1').mesh(tmptag_mesh).feature(tmptag_size).set('hmin', PAL_diameter_mit_patch/10);
    model.component('comp1').mesh(tmptag_mesh).feature(tmptag_size).set('hmaxactive', true);
    model.component('comp1').mesh(tmptag_mesh).feature(tmptag_size).set('hmax', PAL_diameter_mit_patch/2);%max distance allow
    tmptag_ftri=['ftri',num2str(count_mesh_ftri)];count_mesh_ftri=count_mesh_ftri+1;
    model.component('comp1').mesh(tmptag_mesh).create(tmptag_ftri, 'FreeTri');
    model.component('comp1').mesh(tmptag_mesh).feature(tmptag_ftri).selection.named(tmptag_cresel);
    model.component('comp1').mesh(tmptag_mesh).run(tmptag_size);
    model.component('comp1').mesh(tmptag_mesh).run(tmptag_ftri);
    
    [meshstats,meshdata] = mphmeshstats(model,tmptag_mesh);
    pts=meshdata.vertex;
    pts=pts';
    tri=meshdata.elem{2};
    tri=(tri+1)';
    %trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facealpha',0.5);
    %axis equal;hold on
    
    PAL_MIT_center_trimesh_pts{loop_i}=pts;
    PAL_MIT_center_trimesh_tri{loop_i}=tri;
end

% % check distance matrix
% for loop_i=1:numel(TagsCreSel)
%     tmp_pts=PAL_MIT_center_trimesh_pts{loop_i};
%     tmp_tri=PAL_MIT_center_trimesh_tri{loop_i};
%     tmpd1=norm(tmp_pts(tmp_tri(loop_i,1),:)-tmp_pts(tmp_tri(loop_i,2),:));
%     tmpd2=norm(tmp_pts(tmp_tri(loop_i,2),:)-tmp_pts(tmp_tri(loop_i,3),:));
%     tmpd3=norm(tmp_pts(tmp_tri(loop_i,3),:)-tmp_pts(tmp_tri(loop_i,1),:));
%     tmpd_all(loop_i,:)=[tmpd1,tmpd2,tmpd3];
% end

% generate centroid from {pts,tri}
PAL_MIT_center_final_pts={};
for loop_i=1:numel(TagsCreSel)
    tmp_pts=PAL_MIT_center_trimesh_pts{loop_i};
    tmp_tri=PAL_MIT_center_trimesh_tri{loop_i};
    clearvars tmp_flag_pts tmp_idx_sel_pts tmp_sel_pts
    tmp_flag_pts=zeros(size(tmp_pts,1),1);%record pts select & pts within certain distance to pts select
    count_pts=1;
    tmp_flag_pts(1)=1;% add pt1 as initial; 2=dead pts, all neighbour are searched; 1=active pts
    tmp_idx_sel_pts(1)=1;
    tmp_sel_pts(1,:)=tmp_pts(1,:);
    num_sel_pts=1;
    while all(tmp_flag_pts==2,'all')==0
        %search pts connect to active pts in tmp_flag_pts
        idx_active_pts=find(tmp_flag_pts==1);%idx set
        cur_idx_active_pts=idx_active_pts(1);
        [tmp_x,tmp_y]=find(tmp_tri==cur_idx_active_pts);
        tmp_x2=unique(reshape(tmp_tri(tmp_x,:),1,[]));
        tmp_x2(tmp_x2==cur_idx_active_pts)=[];
        tmp_x3=setdiff(tmp_x2,find(tmp_flag_pts>=1));
        if isempty(tmp_x3)==1
            tmp_flag_pts(cur_idx_active_pts)=2;
        else
            tmp_flag=0;
            for loop_pts=1:numel(tmp_x3)
                %calculate distance to selected pts
                %if smaller than patch size, flag that pt=1
                %if larger than, mark first pt=2, and break
                %if all smaller, flag cur_idx_active_pts=2
                tmp_pts_cur=tmp_pts(tmp_x3(loop_pts),:);
                tmp_dis_sq_all=sum((tmp_sel_pts-tmp_pts_cur).^2,2);
                if any(tmp_dis_sq_all<=PAL_diameter_mit_patch^2)
                    tmp_flag_pts(tmp_x3(loop_pts))=1;
                else
                    tmp_flag_pts(tmp_x3(loop_pts))=2;
                    num_sel_pts=num_sel_pts+1;
                    tmp_idx_sel_pts(num_sel_pts)=tmp_x3(loop_pts);
                    tmp_sel_pts(num_sel_pts,:)=tmp_pts_cur;
                    tmp_flag=1;
                    break;
                end
            end
            %%if all smaller, flag cur_idx_active_pts=2
            if tmp_flag==0
                tmp_flag_pts(cur_idx_active_pts)=2;
            end
        end
        %tmp_flag_pts;
    end
    PAL_MIT_center_final_pts{loop_i}=tmp_sel_pts;
end

% plot sphere
for loop_i=1:numel(TagsCreSel)
    tmp_pts=PAL_MIT_center_final_pts{loop_i};
    %create MIT sphere
    tmptag_set={};
    for loop_j=1:size(tmp_pts,1)
        tmptag_sph=['sph',num2str(count_sph)];count_sph=count_sph+1;
        model.component('comp1').geom('geom1').create(tmptag_sph, 'Sphere');
        model.component('comp1').geom('geom1').feature(tmptag_sph).set('r', 'PAL_MS_radius*PAL_radius_mit_rt');
        model.component('comp1').geom('geom1').feature(tmptag_sph).set('pos', string(tmp_pts(loop_j,:)));%{'5e-6' '0' '150e-6'}
        tmptag_set{loop_j}=tmptag_sph;
    end
end
model.component('comp1').geom('geom1').run(tmptag_sph);

save tmp_PAL_MIT_distribute.mat -regexp '^(?!(model|ans)$).'
disp('Potential centroids for mitochondria in Palisade are ready.')
mphsave(model, 'tmp_PAL_MIT_geo.mph')
ModelUtil.remove('Model_PAL_MIT_distribute');
