% eLeaf: 3D model of dicot leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 0.1

% PART ONE
% e_geo_PAL_MS_distribute_v0_1.m
% distribute PAL in PAL box to match porosity
%
%% PART TWO %%
% e_geo_PAL_CHL_distribute_v0_1.m
% generate centroid for chloroplast
%
% PART THREE
% e_geo_PAL_MIT_distribute_v0_1.m
% generate centroid for mitochondria

clear all;
load parainput.mat

%%DEBUG
%FLAG_debug_msg=1;
%%END DEBUG

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model_PAL_CHL_distribute');

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

%% PART TWO
model.param.set('PAL_dis_chl2wall', [num2str(PAL_dis_chl2wall),'[m]']);
model.param.set('PAL_thick_chl_layer', [num2str(PAL_thick_chl_layer),'[m]']);
model.param.set('PAL_diameter_chl_patch', [num2str(PAL_diameter_chl_patch),'[m]']);
model.param.set('PAL_radius_chl_sphere', [num2str(PAL_radius_chl_sphere),'[m]']);

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
        tmpstr_displz=['SPO_box_thick+MODEL_ddis+(PAL_MS_radius)*(1-PAL_MS_lcap_flat)+PAL_dis_chl2wall+',num2str(tmp_rand_pos_mat(loop_i,loop_j))];
        eval(['model.component(''comp1'').geom(''geom1'').feature(tmptag_pal_1).set(''pos'', {''',tmpstr_displx,''' ''',tmpstr_disply,''' ''',tmpstr_displz,'''});'])
        model.component('comp1').geom('geom1').feature(tmptag_pal_1).set('r', 'PAL_MS_radius-PAL_dis_chl2wall');
        eval(['model.component(''comp1'').geom(''geom1'').feature(tmptag_pal_1).set(''h'', ''PAL_MS_height-2*PAL_dis_chl2wall-',num2str(delta_cylinder_d_mat(loop_i,loop_j)),''');'])
        
        tmptag_pal_2=['elp', num2str(count_elp)];
        model.component('comp1').geom('geom1').create(tmptag_pal_2, 'Ellipsoid');
        count_elp=count_elp+1;
        eval(['model.component(''comp1'').geom(''geom1'').feature(tmptag_pal_2).set(''pos'', {''',tmpstr_displx,''' ''',tmpstr_disply,''' ''',tmpstr_displz,'''});'])
        model.component('comp1').geom('geom1').feature(tmptag_pal_2).set('semiaxes', {'PAL_MS_radius-PAL_dis_chl2wall' 'PAL_MS_radius-PAL_dis_chl2wall' '(PAL_MS_radius-PAL_dis_chl2wall)*(1-PAL_MS_lcap_flat)'});
        
        tmptag_pal_3=['elp', num2str(count_elp)];
        count_elp=count_elp+1;
        tmpstr_displz=['SPO_box_thick+MODEL_ddis+(PAL_MS_radius)*(1-PAL_MS_lcap_flat)+PAL_MS_height+',num2str(tmp_rand_pos_mat(loop_i,loop_j)),'-PAL_dis_chl2wall-',num2str(delta_cylinder_d_mat(loop_i,loop_j))];
        model.component('comp1').geom('geom1').create(tmptag_pal_3, 'Ellipsoid');
        eval(['model.component(''comp1'').geom(''geom1'').feature(tmptag_pal_3).set(''pos'', {''',tmpstr_displx,''' ''',tmpstr_disply,''' ''',tmpstr_displz,'''});'])
        model.component('comp1').geom('geom1').feature(tmptag_pal_3).set('semiaxes', {'PAL_MS_radius-PAL_dis_chl2wall' 'PAL_MS_radius-PAL_dis_chl2wall' '(PAL_MS_radius-PAL_dis_chl2wall)*(1-PAL_MS_ucap_flat)'});
        
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
            PAL_CHL_tag_set{count_PAL}=tmptag_int;
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
            PAL_CHL_tag_set{count_PAL}=tmptag_csur;
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
PAL_CHL_center_trimesh_pts={};%store chloroplasts centroid for all PAL cells
PAL_CHL_center_trimesh_tri={};
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
    model.component('comp1').mesh(tmptag_mesh).feature(tmptag_size).set('hmin', PAL_diameter_chl_patch/20);
    model.component('comp1').mesh(tmptag_mesh).feature(tmptag_size).set('hmaxactive', true);
    model.component('comp1').mesh(tmptag_mesh).feature(tmptag_size).set('hmax', PAL_diameter_chl_patch/10);%max distance allow
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
    
    PAL_CHL_center_trimesh_pts{loop_i}=pts;
    PAL_CHL_center_trimesh_tri{loop_i}=tri;
end

% % check distance matrix
% for loop_i=1:numel(TagsCreSel)
%     tmp_pts=PAL_CHL_center_trimesh_pts{loop_i};
%     tmp_tri=PAL_CHL_center_trimesh_tri{loop_i};
%     tmpd1=norm(tmp_pts(tmp_tri(loop_i,1),:)-tmp_pts(tmp_tri(loop_i,2),:));
%     tmpd2=norm(tmp_pts(tmp_tri(loop_i,2),:)-tmp_pts(tmp_tri(loop_i,3),:));
%     tmpd3=norm(tmp_pts(tmp_tri(loop_i,3),:)-tmp_pts(tmp_tri(loop_i,1),:));
%     tmpd_all(loop_i,:)=[tmpd1,tmpd2,tmpd3];
% end

% generate centroid from {pts,tri}
PAL_CHL_center_final_pts={};
for loop_i=1:numel(TagsCreSel)
    tmp_pts=PAL_CHL_center_trimesh_pts{loop_i};
    tmp_tri=PAL_CHL_center_trimesh_tri{loop_i};
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
                if any(tmp_dis_sq_all<=PAL_diameter_chl_patch^2)
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
    PAL_CHL_center_final_pts{loop_i}=tmp_sel_pts;
end

% plot sphere; check maximum coverage?
PAL_CHL_sph_all_objname={};
for loop_i=1:numel(TagsCreSel)
    tmp_pts=PAL_CHL_center_final_pts{loop_i};
    %create CHL sphere
    tmptag_set={};
    for loop_j=1:size(tmp_pts,1)
        tmptag_sph=['sph',num2str(count_sph)];count_sph=count_sph+1;
        model.component('comp1').geom('geom1').create(tmptag_sph, 'Sphere');
        model.component('comp1').geom('geom1').feature(tmptag_sph).set('r', 'PAL_radius_chl_sphere');
        model.component('comp1').geom('geom1').feature(tmptag_sph).set('pos', string(tmp_pts(loop_j,:)));%{'5e-6' '0' '150e-6'}
        tmptag_set{loop_j}=tmptag_sph;
    end
    PAL_CHL_sph_all_objname{loop_i}=tmptag_set;
    model.component('comp1').geom('geom1').run(tmptag_sph);
    %merge CHL sphere
    tmptag_uni=['uni', num2str(count_uni)];%count_uni=count_uni+1;
    model.component('comp1').geom('geom1').create(tmptag_uni, 'Union');
    model.component('comp1').geom('geom1').feature(tmptag_uni).set('keep', true);
    model.component('comp1').geom('geom1').feature(tmptag_uni).set('intbnd', false);
    model.component('comp1').geom('geom1').feature(tmptag_uni).selection('input').set(tmptag_set);
    %calculate surface area
    model.component('comp1').geom('geom1').measure.selection.init(2);
    tmp_NFaces=model.geom('geom1').obj(PAL_CHL_tag_set{loop_i}).getNFaces;
    model.component('comp1').geom('geom1').measure.selection.set(PAL_CHL_tag_set{loop_i},1:tmp_NFaces);%set('int1',[1,2,3]);%[1,2,3] is corresponding boundary
    full_area(loop_i)=model.geom('geom1').measure().getVolume();%Remark: getVolume here actually getArea
    %int operation
    tmptag_int=['int',num2str(count_int)];%count_int=count_int+1;
    model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
    model.component('comp1').geom('geom1').feature(tmptag_int).set('keep', false);
    model.component('comp1').geom('geom1').feature(tmptag_int).set('intbnd', false);
    tmpcellstr=[PAL_CHL_tag_set(loop_i),cellstr(tmptag_uni)];
    model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set(tmpcellstr);
    %calculate projected area
    model.component('comp1').geom('geom1').run(tmptag_int);
    model.component('comp1').geom('geom1').measure.selection.init(2);
    tmp_NFaces=model.geom('geom1').obj(tmptag_int).getNFaces;%%alternative method: meshstats - getNEntities()
    model.component('comp1').geom('geom1').measure.selection.set(tmptag_int,1:tmp_NFaces);%set('int1',[1,2,3]);%[1,2,3] is corresponding boundary
    project_area_max(loop_i)=model.geom('geom1').measure().getVolume();
    coverage_max(loop_i)=project_area_max(loop_i)/full_area(loop_i);
    %remove uni and int operation
    model.component('comp1').geom('geom1').feature.remove(tmptag_uni);
    model.component('comp1').geom('geom1').feature.remove(tmptag_int);
end

% compare max coverage with PAL_chl_ScperSm
% match coverage of each cell to PAL_chl_ScperSm
% using a 'Bisection method'
if min(coverage_max)<PAL_chl_ScperSm
    error('[ErrMsg] Cannot meet input PAL_chl_ScperSm.')
end

PAL_CHL_sph_selected_idx={};% store the index of obj for selected chlo to match ScperSm
PAL_CHL_sph_notselected_idx={};
for loop_i=1:numel(TagsCreSel)
    tmp_CHL_sph_all_objname=PAL_CHL_sph_all_objname{loop_i};
    tmp_CHL_sph_selected_idx=1:numel(tmp_CHL_sph_all_objname);
    tmp_CHL_sph_notselected_idx=[];
    tmp_CHL_sph_fullnum=numel(tmp_CHL_sph_all_objname);
    tmp_step=ceil(numel(tmp_CHL_sph_all_objname)/2);
    FLAG_tmp_step=1;
    cur_coverage_cell=coverage_max(loop_i);%initial current coverage
    tmp_record_cur_coverage=[];
    tmp_record_selected_idx={};

    % Bisection method
    %%% phase 1
    SIGN_cur_coverage_cell=sign(cur_coverage_cell-PAL_chl_ScperSm);
    while FLAG_tmp_step==1
        if cur_coverage_cell > PAL_chl_ScperSm
            %random delete tmp_step spheres from selected_idx/sph
            tmp_idx=randperm(numel(tmp_CHL_sph_selected_idx),tmp_step);
            tmp_CHL_sph_selected_idx=tmp_CHL_sph_selected_idx(setdiff(1:numel(tmp_CHL_sph_selected_idx),tmp_idx));
            tmp_CHL_sph_notselected_idx=setdiff(1:tmp_CHL_sph_fullnum,tmp_CHL_sph_selected_idx);
        else %cur_coverage_cell < PAL_chl_ScperSm
            %random delete tmp_step spheres from selected_idx/sph
            tmp_idx=randperm(numel(tmp_CHL_sph_notselected_idx),tmp_step);
            tmp_CHL_sph_notselected_idx=tmp_CHL_sph_notselected_idx(setdiff(1:numel(tmp_CHL_sph_notselected_idx),tmp_idx));
            tmp_CHL_sph_selected_idx=setdiff(1:tmp_CHL_sph_fullnum,tmp_CHL_sph_notselected_idx);
        end
        %calculate coverage again
        %%merge selected CHL sphere
        tmptag_uni=['uni', num2str(count_uni)];%count_uni=count_uni+1;
        model.component('comp1').geom('geom1').create(tmptag_uni, 'Union');
        model.component('comp1').geom('geom1').feature(tmptag_uni).set('intbnd', false);
        tmptag_set=tmp_CHL_sph_all_objname(tmp_CHL_sph_selected_idx);
        model.component('comp1').geom('geom1').feature(tmptag_uni).selection('input').set(tmptag_set);
        %%int operation
        tmptag_int=['int',num2str(count_int)];%count_int=count_int+1;
        model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
        model.component('comp1').geom('geom1').feature(tmptag_int).set('keep', false);
        model.component('comp1').geom('geom1').feature(tmptag_int).set('intbnd', false);
        tmpcellstr=[PAL_CHL_tag_set(loop_i),cellstr(tmptag_uni)];
        model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set(tmpcellstr);
        %%calculate projected area
        model.component('comp1').geom('geom1').run(tmptag_int);
        model.component('comp1').geom('geom1').measure.selection.init(2);
        tmp_NFaces=model.geom('geom1').obj(tmptag_int).getNFaces;%%alternative method: meshstats - getNEntities()
        model.component('comp1').geom('geom1').measure.selection.set(tmptag_int,1:tmp_NFaces);%set('int1',[1,2,3]);%[1,2,3] is corresponding boundary
        tmp_project_area=model.geom('geom1').measure().getVolume();
        cur_coverage_cell=tmp_project_area/full_area(loop_i);
        tmp_record_cur_coverage=[tmp_record_cur_coverage,cur_coverage_cell];
        tmp_record_selected_idx{numel(tmp_record_cur_coverage)}=tmp_CHL_sph_selected_idx;
        %%remove uni and int operation
        model.component('comp1').geom('geom1').feature.remove(tmptag_uni);
        model.component('comp1').geom('geom1').feature.remove(tmptag_int);
        %%debug
        if FLAG_debug_msg==1
            disp(['step size: ',num2str(tmp_step),' pts. ','current coverage: ',num2str(cur_coverage_cell)]);
        end
        %%half tmp_step
        if tmp_step==1
            FLAG_tmp_step=0;%execute one step=1 here in phase1
        end
        tmp_step=floor(tmp_step/2);
    end

    %%% phase III: search +1 or -1 until SIGN_cur_coverage_cell changes
    tmp_SIGN=SIGN_cur_coverage_cell;
    tmp_step=1;
    while tmp_SIGN*SIGN_cur_coverage_cell==1
        if SIGN_cur_coverage_cell==1
            %random delete tmp_step spheres from selected_idx/sph
            tmp_idx=randperm(numel(tmp_CHL_sph_selected_idx),tmp_step);
            tmp_CHL_sph_selected_idx=tmp_CHL_sph_selected_idx(setdiff(1:numel(tmp_CHL_sph_selected_idx),tmp_idx));
            tmp_CHL_sph_notselected_idx=setdiff(1:tmp_CHL_sph_fullnum,tmp_CHL_sph_selected_idx);
        else
            %random delete tmp_step spheres from selected_idx/sph
            tmp_idx=randperm(numel(tmp_CHL_sph_notselected_idx),tmp_step);
            tmp_CHL_sph_notselected_idx=tmp_CHL_sph_notselected_idx(setdiff(1:numel(tmp_CHL_sph_notselected_idx),tmp_idx));
            tmp_CHL_sph_selected_idx=setdiff(1:tmp_CHL_sph_fullnum,tmp_CHL_sph_notselected_idx);
        end
        %calculate coverage again
        %%merge selected CHL sphere
        tmptag_uni=['uni', num2str(count_uni)];%count_uni=count_uni+1;
        model.component('comp1').geom('geom1').create(tmptag_uni, 'Union');
        model.component('comp1').geom('geom1').feature(tmptag_uni).set('intbnd', false);
        tmptag_set=tmp_CHL_sph_all_objname(tmp_CHL_sph_selected_idx);
        model.component('comp1').geom('geom1').feature(tmptag_uni).selection('input').set(tmptag_set);
        %%int operation
        tmptag_int=['int',num2str(count_int)];%count_int=count_int+1;
        model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
        model.component('comp1').geom('geom1').feature(tmptag_int).set('keep', false);
        model.component('comp1').geom('geom1').feature(tmptag_int).set('intbnd', false);
        tmpcellstr=[PAL_CHL_tag_set(loop_i),cellstr(tmptag_uni)];
        model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set(tmpcellstr);
        %%calculate projected area
        model.component('comp1').geom('geom1').run(tmptag_int);
        model.component('comp1').geom('geom1').measure.selection.init(2);
        tmp_NFaces=model.geom('geom1').obj(tmptag_int).getNFaces;%%alternative method: meshstats - getNEntities()
        model.component('comp1').geom('geom1').measure.selection.set(tmptag_int,1:tmp_NFaces);%set('int1',[1,2,3]);%[1,2,3] is corresponding boundary
        tmp_project_area=model.geom('geom1').measure().getVolume();
        cur_coverage_cell=tmp_project_area/full_area(loop_i);
        tmp_SIGN=sign(cur_coverage_cell-PAL_chl_ScperSm);
        tmp_record_cur_coverage=[tmp_record_cur_coverage,cur_coverage_cell];
        tmp_record_selected_idx{numel(tmp_record_cur_coverage)}=tmp_CHL_sph_selected_idx;
        %%remove uni and int operation
        model.component('comp1').geom('geom1').feature.remove(tmptag_uni);
        model.component('comp1').geom('geom1').feature.remove(tmptag_int);
        %%debug
        if FLAG_debug_msg==1
            disp(['step size: ',num2str(tmp_step),' pts. ','current coverage: ',num2str(cur_coverage_cell)]);
        end
    end

    % Compare add/remove last chlo, take the case with closer value to ScperSm
    if FLAG_debug_msg==1
        tmp_record_cur_coverage
    end
    [~,tmp_idx]=min(abs(tmp_record_cur_coverage-PAL_chl_ScperSm));
    PAL_CHL_sph_selected_idx{loop_i}=tmp_record_selected_idx{tmp_idx};
    PAL_CHL_sph_notselected_idx{loop_i}=setdiff(1:tmp_CHL_sph_fullnum,tmp_record_selected_idx{tmp_idx});
    coverage_final(loop_i)=tmp_record_cur_coverage(tmp_idx);
end

% delete sphs, make a clean uni and int for each mesophyll cell
for loop_i=1:numel(TagsCreSel)
    tmp_CHL_sph_all_objname=PAL_CHL_sph_all_objname{loop_i};
    tmp_CHL_sph_selected_idx=PAL_CHL_sph_selected_idx{loop_i};
    tmp_CHL_sph_notselected_idx=PAL_CHL_sph_notselected_idx{loop_i};
    if isempty(tmp_CHL_sph_notselected_idx)==0
        %%remove nonselected spheres
        tmptag_del=['del', num2str(count_del)];count_del=count_del+1;
        model.component('comp1').geom('geom1').create(tmptag_del, 'Delete');
        model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').init;
        model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(tmp_CHL_sph_all_objname(tmp_CHL_sph_notselected_idx));
    end
    %%uni selected CHL spheres
    tmptag_uni=['uni', num2str(count_uni)];count_uni=count_uni+1;
    model.component('comp1').geom('geom1').create(tmptag_uni, 'Union');
    model.component('comp1').geom('geom1').feature(tmptag_uni).set('intbnd', false);
    tmptag_set=tmp_CHL_sph_all_objname(tmp_CHL_sph_selected_idx);
    model.component('comp1').geom('geom1').feature(tmptag_uni).selection('input').set(tmptag_set);
    %%int operation
    tmptag_int=['int',num2str(count_int)];count_int=count_int+1;
    model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
    model.component('comp1').geom('geom1').feature(tmptag_int).set('keep', false);
    model.component('comp1').geom('geom1').feature(tmptag_int).set('intbnd', false);
    tmpcellstr=[PAL_CHL_tag_set(loop_i),cellstr(tmptag_uni)];
    model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set(tmpcellstr);

    model.component('comp1').geom('geom1').run(tmptag_int);
end

% save mat and mph files
save tmp_PAL_CHL_distribute.mat -regexp '^(?!(model|ans)$).'
disp('Potential centroids for chloroplasts in Palisade are ready.')
mphsave(model, 'tmp_PAL_CHL_geo.mph')
ModelUtil.remove('Model_PAL_CHL_distribute');
