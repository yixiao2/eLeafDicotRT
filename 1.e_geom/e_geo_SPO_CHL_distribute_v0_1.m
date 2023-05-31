% eLeaf: 3D model of dicot leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 0.1

% PART ONE
% e_geo_SPO_MS_distribute_v0_1.m
% distribute SPO in SPO box to match porosity
%
%% PART TWO %%
% e_geo_SPO_CHL_distribute_v0_1.m
% generate centroid for chloroplast
%
% PART THREE
% e_geo_SPO_MT_distribute_v0_1.m
% generate centroid for mitochondria

clear all;
load parainput.mat

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model_SPO_CHL_distribute');

model.param.set('SPO_box_thick', [num2str(SPO_box_thick),'[m]']);
model.param.set('SPO_box_length', [num2str(SPO_box_length),'[m]']);
model.param.set('SPO_MS_radius', [num2str(SPO_MS_radius),'[m]']);

model.component.create('comp1', true);
model.component('comp1').geom.create('geom1', 3);

model.component('comp1').geom('geom1').create('blk2', 'Block');
model.component('comp1').geom('geom1').feature('blk2').set('pos', [0 0 0]);
model.component('comp1').geom('geom1').feature('blk2').set('size', {'SPO_box_length' 'SPO_box_length' 'SPO_box_thick'});

model.param.set('SPO_dis_chl2wall', [num2str(SPO_dis_chl2wall),'[m]']);
model.param.set('SPO_thick_chl_layer', [num2str(SPO_thick_chl_layer),'[m]']);
model.param.set('SPO_diameter_chl_patch', [num2str(SPO_diameter_chl_patch),'[m]']);
model.param.set('SPO_radius_chl_sphere', [num2str(SPO_radius_chl_sphere),'[m]']);

load('tmp_SPO_MS_distribute.mat', 'SPO_MS_center_set_X', 'SPO_MS_center_set_Y', 'SPO_MS_center_set_Z', 'dis_SPO_MS_center')

Num_spo=numel(SPO_MS_center_set_X);
count_sph=1;
count_dif=1;
count_list_SPO_MS=0;
for loop=1:Num_spo
    tmp_center_X=SPO_MS_center_set_X(loop);
    tmp_center_Y=SPO_MS_center_set_Y(loop);
    tmp_center_Z=SPO_MS_center_set_Z(loop);
    tmptag_spo=['sph',num2str(count_sph)];
    model.component('comp1').geom('geom1').create(tmptag_spo, 'Sphere');
    count_sph=count_sph+1;
    model.component('comp1').geom('geom1').feature(tmptag_spo).set('r', 'SPO_MS_radius-SPO_dis_chl2wall');
    tmpstr_pos={num2str(tmp_center_X), num2str(tmp_center_Y), num2str(tmp_center_Z)};
    model.component('comp1').geom('geom1').feature(tmptag_spo).set('pos', tmpstr_pos);
    
    count_list_SPO_MS=count_list_SPO_MS+1;
    SPO_tag_set{count_list_SPO_MS}=tmptag_spo;
    list_SPO_MS(count_list_SPO_MS)=cellstr(tmptag_spo);%list_SPO_MS record obj name of SPO MS
    %%[to update] SPO_tag_set is the same to list_SPO_MS?
end
% tmptag_spo=['dif',num2str(count_dif)];
% model.component('comp1').geom('geom1').create(tmptag_spo, 'Difference');
% %count_dif=count_dif+1;
% model.component('comp1').geom('geom1').feature(tmptag_spo).selection('input').set({'blk2'});
% model.component('comp1').geom('geom1').feature(tmptag_spo).selection('input2').set(SPO_tag_set);
% model.component('comp1').geom('geom1').feature(tmptag_spo).set('intbnd', false);
% model.component('comp1').geom('geom1').feature(tmptag_spo).set('keep', true);
% %%calculate tmp_vol_IAS_SPO
% model.component('comp1').geom('geom1').run(tmptag_spo);
% model.component('comp1').geom('geom1').measure.selection.init(3);
% model.component('comp1').geom('geom1').measure.selection.set(tmptag_spo, 1);
% tmp_vol_IAS_SPO = model.geom('geom1').measure().getVolume()/(SPO_box_thick*SPO_box_length*SPO_box_length)


%%%create MS with contact area, based on center_set and radius_set
%%%%% store contact information into a matrix mat_contact_cell_idx;
%%%%% if ith cell overlaps jth cell, then mat_contact_cell_idx(i,j)=1
%%%%% if ith cell overlaps with box, then mat_contact_cell_idx(i,i)=1
%%%%% mat_contact_cell_idx is symetric
mat_contact_cell_idx=zeros(Num_spo);% initial
tmp_radius_2=(SPO_MS_radius)^2;
for loop_i=1:Num_spo
    for loop_j=loop_i:Num_spo
        if loop_j==loop_i
            % contact with box?
            if (SPO_MS_center_set_X(loop_i)<SPO_MS_radius) ...,
                    || (SPO_MS_center_set_X(loop_i)>SPO_box_length-(SPO_MS_radius)) ...,
                    || (SPO_MS_center_set_Y(loop_i)<SPO_MS_radius) ...,
                    || (SPO_MS_center_set_Y(loop_i)>SPO_box_length-(SPO_MS_radius)) ...,
                    || (SPO_MS_center_set_Z(loop_i)<SPO_MS_radius) ...,
                    || (SPO_MS_center_set_Z(loop_i)>SPO_box_thick-(SPO_MS_radius))
                mat_contact_cell_idx(loop_i,loop_i)=1;
                %%97/131 contact with box in test_demo
            end
        else
            if ((SPO_MS_center_set_X(loop_j)-SPO_MS_center_set_X(loop_i))^2+(SPO_MS_center_set_Y(loop_j)-SPO_MS_center_set_Y(loop_i))^2 ...,
                +(SPO_MS_center_set_Z(loop_j)-SPO_MS_center_set_Z(loop_i))^2) < tmp_radius_2*4 % (2*radius)^2
                mat_contact_cell_idx(loop_i,loop_j)=1;
                mat_contact_cell_idx(loop_j,loop_i)=1;
            end
        end
    end
end

%%%process each cell with information from mat_contact_cell_idx
count_wp=1;%work plane
count_par=1;%operation partition
count_spl=1;%
count_del=1;
count_uni=1;
count_csur=1;
count_int=1;
%count_sph=1;

for loop_i=1:Num_spo
    tmp_idx_contact=find(mat_contact_cell_idx(loop_i,:)==1);
    for tmp_loop_j=1:numel(tmp_idx_contact)
        loop_j=tmp_idx_contact(tmp_loop_j);
        if loop_j > loop_i
            tmp_center_a=[SPO_MS_center_set_X(loop_i),SPO_MS_center_set_Y(loop_i),SPO_MS_center_set_Z(loop_i)];
            tmp_center_b=[SPO_MS_center_set_X(loop_j),SPO_MS_center_set_Y(loop_j),SPO_MS_center_set_Z(loop_j)];
            tmp_pt1=(tmp_center_a+tmp_center_b)/2;
            tmp_vec1=tmp_center_a-tmp_center_b;% (a,b,c) --> (c,c,-a-b) or (-b-c,a,a)
            if tmp_vec1(3) ~= 0
                tmp_vec2=[tmp_vec1(3),tmp_vec1(3),-tmp_vec1(1)-tmp_vec1(2)];
            else
                tmp_vec2=[-tmp_vec1(2)-tmp_vec1(3),tmp_vec1(1),tmp_vec1(1)];
            end
            tmp_vec3=cross(tmp_vec1/norm(tmp_vec1),tmp_vec2/norm(tmp_vec2))*norm(tmp_vec2);
            tmp_pt2=tmp_pt1+tmp_vec2;
            tmp_pt3=tmp_pt1+tmp_vec3;
            %%%for SPO CHL
            %%%shift wp toward/backward by distance SPO_dis_chl2wall
            tmp_pt1_a=tmp_pt1+tmp_vec1/norm(tmp_vec1)*SPO_dis_chl2wall;
            tmp_pt1_b=tmp_pt1-tmp_vec1/norm(tmp_vec1)*SPO_dis_chl2wall;
            tmp_pt2_a=tmp_pt1_a+tmp_vec2;
            tmp_pt2_b=tmp_pt1_b+tmp_vec2;
            tmp_pt3_a=tmp_pt1_a+tmp_vec3;
            tmp_pt3_b=tmp_pt1_b+tmp_vec3;
            
            tmptag_wp_a=['wp', num2str(count_wp)]; count_wp=count_wp+1;
            model.component('comp1').geom('geom1').create(tmptag_wp_a, 'WorkPlane');
            model.component('comp1').geom('geom1').feature(tmptag_wp_a).set('planetype', 'coordinates');
            model.component('comp1').geom('geom1').feature(tmptag_wp_a).set('genpoints', [tmp_pt1_a; tmp_pt2_a; tmp_pt3_a]);
            model.component('comp1').geom('geom1').feature(tmptag_wp_a).set('unite', true);
            tmptag_wp_b=['wp', num2str(count_wp)]; count_wp=count_wp+1;
            model.component('comp1').geom('geom1').create(tmptag_wp_b, 'WorkPlane');
            model.component('comp1').geom('geom1').feature(tmptag_wp_b).set('planetype', 'coordinates');
            model.component('comp1').geom('geom1').feature(tmptag_wp_b).set('genpoints', [tmp_pt1_b; tmp_pt2_b; tmp_pt3_b]);
            model.component('comp1').geom('geom1').feature(tmptag_wp_b).set('unite', true);
            %partition both obj with workplane
            %target = obj loop_i
            tmptag_obj=list_SPO_MS{loop_i};%{} since it's a cell array
            tmptag_par=['par', num2str(count_par)];count_par=count_par+1;
            model.component('comp1').geom('geom1').create(tmptag_par, 'Partition');
            model.component('comp1').geom('geom1').feature(tmptag_par).set('partitionwith', 'workplane');
            model.component('comp1').geom('geom1').feature(tmptag_par).set('workplane', tmptag_wp_a);
            model.component('comp1').geom('geom1').feature(tmptag_par).selection('input').set(tmptag_obj);
            %split par
            tmptag_spl=['spl', num2str(count_spl)];count_spl=count_spl+1;
            model.component('comp1').geom('geom1').create(tmptag_spl, 'Split');
            model.component('comp1').geom('geom1').feature(tmptag_spl).selection('input').set({tmptag_par});
            model.component('comp1').geom('geom1').run(tmptag_spl);
            if(numel(model.geom('geom1').feature(tmptag_spl).objectNames())==2)
                %select domain, measure vol
                tmptag_spl_dom1=[tmptag_spl,'(1)'];
                model.component('comp1').geom('geom1').measure.selection.init(3);
                model.component('comp1').geom('geom1').measure.selection.set(tmptag_spl_dom1, 1);
                tmp_vol_spl_dom1 = model.geom('geom1').measure().getVolume();
                tmptag_spl_dom2=[tmptag_spl,'(2)'];
                model.component('comp1').geom('geom1').measure.selection.init(3);
                model.component('comp1').geom('geom1').measure.selection.set(tmptag_spl_dom2, 1);
                tmp_vol_spl_dom2 = model.geom('geom1').measure().getVolume();
                %delete the smaller domain, record the remain obj for this MS
                tmptag_del=['del', num2str(count_del)];count_del=count_del+1;
                model.component('comp1').geom('geom1').create(tmptag_del, 'Delete');
                model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').init;
                if tmp_vol_spl_dom1 < tmp_vol_spl_dom2
                    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(cellstr(tmptag_spl_dom1));
                    list_SPO_MS(loop_i)=cellstr(tmptag_spl_dom2);
                else
                    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(cellstr(tmptag_spl_dom2));
                    list_SPO_MS(loop_i)=cellstr(tmptag_spl_dom1);
                end
            else
                list_SPO_MS(loop_i)=cellstr(tmptag_spl);
            end
            
            %target = obj loop_j
            tmptag_obj=list_SPO_MS{loop_j};%{} since it's a cell array
            tmptag_par=['par', num2str(count_par)];count_par=count_par+1;
            model.component('comp1').geom('geom1').create(tmptag_par, 'Partition');
            model.component('comp1').geom('geom1').feature(tmptag_par).set('partitionwith', 'workplane');
            model.component('comp1').geom('geom1').feature(tmptag_par).set('workplane', tmptag_wp_b);
            model.component('comp1').geom('geom1').feature(tmptag_par).selection('input').set(tmptag_obj);
            %split par
            tmptag_spl=['spl', num2str(count_spl)];count_spl=count_spl+1;
            model.component('comp1').geom('geom1').create(tmptag_spl, 'Split');
            model.component('comp1').geom('geom1').feature(tmptag_spl).selection('input').set({tmptag_par});
            model.component('comp1').geom('geom1').run(tmptag_spl);
            if(numel(model.geom('geom1').feature(tmptag_spl).objectNames())==2)
                %select domain, measure vol
                tmptag_spl_dom1=[tmptag_spl,'(1)'];
                model.component('comp1').geom('geom1').measure.selection.init(3);
                model.component('comp1').geom('geom1').measure.selection.set(tmptag_spl_dom1, 1);
                tmp_vol_spl_dom1 = model.geom('geom1').measure().getVolume();
                tmptag_spl_dom2=[tmptag_spl,'(2)'];
                model.component('comp1').geom('geom1').measure.selection.init(3);
                model.component('comp1').geom('geom1').measure.selection.set(tmptag_spl_dom2, 1);
                tmp_vol_spl_dom2 = model.geom('geom1').measure().getVolume();
                %delete the smaller domain, record the remain obj for this MS
                tmptag_del=['del', num2str(count_del)];count_del=count_del+1;
                model.component('comp1').geom('geom1').create(tmptag_del, 'Delete');
                model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').init;
                if tmp_vol_spl_dom1 < tmp_vol_spl_dom2
                    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(tmptag_spl_dom1);
                    list_SPO_MS(loop_j)=cellstr(tmptag_spl_dom2);
                else
                    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(tmptag_spl_dom2);
                    list_SPO_MS(loop_j)=cellstr(tmptag_spl_dom1);
                end
            else
                list_SPO_MS(loop_j)=cellstr(tmptag_spl);
            end
        end
    end
    
%     %contact with box?
%     if mat_contact_cell_idx(loop_i,loop_i)~=0
%         
%     end
end

%%% for SPO CHL;
%%% convert to surface; if contact with box, then intersect with box;
SPO_CHL_tag_set={};
for loop_i=1:Num_spo
    if mat_contact_cell_idx(loop_i,loop_i)~=0 % not contact with box
        tmptag_csur=['csur', num2str(count_csur)];count_csur=count_csur+1;
        model.component('comp1').geom('geom1').create(tmptag_csur, 'ConvertToSurface');
        %model.component('comp1').geom('geom1').feature(tmptag_csur).set('selresult', true);
        model.component('comp1').geom('geom1').feature(tmptag_csur).selection('input').set(list_SPO_MS{loop_i});%[??] tmptag_uni is str; cellstr() is cell;
        %%% intersect csur with box
        tmptag_int=['int', num2str(count_int)];count_int=count_int+1;
        model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
        model.component('comp1').geom('geom1').feature(tmptag_int).set('selresult', true);
        model.component('comp1').geom('geom1').feature(tmptag_int).set('selresultshow', 'bnd');
        model.component('comp1').geom('geom1').feature(tmptag_int).set('keep', true);
        model.component('comp1').geom('geom1').feature(tmptag_int).set('intbnd', false);
        tmpcellstr=[cellstr('blk2'),cellstr(tmptag_csur)];
        model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set(tmpcellstr);%set({'blk1' 'csur1'})
        SPO_CHL_tag_set{loop_i}=tmptag_int;
        %%% delete csur
        tmptag_del=['del', num2str(count_del)];count_del=count_del+1;
        model.component('comp1').geom('geom1').create(tmptag_del, 'Delete');
        model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').init;
        model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(tmptag_csur);
    else
        tmptag_csur=['csur', num2str(count_csur)];count_csur=count_csur+1;
        model.component('comp1').geom('geom1').create(tmptag_csur, 'ConvertToSurface');
        model.component('comp1').geom('geom1').feature(tmptag_csur).set('selresult', true);
        model.component('comp1').geom('geom1').feature(tmptag_csur).selection('input').set(list_SPO_MS{loop_i});
        SPO_CHL_tag_set{loop_i}=tmptag_csur;
    end
end
tmptag_del=['del', num2str(count_del)];count_del=count_del+1;
model.component('comp1').geom('geom1').create(tmptag_del, 'Delete');
model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').init;
model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set({'blk2'});

model.geom('geom1').feature('fin').name('Form Assembly');
model.geom('geom1').feature('fin').set('action', 'assembly');
model.geom('geom1').feature('fin').set('createpairs', 'on');
model.geom('geom1').feature('fin').set('imprint', 'on');
%model.geom('geom1').feature('fin').set('repairtol', '1e-5');
%model.geom('geom1').feature(tmptag_IAS).set('repairtoltype', 'auto');
model.geom('geom1').run;

%NBnd=model.geom('geom1').getNBoundaries();% number of boundaries
%NCreSel=model.selection.size();% number of created selection, i.e. MS surface
TagsCreSel=model.selection.tags();% almost equal to SPO_CHL_tag_set

count_mesh=1;
count_mesh_size=1;
count_mesh_ftri=1;
SPO_CHL_center_trimesh_pts={};%store chloroplasts centroid for all SPO cells
SPO_CHL_center_trimesh_tri={};
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
    model.component('comp1').mesh(tmptag_mesh).feature(tmptag_size).set('hmin', SPO_diameter_chl_patch/50);
    model.component('comp1').mesh(tmptag_mesh).feature(tmptag_size).set('hmaxactive', true);
    model.component('comp1').mesh(tmptag_mesh).feature(tmptag_size).set('hmax', SPO_diameter_chl_patch/20);%max distance allow
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
    
    SPO_CHL_center_trimesh_pts{loop_i}=pts;
    SPO_CHL_center_trimesh_tri{loop_i}=tri;
end

% % check distance matrix
% for loop_i=1:numel(TagsCreSel)
%     tmp_pts=SPO_CHL_center_trimesh_pts{loop_i};
%     tmp_tri=SPO_CHL_center_trimesh_tri{loop_i};
%     tmpd1=norm(tmp_pts(tmp_tri(loop_i,1),:)-tmp_pts(tmp_tri(loop_i,2),:));
%     tmpd2=norm(tmp_pts(tmp_tri(loop_i,2),:)-tmp_pts(tmp_tri(loop_i,3),:));
%     tmpd3=norm(tmp_pts(tmp_tri(loop_i,3),:)-tmp_pts(tmp_tri(loop_i,1),:));
%     tmpd_all(loop_i,:)=[tmpd1,tmpd2,tmpd3];
% end

% generate centroid from {pts,tri}
SPO_CHL_center_final_pts={};
for loop_i=1:numel(TagsCreSel)
    tmp_pts=SPO_CHL_center_trimesh_pts{loop_i};
    tmp_tri=SPO_CHL_center_trimesh_tri{loop_i};
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
                if any(tmp_dis_sq_all<=SPO_diameter_chl_patch^2)
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
    SPO_CHL_center_final_pts{loop_i}=tmp_sel_pts;
end

% plot sphere; check maximum coverage?
for loop_i=1:numel(TagsCreSel)
    tmp_pts=SPO_CHL_center_final_pts{loop_i};
    %create CHL sphere
    tmptag_set={};
    for loop_j=1:size(tmp_pts,1)
        tmptag_sph=['sph',num2str(count_sph)];count_sph=count_sph+1;
        model.component('comp1').geom('geom1').create(tmptag_sph, 'Sphere');
        model.component('comp1').geom('geom1').feature(tmptag_sph).set('r', 'SPO_radius_chl_sphere');
        model.component('comp1').geom('geom1').feature(tmptag_sph).set('pos', string(tmp_pts(loop_j,:)));%{'5e-6' '0' '150e-6'}
        tmptag_set{loop_j}=tmptag_sph;
    end
    model.component('comp1').geom('geom1').run(tmptag_sph);
    %merge CHL sphere
    tmptag_uni=['uni', num2str(count_uni)];count_uni=count_uni+1;
    model.component('comp1').geom('geom1').create(tmptag_uni, 'Union');
    model.component('comp1').geom('geom1').feature(tmptag_uni).set('intbnd', false);
    model.component('comp1').geom('geom1').feature(tmptag_uni).selection('input').set(tmptag_set);
    %calculate surface area
    model.component('comp1').geom('geom1').measure.selection.init(2);
    tmp_NFaces=model.geom('geom1').obj(SPO_CHL_tag_set{loop_i}).getNFaces;
    model.component('comp1').geom('geom1').measure.selection.set(SPO_CHL_tag_set{loop_i},1:tmp_NFaces);%set('int1',[1,2,3]);%[1,2,3] is corresponding boundary
    full_area(loop_i)=model.geom('geom1').measure().getVolume();%Remark: getVolume here actually getArea
    %int operation
    tmptag_int=['int',num2str(count_int)];count_int=count_int+1;
    model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
    model.component('comp1').geom('geom1').feature(tmptag_int).set('keep', false);
    model.component('comp1').geom('geom1').feature(tmptag_int).set('intbnd', false);
    tmpcellstr=[SPO_CHL_tag_set(loop_i),cellstr(tmptag_uni)];
    model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set(tmpcellstr);
    %calculate projected area
    model.component('comp1').geom('geom1').run(tmptag_int);
    model.component('comp1').geom('geom1').measure.selection.init(2);
    tmp_NFaces=model.geom('geom1').obj(tmptag_int).getNFaces;%%alternative method: meshstats - getNEntities()
    model.component('comp1').geom('geom1').measure.selection.set(tmptag_int,1:tmp_NFaces);%set('int1',[1,2,3]);%[1,2,3] is corresponding boundary
    project_area(loop_i)=model.geom('geom1').measure().getVolume();
    coverage(loop_i)=project_area(loop_i)/full_area(loop_i);
end

save tmp_SPO_CHL_distribute.mat -regexp '^(?!(model|ans)$).'
disp('Potential centroids for chloroplasts in Spongy are ready.')
mphsave(model, 'tmp_SPO_CHL_geo.mph')
ModelUtil.remove('Model_SPO_CHL_distribute');
