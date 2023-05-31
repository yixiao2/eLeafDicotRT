% eLeaf: 3D model of dicot leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 0.1

%% calculate volume and surface area information
%% post-processing after RT

% 1. volume of PAL, SPO
load tmp_PAL_MS_distribute.mat tmp_vol_IAS_PAL_new PAL_box_thick PAL_box_length PAL_box_length
INFO_vol_pal_msc=PAL_box_thick*PAL_box_length*PAL_box_length*(1-tmp_vol_IAS_PAL_new);
load tmp_SPO_MS_distribute.mat tmp_vol_IAS_SPO SPO_box_thick SPO_box_length SPO_box_length
INFO_vol_spo_msc=SPO_box_thick*SPO_box_length*SPO_box_length*(1-tmp_vol_IAS_SPO);

% 2. volume of chloroplasts in PAL and SPO
load SAVE_e_geom4RT.mat chlo_vol_profile GLB_count_PAL_chl GLB_count_SPO_chl
numofchlo_allpal=sum(GLB_count_PAL_chl);
numofchlo_allspo=sum(GLB_count_SPO_chl);
INFO_vol_pal_chlo=sum(chlo_vol_profile(1:numofchlo_allpal));
INFO_vol_spo_chlo=sum(chlo_vol_profile((numofchlo_allpal+1):end));

% 3. surface area of PAL, SPO
% continue from e_geo_PAL_MS_distribute_v0_1.m
% remark: exclude box boundaries
import com.comsol.model.*
import com.comsol.model.util.*
ModelUtil.remove('Model_PAL_MS_distribute')
load('tmp_PAL_MS_distribute.mat')
model=mphopen('tmp_PAL_MS_geo.mph','Model_PAL_MS_distribute');
model.component('comp1').geom('geom1').feature.remove('dif1');
num_pal=numel(PAL_tag_set);
count_csur=1;
count_int=1;
count_del=1;
PAL_csur_tag_set={};
for loop_cell=1:num_pal
    %%% convert to surface
    tmptag_csur=['csur', num2str(count_csur)];count_csur=count_csur+1;
    model.component('comp1').geom('geom1').create(tmptag_csur, 'ConvertToSurface');
    model.component('comp1').geom('geom1').feature(tmptag_csur).selection('input').set(PAL_tag_set(loop_cell));
    %%% intesect with box
    tmptag_int=['int', num2str(count_int)];count_int=count_int+1;
    model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
    model.component('comp1').geom('geom1').feature(tmptag_int).set('keep', true);
    model.component('comp1').geom('geom1').feature(tmptag_int).set('intbnd', false);
    tmpcellstr=[cellstr('blk1'),cellstr(tmptag_csur)];
    model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set(tmpcellstr);
    %%% delete csur
    tmptag_del=['del', num2str(count_del)];count_del=count_del+1;
    model.component('comp1').geom('geom1').create(tmptag_del, 'Delete');
    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').init;
    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(tmptag_csur);
    PAL_csur_tag_set{loop_cell}=tmptag_int;
end
%%% delete box
tmptag_del=['del', num2str(count_del)];count_del=count_del+1;
model.component('comp1').geom('geom1').create(tmptag_del, 'Delete');
model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').init;
model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set({'blk1'});
%%% merge and measure surface area
tmptag_uni=['uni', num2str(count_uni)];count_uni=count_uni+1;
model.component('comp1').geom('geom1').create(tmptag_uni, 'Union');
model.component('comp1').geom('geom1').feature(tmptag_uni).set('intbnd', false);
model.component('comp1').geom('geom1').feature(tmptag_uni).selection('input').set(PAL_csur_tag_set);

model.component('comp1').geom('geom1').run(tmptag_uni);
model.component('comp1').geom('geom1').measure.selection.init(2);
tmp_NFaces=model.geom('geom1').obj(tmptag_int).getNFaces;%%alternative method: meshstats - getNEntities()
model.component('comp1').geom('geom1').measure.selection.set(tmptag_uni,1:tmp_NFaces);%set('int1',[1,2,3]);%[1,2,3] is corresponding boundary
INFO_sur_pal_msc=model.geom('geom1').measure().getVolume();

mphsave(model, 'tmp_PAL_MS_geo_forSm.mph')
%ModelUtil.remove('Model_PAL_MS_distribute')%% [?] crash matlab

%%% 3-SPO MSC
ModelUtil.remove('Model_SPO_MS_distribute')
load tmp_SPO_MS_distribute.mat
model=mphopen('tmp_SPO_MS_geo.mph','Model_SPO_MS_distribute');

model.component('comp1').geom('geom1').feature.remove('dif1');

%%%create extended blk1, i.e. blk2 for retaining surface areas intersect with upper and
%%%lower boundaries of box
model.component('comp1').geom('geom1').create('blk2', 'Block');
model.component('comp1').geom('geom1').feature('blk2').set('pos', [-SPO_MS_radius*2 -SPO_MS_radius*2 0]);
model.component('comp1').geom('geom1').feature('blk2').set('size', {'SPO_box_length+SPO_MS_radius*4' 'SPO_box_length+SPO_MS_radius*4' 'SPO_box_thick'});

%%%process spheric MS to create MS with contact area, based on center_set and radius_set
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
                GLB_mat_contact_SPO_MSVAC_idx(loop_i)=1;
            else
                GLB_mat_contact_SPO_MSVAC_idx(loop_i)=0;
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

%%% cut some spheres with blk2 first, then convert to surfaces, then cut
%%% with blk1
num_spo=numel(SPO_tag_set);
count_csur=1;
count_int=1;
count_del=1;
count_wp=1;
count_par=1;
count_spl=1;
SPO_csur_tag_set={};
GLB_tag_SPO_MS=SPO_tag_set;
for loop_i=1:num_spo
    %%% IF - for spheres intersect with upper or lower, cut with blk2
    %%% then, for spheres intersect with any box boundaries, convert to
    %%% surface and cut with blk1
    if SPO_MS_center_set_Z(loop_i)<SPO_MS_radius||SPO_MS_center_set_Z(loop_i)>(SPO_box_thick-SPO_MS_radius)
        tmptag_int=['int', num2str(count_int)];count_int=count_int+1;
        model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
        model.component('comp1').geom('geom1').feature(tmptag_int).set('keep', true);
        model.component('comp1').geom('geom1').feature(tmptag_int).set('intbnd', false);
        tmpcellstr=[cellstr('blk2'),cellstr(GLB_tag_SPO_MS{loop_i})];
        model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set(tmpcellstr);

        tmptag_del=['del', num2str(count_del)];count_del=count_del+1;
        model.component('comp1').geom('geom1').create(tmptag_del, 'Delete');
        model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').init;
        model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(GLB_tag_SPO_MS{loop_i});
        GLB_tag_SPO_MS{loop_i}=tmptag_int;
    end

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

            %%% (1/4)process MS shape
            tmptag_wp=['wp', num2str(count_wp)]; count_wp=count_wp+1;
            model.component('comp1').geom('geom1').create(tmptag_wp, 'WorkPlane');
            model.component('comp1').geom('geom1').feature(tmptag_wp).set('planetype', 'coordinates');
            model.component('comp1').geom('geom1').feature(tmptag_wp).set('genpoints', [tmp_pt1; tmp_pt2; tmp_pt3]);
            model.component('comp1').geom('geom1').feature(tmptag_wp).set('unite', true);
            %partition both obj with workplane
            %target = obj loop_i
            tmptag_obj=GLB_tag_SPO_MS{loop_i};%{} since it's a cell array
            tmptag_par=['par', num2str(count_par)];count_par=count_par+1;
            model.component('comp1').geom('geom1').create(tmptag_par, 'Partition');
            model.component('comp1').geom('geom1').feature(tmptag_par).set('partitionwith', 'workplane');
            model.component('comp1').geom('geom1').feature(tmptag_par).set('workplane', tmptag_wp);
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
                    GLB_tag_SPO_MS(loop_i)=cellstr(tmptag_spl_dom2);
                else
                    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(tmptag_spl_dom2);
                    GLB_tag_SPO_MS(loop_i)=cellstr(tmptag_spl_dom1);
                end
            else
                GLB_tag_SPO_MS(loop_i)=cellstr(tmptag_spl);
            end
            %target = obj loop_j
            tmptag_obj=GLB_tag_SPO_MS{loop_j};%{} since it's a cell array
            tmptag_par=['par', num2str(count_par)];count_par=count_par+1;
            model.component('comp1').geom('geom1').create(tmptag_par, 'Partition');
            model.component('comp1').geom('geom1').feature(tmptag_par).set('partitionwith', 'workplane');
            model.component('comp1').geom('geom1').feature(tmptag_par).set('workplane', tmptag_wp);
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
                    GLB_tag_SPO_MS(loop_j)=cellstr(tmptag_spl_dom2);
                else
                    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(tmptag_spl_dom2);
                    GLB_tag_SPO_MS(loop_j)=cellstr(tmptag_spl_dom1);
                end
            else
                GLB_tag_SPO_MS(loop_j)=cellstr(tmptag_spl);
            end
        end
    end
end

for loop_cell=1:num_spo
    %%% convert to surface
    tmptag_csur=['csur', num2str(count_csur)];count_csur=count_csur+1;
    model.component('comp1').geom('geom1').create(tmptag_csur, 'ConvertToSurface');
    model.component('comp1').geom('geom1').feature(tmptag_csur).selection('input').set(GLB_tag_SPO_MS(loop_cell));
    %%% intesect with box
    tmptag_int=['int', num2str(count_int)];count_int=count_int+1;
    model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
    model.component('comp1').geom('geom1').feature(tmptag_int).set('keep', true);
    model.component('comp1').geom('geom1').feature(tmptag_int).set('intbnd', false);
    tmpcellstr=[cellstr('blk1'),cellstr(tmptag_csur)];
    model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set(tmpcellstr);
    %%% delete csur
    tmptag_del=['del', num2str(count_del)];count_del=count_del+1;
    model.component('comp1').geom('geom1').create(tmptag_del, 'Delete');
    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').init;
    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(tmptag_csur);
    SPO_csur_tag_set{loop_cell}=tmptag_int;

    model.component('comp1').geom('geom1').run(tmptag_del);
    model.component('comp1').geom('geom1').measure.selection.init(2);
    tmp_NFaces=model.geom('geom1').obj(tmptag_int).getNFaces;%%alternative method: meshstats - getNEntities()
    model.component('comp1').geom('geom1').measure.selection.set(tmptag_int,1:tmp_NFaces);%set('int1',[1,2,3]);%[1,2,3] is corresponding boundary
    sur_spo_msc_vec(loop_cell)=model.geom('geom1').measure().getVolume();
end
%%% delete box
tmptag_del=['del', num2str(count_del)];count_del=count_del+1;
model.component('comp1').geom('geom1').create(tmptag_del, 'Delete');
model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').init;
model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set({'blk1','blk2'});
model.component('comp1').geom('geom1').run(tmptag_del);
INFO_sur_spo_msc=sum(sur_spo_msc_vec);
mphsave(model, 'tmp_SPO_MS_geo_forSm.mph')


% 4. surface area of SPO in contact with IAS
ModelUtil.remove('Model_SPO_MS_distribute')
load tmp_SPO_MS_distribute.mat
model=mphopen('tmp_SPO_MS_geo.mph','Model_SPO_MS_distribute');

model.component('comp1').geom('geom1').feature.remove('dif1');

%%%create extended blk1, i.e. blk2 for retaining surface areas intersect with upper and
%%%lower boundaries of box
model.component('comp1').geom('geom1').create('blk2', 'Block');
model.component('comp1').geom('geom1').feature('blk2').set('pos', [-SPO_MS_radius*2 -SPO_MS_radius*2 0]);
model.component('comp1').geom('geom1').feature('blk2').set('size', {'SPO_box_length+SPO_MS_radius*4' 'SPO_box_length+SPO_MS_radius*4' 'SPO_box_thick'});

%%%mat_contact_cell_idx is already calculated above

%%% cut some spheres with blk2 first, then convert to surfaces, then cut
%%% with blk1
num_spo=numel(SPO_tag_set);
count_csur=1;
count_int=1;
count_del=1;
count_wp=1;
count_par=1;
count_spl=1;
count_cyl=1;
SPO_csur_tag_set={};
GLB_tag_SPO_MS=SPO_tag_set;

for loop_i=1:num_spo
    %%% IF - for spheres intersect with upper or lower, cut with blk2
    %%% then, for spheres intersect with any box boundaries, convert to
    %%% surface and cut with blk1
    if SPO_MS_center_set_Z(loop_i)<SPO_MS_radius||SPO_MS_center_set_Z(loop_i)>(SPO_box_thick-SPO_MS_radius)
        tmptag_int=['int', num2str(count_int)];count_int=count_int+1;
        model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
        model.component('comp1').geom('geom1').feature(tmptag_int).set('keep', true);
        model.component('comp1').geom('geom1').feature(tmptag_int).set('intbnd', false);
        tmpcellstr=[cellstr('blk2'),cellstr(GLB_tag_SPO_MS{loop_i})];
        model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set(tmpcellstr);

        tmptag_del=['del', num2str(count_del)];count_del=count_del+1;
        model.component('comp1').geom('geom1').create(tmptag_del, 'Delete');
        model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').init;
        model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(GLB_tag_SPO_MS{loop_i});
        GLB_tag_SPO_MS{loop_i}=tmptag_int;
    end

    %convert sph to surface; prepare for later cutting with box
    tmptag_csur=['csur', num2str(count_csur)];count_csur=count_csur+1;
    model.component('comp1').geom('geom1').create(tmptag_csur, 'ConvertToSurface');
    model.component('comp1').geom('geom1').feature(tmptag_csur).selection('input').set(GLB_tag_SPO_MS{loop_i});

    SPO_csur_tag_set{loop_i}=tmptag_csur;
end

for loop_i=1:num_spo
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

            %%%use cylinder to cut surface, instead using plane to split sphere.
            %%%cylinder radius SPO_MS_radius*2; cylinder height
            %%%SPO_MS_radius*3
            tmptag_cyl_a=['cyl', num2str(count_cyl)]; count_cyl=count_cyl+1;
            model.component('comp1').geom('geom1').create(tmptag_cyl_a, 'Cylinder');
            model.component('comp1').geom('geom1').feature(tmptag_cyl_a).set('r', SPO_MS_radius*2);
            model.component('comp1').geom('geom1').feature(tmptag_cyl_a).set('h', SPO_MS_radius*3);
            model.component('comp1').geom('geom1').feature(tmptag_cyl_a).set('pos', tmp_pt1);
            model.component('comp1').geom('geom1').feature(tmptag_cyl_a).set('axistype', 'cartesian');
            model.component('comp1').geom('geom1').feature(tmptag_cyl_a).set('axis', tmp_vec1);
            tmptag_cyl_b=['cyl', num2str(count_cyl)]; count_cyl=count_cyl+1;
            model.component('comp1').geom('geom1').create(tmptag_cyl_b, 'Cylinder');
            model.component('comp1').geom('geom1').feature(tmptag_cyl_b).set('r', SPO_MS_radius*2);
            model.component('comp1').geom('geom1').feature(tmptag_cyl_b).set('h', SPO_MS_radius*3);
            model.component('comp1').geom('geom1').feature(tmptag_cyl_b).set('pos', tmp_pt1);
            model.component('comp1').geom('geom1').feature(tmptag_cyl_b).set('axistype', 'cartesian');
            model.component('comp1').geom('geom1').feature(tmptag_cyl_b).set('axis', -tmp_vec1);
            
            %cut surface with cyl correspondingly
            %target = obj loop_i
            tmptag_obj=SPO_csur_tag_set{loop_i};%{} since it's a cell array
            tmptag_int=['int', num2str(count_int)];count_int=count_int+1;
            model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
            model.component('comp1').geom('geom1').feature(tmptag_int).set('keep', false);
            model.component('comp1').geom('geom1').feature(tmptag_int).set('intbnd', false);
            tmpcellstr=[cellstr(tmptag_cyl_a),cellstr(tmptag_obj)];
            model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set(tmpcellstr);%set({'blk1' 'csur1'})
            SPO_csur_tag_set{loop_i}=tmptag_int;
            
            %target = obj loop_j
            tmptag_obj=SPO_csur_tag_set{loop_j};%{} since it's a cell array
            tmptag_int=['int', num2str(count_int)];count_int=count_int+1;
            model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
            model.component('comp1').geom('geom1').feature(tmptag_int).set('keep', false);
            model.component('comp1').geom('geom1').feature(tmptag_int).set('intbnd', false);
            tmpcellstr=[cellstr(tmptag_cyl_b),cellstr(tmptag_obj)];
            model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set(tmpcellstr);%set({'blk1' 'csur1'})
            SPO_csur_tag_set{loop_j}=tmptag_int;
            
        end
    end
end

for loop_cell=1:num_spo
%     %%% convert to surface
%     tmptag_csur=['csur', num2str(count_csur)];count_csur=count_csur+1;
%     model.component('comp1').geom('geom1').create(tmptag_csur, 'ConvertToSurface');
%     model.component('comp1').geom('geom1').feature(tmptag_csur).selection('input').set(GLB_tag_SPO_MS(loop_cell));
    %%% intesect with box
    tmptag_int=['int', num2str(count_int)];count_int=count_int+1;
    model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
    model.component('comp1').geom('geom1').feature(tmptag_int).set('keep', true);
    model.component('comp1').geom('geom1').feature(tmptag_int).set('intbnd', false);
    tmpcellstr=[cellstr('blk1'),cellstr(SPO_csur_tag_set{loop_cell})];
    model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set(tmpcellstr);
    %%% delete csur
    tmptag_del=['del', num2str(count_del)];count_del=count_del+1;
    model.component('comp1').geom('geom1').create(tmptag_del, 'Delete');
    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').init;
    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(SPO_csur_tag_set(loop_cell));
    SPO_csur_tag_set{loop_cell}=tmptag_int;

    model.component('comp1').geom('geom1').run(tmptag_del);
    model.component('comp1').geom('geom1').measure.selection.init(2);
    tmp_NFaces=model.geom('geom1').obj(tmptag_int).getNFaces;%%alternative method: meshstats - getNEntities()
    model.component('comp1').geom('geom1').measure.selection.set(tmptag_int,1:tmp_NFaces);%set('int1',[1,2,3]);%[1,2,3] is corresponding boundary
    sur_spo_msc_withias_vec(loop_cell)=model.geom('geom1').measure().getVolume();
end
%%% delete box
tmptag_del=['del', num2str(count_del)];count_del=count_del+1;
model.component('comp1').geom('geom1').create(tmptag_del, 'Delete');
model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').init;
model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set({'blk1','blk2'});
model.component('comp1').geom('geom1').run(tmptag_del);
INFO_sur_spo_msc_withias=sum(sur_spo_msc_withias_vec);
mphsave(model, 'tmp_SPO_MS_geo_forSmwithIAS.mph')

save SAVE_info_surfvol.mat -regexp '^(?!(model|ans)$).'
