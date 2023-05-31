% eLeaf: 3D model of dicot leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 0.1

%% PART ONE %%
% e_geo_SPO_MS_distribute_v0_1.m
% distribute SPO in SPO box to match porosity
%
% PART TWO
% generate centroid for chloroplast
%
% PART THREE
% e_geo_SPO_MT_distribute_v0_1.m
% generate centroid for mitochondria

%%%[To update]
%%%adjust some SPO_MS radius to match MS_porosity exactly
%%%maybe shrink 0.9^3 at most; iterate algorithm

clear all;
load parainput.mat

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model_SPO_MS_distribute');

model.param.set('SPO_box_thick', [num2str(SPO_box_thick),'[m]']);
model.param.set('SPO_box_length', [num2str(SPO_box_length),'[m]']);
model.param.set('SPO_MS_radius', [num2str(SPO_MS_radius),'[m]']);

model.component.create('comp1', true);
model.component('comp1').geom.create('geom1', 3);

model.component('comp1').geom('geom1').create('blk1', 'Block');
model.component('comp1').geom('geom1').feature('blk1').set('pos', [0 0 0]);
model.component('comp1').geom('geom1').feature('blk1').set('size', {'SPO_box_length' 'SPO_box_length' 'SPO_box_thick'});

%%% distribute center of spheric spongy evenly and randomly
%dis_SPO_MS_center = 13e-6;% minimum distance between two centers
tmp_dis_2=dis_SPO_MS_center^2;
Num_spo=1;% increase number of chloroplast when tmp_IAS<MS_porosity
count_sph=1;
count_dif=1;
% random x from [dis_SPO_MS_center, SPO_box_length-dis_SPO_MS_center]
% random y from [dis_SPO_MS_center, SPO_box_length-dis_SPO_MS_center]
% random z from [dis_SPO_MS_center, SPO_box_thick-dis_SPO_MS_center]
% SPO_MS_center_set_X=dis_SPO_MS_center+(SPO_box_length-dis_SPO_MS_center*2)*rand;
% SPO_MS_center_set_Y=dis_SPO_MS_center+(SPO_box_length-dis_SPO_MS_center*2)*rand;
% SPO_MS_center_set_Z=dis_SPO_MS_center+(SPO_box_thick-dis_SPO_MS_center*2)*rand;
SPO_MS_center_set_X=SPO_box_length*rand;
SPO_MS_center_set_Y=SPO_box_length*rand;
SPO_MS_center_set_Z=SPO_MS_radius/2+(SPO_box_thick-SPO_MS_radius/2*2)*rand;
tmptag_spo=['sph',num2str(count_sph)];
model.component('comp1').geom('geom1').create(tmptag_spo, 'Sphere');
SPO_tag_set{1}=tmptag_spo;
count_sph=count_sph+1;
model.component('comp1').geom('geom1').feature(tmptag_spo).set('r', 'SPO_MS_radius');
tmpstr_pos={num2str(SPO_MS_center_set_X), num2str(SPO_MS_center_set_Y), num2str(SPO_MS_center_set_Z)};
model.component('comp1').geom('geom1').feature(tmptag_spo).set('pos', tmpstr_pos);
model.component('comp1').geom('geom1').run(tmptag_spo);
tmp_vol_IAS_SPO=1;
count_list_SPO_MS=1;
list_SPO_MS(count_list_SPO_MS)=cellstr(tmptag_spo);%list_SPO_MS record obj name of SPO MS
%%[to update] SPO_tag_set is the same to list_SPO_MS? -> yes

numloop=1;
while tmp_vol_IAS_SPO > MS_porosity && numloop < 5e6
    %tmp_center_newX=dis_SPO_MS_center+(SPO_box_length-dis_SPO_MS_center*2)*rand;
    %tmp_center_newY=dis_SPO_MS_center+(SPO_box_length-dis_SPO_MS_center*2)*rand;
    %tmp_center_newZ=dis_SPO_MS_center+(SPO_box_thick-dis_SPO_MS_center*2)*rand;
    tmp_center_newX=SPO_box_length*rand;
    tmp_center_newY=SPO_box_length*rand;
    tmp_center_newZ=SPO_MS_radius/2+(SPO_box_thick-SPO_MS_radius/2*2)*rand;
    if all(((SPO_MS_center_set_X(1:Num_spo)-tmp_center_newX).^2+(SPO_MS_center_set_Y(1:Num_spo)-tmp_center_newY).^2+(SPO_MS_center_set_Z(1:Num_spo)-tmp_center_newZ).^2) > tmp_dis_2)
        Num_spo=Num_spo+1;
        SPO_MS_center_set_X(Num_spo)=tmp_center_newX;
        SPO_MS_center_set_Y(Num_spo)=tmp_center_newY;
        SPO_MS_center_set_Z(Num_spo)=tmp_center_newZ;
        %%build sphere
        tmptag_spo=['sph',num2str(count_sph)];
        model.component('comp1').geom('geom1').create(tmptag_spo, 'Sphere');
        SPO_tag_set{Num_spo}=tmptag_spo;
        count_sph=count_sph+1;
        model.component('comp1').geom('geom1').feature(tmptag_spo).set('r', 'SPO_MS_radius');
        tmpstr_pos={num2str(tmp_center_newX), num2str(tmp_center_newY), num2str(tmp_center_newZ)};
        model.component('comp1').geom('geom1').feature(tmptag_spo).set('pos', tmpstr_pos);
        %%generate IAS_SPO
        tmptag_spo=['dif',num2str(count_dif)];
        if Num_spo~=2
            model.component('comp1').geom('geom1').feature.remove(tmptag_spo);
        end
        model.component('comp1').geom('geom1').create(tmptag_spo, 'Difference');
        %count_dif=count_dif+1;
        model.component('comp1').geom('geom1').feature(tmptag_spo).selection('input').set({'blk1'});
        model.component('comp1').geom('geom1').feature(tmptag_spo).selection('input2').set(SPO_tag_set);
        model.component('comp1').geom('geom1').feature(tmptag_spo).set('intbnd', false);
        %%calculate tmp_vol_IAS_SPO
        model.component('comp1').geom('geom1').run(tmptag_spo);
        model.component('comp1').geom('geom1').measure.selection.init(3);
        model.component('comp1').geom('geom1').measure.selection.set(tmptag_spo, 1);
        tmp_vol_IAS_SPO = model.geom('geom1').measure().getVolume()/(SPO_box_thick*SPO_box_length*SPO_box_length);
        if FLAG_debug_msg==1
            tmp_vol_IAS_SPO
        end
        count_list_SPO_MS=count_list_SPO_MS+1;
        list_SPO_MS(count_list_SPO_MS)=cellstr(tmptag_spo);%list_SPO_MS record obj name of SPO MS
    end
    numloop=numloop+1;
end

save tmp_SPO_MS_distribute.mat -regexp '^(?!(model|ans)$).'
disp('Spongy IAS porosity matched')
mphsave(model, 'tmp_SPO_MS_geo.mph')
ModelUtil.remove('Model_SPO_MS_distribute');
