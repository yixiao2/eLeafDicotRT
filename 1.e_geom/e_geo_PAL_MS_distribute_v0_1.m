% eLeaf: 3D model of dicot leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 0.1

%% PART ONE %%
% e_geo_PAL_MS_distribute_v0_1.m
% distribute PAL in PAL box to match porosity
%
% PART TWO 
% e_geo_PAL_CH_distribute_v0_1.m
% generate centroid for chloroplast
%
% PART THREE
% e_geo_PAL_MT_distribute_v0_1.m
% generate centroid for mitochondria

%%%[to update]
%%%compact distribution of PAL and allow contact area?

clear all;
load parainput.mat

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model_PAL_MS_distribute');

model.param.set('PAL_box_thick', [num2str(PAL_box_thick),'[m]']);
model.param.set('PAL_box_length', [num2str(PAL_box_length),'[m]']);
model.param.set('SPO_box_thick', [num2str(SPO_box_thick),'[m]']);
model.param.set('MODEL_ddis', [num2str(MODEL_ddis),'[m]']);
model.param.set('PAL_MS_radius', [num2str(PAL_MS_radius),'[m]']);
model.param.set('PAL_MS_lcap_flat', [num2str(PAL_MS_lcap_flat),'[m]']);
model.param.set('PAL_MS_ucap_flat', [num2str(PAL_MS_ucap_flat),'[m]']);
model.param.set('PAL_MS_height', 'PAL_box_thick-MODEL_ddis*2-(PAL_MS_radius)*(1-PAL_MS_ucap_flat)-(PAL_MS_radius)*(1-PAL_MS_lcap_flat)');

model.component.create('comp1', true);
model.component('comp1').geom.create('geom1', 3);

model.component('comp1').geom('geom1').create('blk1', 'Block');
model.component('comp1').geom('geom1').feature('blk1').set('pos', {'0' '0' 'SPO_box_thick'});
model.component('comp1').geom('geom1').feature('blk1').set('size', {'PAL_box_length' 'PAL_box_length' 'PAL_box_thick'});

%% PART ONE
%%% distribute PAL and then adjust height to match porosity
count_cyl=1;
count_elp=1;
count_uni=1;
count_PAL=1;
for loop_i=1:N_PAL_rows
    for loop_j=1:N_PAL_cols
        tmptag_pal_cp1=['cyl', num2str(count_cyl)];
        model.component('comp1').geom('geom1').create(tmptag_pal_cp1, 'Cylinder');
        count_cyl=count_cyl+1;
        tmpstr_displx=['(PAL_MS_radius*2+MODEL_ddis)*',num2str(loop_j-1)];
        tmpstr_disply=['(PAL_MS_radius*2+MODEL_ddis)*',num2str(loop_i-1)];
        eval(['model.component(''comp1'').geom(''geom1'').feature(tmptag_pal_cp1).set(''pos'', {''',tmpstr_displx,''' ''',tmpstr_disply,''' ''SPO_box_thick+MODEL_ddis+(PAL_MS_radius)*(1-PAL_MS_lcap_flat)''});'])
        model.component('comp1').geom('geom1').feature(tmptag_pal_cp1).set('r', 'PAL_MS_radius');
        model.component('comp1').geom('geom1').feature(tmptag_pal_cp1).set('h', 'PAL_MS_height');
        
        tmptag_pal_cp2=['elp', num2str(count_elp)];
        model.component('comp1').geom('geom1').create(tmptag_pal_cp2, 'Ellipsoid');
        count_elp=count_elp+1;
        eval(['model.component(''comp1'').geom(''geom1'').feature(tmptag_pal_cp2).set(''pos'', {''',tmpstr_displx,''' ''',tmpstr_disply,''' ''SPO_box_thick+MODEL_ddis+(PAL_MS_radius)*(1-PAL_MS_lcap_flat)''});'])
        model.component('comp1').geom('geom1').feature(tmptag_pal_cp2).set('semiaxes', {'PAL_MS_radius' 'PAL_MS_radius' '(PAL_MS_radius)*(1-PAL_MS_lcap_flat)'});
        
        tmptag_pal_cp3=['elp', num2str(count_elp)];
        count_elp=count_elp+1;
        model.component('comp1').geom('geom1').create(tmptag_pal_cp3, 'Ellipsoid');
        eval(['model.component(''comp1'').geom(''geom1'').feature(tmptag_pal_cp3).set(''pos'', {''',tmpstr_displx,''' ''',tmpstr_disply,''' ''SPO_box_thick+MODEL_ddis+(PAL_MS_radius)*(1-PAL_MS_lcap_flat)+PAL_MS_height''});'])
        model.component('comp1').geom('geom1').feature(tmptag_pal_cp3).set('semiaxes', {'PAL_MS_radius' 'PAL_MS_radius' '(PAL_MS_radius)*(1-PAL_MS_ucap_flat)'});
        
        tmptag_uni=['uni', num2str(count_uni)];
        model.component('comp1').geom('geom1').create(tmptag_uni, 'Union');
        count_uni=count_uni+1;
        model.component('comp1').geom('geom1').feature(tmptag_uni).set('intbnd', false);
        tmptag_set={tmptag_pal_cp1,tmptag_pal_cp2,tmptag_pal_cp3};
        model.component('comp1').geom('geom1').feature(tmptag_uni).selection('input').set(tmptag_set);
        
        PAL_tag_set{count_PAL}=tmptag_uni;
        count_PAL=count_PAL+1;
    end
end

model.component('comp1').geom('geom1').create('dif1', 'Difference');
model.component('comp1').geom('geom1').feature('dif1').set('intbnd', false);
model.component('comp1').geom('geom1').feature('dif1').selection('input').set({'blk1'});
model.component('comp1').geom('geom1').feature('dif1').selection('input2').set(PAL_tag_set);

model.component('comp1').geom('geom1').run('dif1');%% IMPORTANT: run before measure
model.component('comp1').geom('geom1').measure.selection.init(3);
model.component('comp1').geom('geom1').measure.selection.set('dif1', 1);
tmp_vol_IAS_PAL = model.geom('geom1').measure().getVolume()/(PAL_box_thick*PAL_box_length*PAL_box_length)
%MS_porosity = 0.328;
% if tmp_vol_IAS_PAL is not close to MS_porosity
if tmp_vol_IAS_PAL>MS_porosity
    error('Unable to reach that porosity in PAL MS.')
else
    delta_cylinder_d_all = (PAL_box_thick*PAL_box_length*PAL_box_length)*(MS_porosity-tmp_vol_IAS_PAL)/(pi*(PAL_MS_radius)^2/4);
    % divide delta_cylinder_d to all PAL
    N_rand=2+2*(N_PAL_cols-2);
    tmp_rand = rand(N_rand,N_rand);%1 unit = 1/4 of a cylinder; take averge of 2 or 4 units for 1/2 or 1 cylinder respectively.
    %% to improve; use matrix operation
    tmp_delta_d_mat = delta_cylinder_d_all*tmp_rand/sum(tmp_rand,'All');
    %1/4 cylinder on the corners, i.e. 1 unit
    delta_cylinder_d_mat(1,1)=tmp_delta_d_mat(1,1);
    delta_cylinder_d_mat(1,N_PAL_cols)=tmp_delta_d_mat(1,N_rand);
    delta_cylinder_d_mat(N_PAL_cols,1)=tmp_delta_d_mat(N_rand,1);
    delta_cylinder_d_mat(N_PAL_cols,N_PAL_cols)=tmp_delta_d_mat(N_rand,N_rand);
    if N_PAL_cols >=3
        %1/2 cylinder on the boundary edge, i.e. 2 units
        delta_cylinder_d_mat(1,2:(N_PAL_cols-1))=(tmp_delta_d_mat(1,2:2:(N_rand-1))+tmp_delta_d_mat(1,3:2:(N_rand-1)))/2;
        delta_cylinder_d_mat(N_PAL_cols,2:(N_PAL_cols-1))=(tmp_delta_d_mat(N_rand,2:2:(N_rand-1))+tmp_delta_d_mat(N_rand,3:2:(N_rand-1)))/2;
        delta_cylinder_d_mat(2:(N_PAL_cols-1),1)=(tmp_delta_d_mat(2:2:(N_rand-1),1)+tmp_delta_d_mat(3:2:(N_rand)-1,1))/2;
        delta_cylinder_d_mat(2:(N_PAL_cols-1),N_PAL_cols)=(tmp_delta_d_mat(2:2:(N_rand-1),N_rand)+tmp_delta_d_mat(3:2:(N_rand)-1,N_rand))/2;
        
        %1 cylinder in the middle, i.e. 4 units
        tmp_int1=N_PAL_cols-2;
        tmp_int2=2*tmp_int1;
        tmp_mat1=diag(ones(1,tmp_int1));
        tmp_mat2=zeros(tmp_int1,tmp_int2);%initial in case size error
        tmp_mat2(:,1:2:(tmp_int2-1))=tmp_mat1;
        tmp_mat2(:,2:2:tmp_int2)=tmp_mat1;
        delta_cylinder_d_mat(2:(N_PAL_cols-1),2:(N_PAL_cols-1))=tmp_mat2*(tmp_delta_d_mat(2:(N_rand-1),2:(N_rand-1)))*(tmp_mat2')*1/4;
    end
    
    % random pos
    tmp_rand_pos_mat = rand(N_PAL_cols,N_PAL_cols).*delta_cylinder_d_mat;
    
    count_cyl=1; %% need to update: record components of PAL somewhere and load here.
    count_elp=1;
    for loop_i=1:N_PAL_rows
        for loop_j=1:N_PAL_cols
            tmpstr_displx=['(PAL_MS_radius*2+MODEL_ddis)*',num2str(loop_j-1)];
            tmpstr_disply=['(PAL_MS_radius*2+MODEL_ddis)*',num2str(loop_i-1)];
        
            tmptag_pal_cp1=['cyl', num2str(count_cyl)];
            count_cyl=count_cyl+1;
            tmpstr_displz=['SPO_box_thick+MODEL_ddis+(PAL_MS_radius)*(1-PAL_MS_lcap_flat)+',num2str(tmp_rand_pos_mat(loop_i,loop_j))];
            eval(['model.component(''comp1'').geom(''geom1'').feature(tmptag_pal_cp1).set(''pos'', {''',tmpstr_displx,''' ''',tmpstr_disply,''' ''',tmpstr_displz,'''});'])
            eval(['model.component(''comp1'').geom(''geom1'').feature(tmptag_pal_cp1).set(''h'', ''PAL_MS_height-',num2str(delta_cylinder_d_mat(loop_i,loop_j)),''');'])
            
            tmptag_pal_cp2=['elp', num2str(count_elp)];
            count_elp=count_elp+1;
            eval(['model.component(''comp1'').geom(''geom1'').feature(tmptag_pal_cp2).set(''pos'', {''',tmpstr_displx,''' ''',tmpstr_disply,''' ''',tmpstr_displz,'''});'])
            
            tmptag_pal_cp3=['elp', num2str(count_elp)];
            count_elp=count_elp+1;
            tmpstr_displz=['SPO_box_thick+MODEL_ddis+(PAL_MS_radius)*(1-PAL_MS_lcap_flat)+PAL_MS_height+',num2str(tmp_rand_pos_mat(loop_i,loop_j)),'-',num2str(delta_cylinder_d_mat(loop_i,loop_j))];
            eval(['model.component(''comp1'').geom(''geom1'').feature(tmptag_pal_cp3).set(''pos'', {''',tmpstr_displx,''' ''',tmpstr_disply,''' ''',tmpstr_displz,'''});'])
        end
    end
end
model.component('comp1').geom('geom1').run('dif1');%% IMPORTANT: run before measure
model.component('comp1').geom('geom1').measure.selection.init(3);
model.component('comp1').geom('geom1').measure.selection.set('dif1', 1);
model.component('comp1').geom('geom1').feature('dif1').set('keep', true);
tmp_vol_IAS_PAL_new = model.geom('geom1').measure().getVolume()/(PAL_box_thick*PAL_box_length*PAL_box_length)

save tmp_PAL_MS_distribute.mat -regexp '^(?!(model|ans)$).'
disp('Palisade IAS porosity matched')
mphsave(model, 'tmp_PAL_MS_geo.mph')
ModelUtil.remove('Model_PAL_MS_distribute');
