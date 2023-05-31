% eLeaf: 3D model of dicot leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 0.1b

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update 2020-Jan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output mph for e_RT_geo_export
% add more output Msg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update 2020-Oct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to develop:
% - {simple, debug} modes for output messaqes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
load parainput.mat

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model_3Dleaf_dicot');

model.param.set('EPL_u_thick', [num2str(EPL_u_thick),'[m]']);
model.param.set('EPL_u_radius', [num2str(EPL_u_radius), '[m]']);
model.param.set('EPL_u_flat', num2str(EPL_u_flat));
model.param.set('EPL_u_z', 'SPO_box_thick+PAL_box_thick+EPL_u_thick/2');
model.param.set('EPL_u_rand_displacex',[num2str(EPL_u_rand_displacex),'[m]']);
model.param.set('EPL_u_rand_displacey',[num2str(EPL_u_rand_displacey),'[m]']);
model.param.set('EPL_l_thick', [num2str(EPL_l_thick),'[m]']);
model.param.set('EPL_l_radius', [num2str(EPL_l_radius), '[m]']);
model.param.set('EPL_l_flat', num2str(EPL_l_flat));
model.param.set('EPL_l_z', '-EPL_l_thick/2');
model.param.set('EPL_l_rand_displacex',[num2str(EPL_l_rand_displacex),'[m]']);
model.param.set('EPL_l_rand_displacey',[num2str(EPL_l_rand_displacey),'[m]']);

model.param.set('PAL_box_thick', [num2str(PAL_box_thick),'[m]']);
model.param.set('PAL_MS_radius', [num2str(PAL_MS_radius),'[m]']);
%N_PAL_cols=4;
%N_PAL_rows=4;
model.param.set('PAL_box_length', [num2str(PAL_box_length),'[m]'], '4*4 and cut by box');%PAL_box_length=PAL_MS_radius*6;
model.param.set('PAL_MS_ucap_flat', num2str(PAL_MS_ucap_flat));
model.param.set('PAL_MS_lcap_flat', num2str(PAL_MS_lcap_flat));
model.param.set('MODEL_ddis', [num2str(MODEL_ddis),'[m]']);
model.param.set('PAL_MS_height', 'PAL_box_thick-MODEL_ddis*2-(PAL_MS_radius)*(1-PAL_MS_ucap_flat)-(PAL_MS_radius)*(1-PAL_MS_lcap_flat)');
model.param.set('SPO_box_thick', [num2str(SPO_box_thick),'[m]']);
model.param.set('SPO_box_length', [num2str(SPO_box_length),'[m]']);%SPO_box_length=PAL_box_length;

model.param.set('SPO_MS_radius', [num2str(SPO_MS_radius),'[m]']);

%%%for PAL CHL
model.param.set('PAL_dis_chl2wall', [num2str(PAL_dis_chl2wall),'[m]']);
model.param.set('PAL_thick_chl_layer', [num2str(PAL_thick_chl_layer),'[m]']);
%model.param.set('PAL_diameter_chl_patch', [num2str(PAL_diameter_chl_patch),'[m]']);
model.param.set('PAL_radius_chl_sphere', [num2str(PAL_radius_chl_sphere),'[m]']);
%%%for PAL MIT & VAC
model.param.set('PAL_dis_mit2chl_rt', num2str(PAL_dis_mit2chl_rt));%relative to PAL_MS_radius
model.param.set('PAL_radius_mit_rt', num2str(PAL_radius_mit_rt));
model.param.set('PAL_dis_vac2mit_rt', num2str(PAL_dis_vac2mit_rt));
%model.param.set('PAL_diameter_mit_patch', [num2str(PAL_diameter_mit_patch),'[m]']);
%%%for SPO CHL
model.param.set('SPO_dis_chl2wall', [num2str(SPO_dis_chl2wall),'[m]']);
model.param.set('SPO_thick_chl_layer', [num2str(SPO_thick_chl_layer),'[m]']);
%model.param.set('SPO_diameter_chl_patch', [num2str(SPO_diameter_chl_patch),'[m]']);
model.param.set('SPO_radius_chl_sphere', [num2str(SPO_radius_chl_sphere),'[m]']);
%%%for SPO MIT & VAC
model.param.set('SPO_dis_mit2chl_rt', num2str(SPO_dis_mit2chl_rt));%relative to SPO_MS_radius
model.param.set('SPO_radius_mit_rt', num2str(SPO_radius_mit_rt));
model.param.set('SPO_dis_vac2mit_rt', num2str(SPO_dis_vac2mit_rt));
%model.param.set('SPO_diameter_mit_patch', [num2str(SPO_diameter_mit_patch),'[m]']);

model.component.create('comp1', true);
model.component('comp1').geom.create('geom1', 3);

%% initialization
%N_EPL_cols=5;
%N_EPL_rows=4;
count_copy=1;
count_del=1;
count_cyl=1;
count_elp=1;
count_uni=1;
count_sph=1;
count_dif=1;
count_wp=1;%work plane
count_par=1;%operation partition
count_spl=1;%
count_int=1;
count_blk=1;

%% lower epidermis
tmptag_blk=['blk',num2str(count_blk)];count_blk=count_blk+1;
model.component('comp1').geom('geom1').create(tmptag_blk, 'Block');
model.component('comp1').geom('geom1').feature(tmptag_blk).set('pos', {'0' '0' '-EPL_l_thick*1.25'});
model.component('comp1').geom('geom1').feature(tmptag_blk).set('size', {'PAL_box_length' 'PAL_box_length' 'EPL_l_thick*1.5'});
GLB_tag_EPL_l_blk=tmptag_blk;
GLB_tag_EPL_l_set={};tmp_count=1;

tmptag_wp=['wp',num2str(count_wp)];count_wp=count_wp+1;
model.component('comp1').geom('geom1').create(tmptag_wp, 'WorkPlane');
model.component('comp1').geom('geom1').feature(tmptag_wp).set('quickz', '-EPL_l_thick');
model.component('comp1').geom('geom1').feature(tmptag_wp).geom.create('pol1', 'Polygon');
model.component('comp1').geom('geom1').feature(tmptag_wp).geom.feature('pol1').set('source', 'table');
model.component('comp1').geom('geom1').feature(tmptag_wp).geom.feature('pol1').set('table', {'EPL_l_radius' '0';  ...
'1/2*EPL_l_radius' '1/2*sqrt(3)*EPL_l_radius';  ...
'-1/2*EPL_l_radius' '1/2*sqrt(3)*EPL_l_radius';  ...
'-EPL_l_radius' '0';  ...
'-1/2*EPL_l_radius' '-1/2*sqrt(3)*EPL_l_radius';  ...
'1/2*EPL_l_radius' '-1/2*sqrt(3)*EPL_l_radius'});
model.component('comp1').geom('geom1').create('ext1', 'Extrude');
model.component('comp1').geom('geom1').feature('ext1').setIndex('distance', 'EPL_l_thick*2', 0);
model.component('comp1').geom('geom1').feature('ext1').selection('input').set({'wp1'});
tmptag_elp=['elp',num2str(count_elp)];count_elp=count_elp+1;
model.component('comp1').geom('geom1').create(tmptag_elp, 'Ellipsoid');
model.component('comp1').geom('geom1').feature(tmptag_elp).set('semiaxes', {'EPL_l_thick/2/(1-EPL_l_flat)' 'EPL_l_thick/2/(1-EPL_l_flat)' 'EPL_l_thick/2'});
tmptag_int=['int',num2str(count_int)];count_int=count_int+1;
model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set([cellstr(tmptag_elp), cellstr('ext1')]);
TMP_tag_EPL_l=tmptag_int;
TMP_tag_EPL_l_set_final_del={};tmp_count4TMP_tag=1;

for loop_i=1:N_EPL_l_rows
    for loop_j=1:N_EPL_l_cols
        tmptag_epl=['copy', num2str(count_copy)];
        model.component('comp1').geom('geom1').create(tmptag_epl, 'Copy');
        count_copy=count_copy+1;
        eval(['model.component(''comp1'').geom(''geom1'').feature(tmptag_epl).set(''displx'', ''1.5*',num2str(loop_j-1),'*EPL_l_radius-EPL_l_rand_displacex'');'])
        eval(['model.component(''comp1'').geom(''geom1'').feature(tmptag_epl).set(''disply'', ''sqrt(3)*(',num2str((loop_i-1)+0.5*mod(loop_j-1,2)),')*EPL_l_radius-EPL_l_rand_displacey'');'])
        model.component('comp1').geom('geom1').feature(tmptag_epl).set('displz', 'EPL_l_z');
        model.component('comp1').geom('geom1').feature(tmptag_epl).selection('input').set(TMP_tag_EPL_l);
        
        tmptag_int=['int',num2str(count_int)];count_int=count_int+1;
        model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
        model.component('comp1').geom('geom1').feature(tmptag_int).set('keep', true);
        model.component('comp1').geom('geom1').feature(tmptag_int).set('intbnd', false);
        model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set([cellstr(GLB_tag_EPL_l_blk), cellstr(tmptag_epl)]);
        %check void domain
        model.component('comp1').geom('geom1').run(tmptag_int);
        if model.geom('geom1').obj(tmptag_int).getNDomains==1
            GLB_tag_EPL_l_set(tmp_count)=cellstr(tmptag_int);
            tmp_count=tmp_count+1;
            TMP_tag_EPL_l_set_final_del(tmp_count4TMP_tag)=cellstr(tmptag_epl);
            tmp_count4TMP_tag=tmp_count4TMP_tag+1;
        else
            model.component('comp1').geom('geom1').feature.remove(tmptag_int);
            model.component('comp1').geom('geom1').feature.remove(tmptag_epl);
        end
    end
end
tmptag_del=['del',num2str(count_del)];count_del=count_del+1;
model.component('comp1').geom('geom1').create(tmptag_del, 'Delete');
model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').init;
model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set([cellstr(TMP_tag_EPL_l),TMP_tag_EPL_l_set_final_del,cellstr(GLB_tag_EPL_l_blk)]);

%% upper epidermis
tmptag_blk=['blk',num2str(count_blk)];count_blk=count_blk+1;
model.component('comp1').geom('geom1').create(tmptag_blk, 'Block');
model.component('comp1').geom('geom1').feature(tmptag_blk).set('pos', {'0' '0' 'SPO_box_thick+PAL_box_thick'});
model.component('comp1').geom('geom1').feature(tmptag_blk).set('size', {'PAL_box_length' 'PAL_box_length' 'EPL_u_thick*1.25'});
GLB_tag_EPL_u_blk=tmptag_blk;
GLB_tag_EPL_u_set={};tmp_count=1;
TMP_tag_EPL_u_set_final_del={};tmp_count4TMP_tag=1;

tmptag_wp=['wp',num2str(count_wp)];count_wp=count_wp+1;
model.component('comp1').geom('geom1').create(tmptag_wp, 'WorkPlane');
model.component('comp1').geom('geom1').feature(tmptag_wp).set('quickz', '-EPL_u_thick');
model.component('comp1').geom('geom1').feature(tmptag_wp).geom.create('pol1', 'Polygon');
model.component('comp1').geom('geom1').feature(tmptag_wp).geom.feature('pol1').set('source', 'table');
model.component('comp1').geom('geom1').feature(tmptag_wp).geom.feature('pol1').set('table', {'EPL_u_radius' '0';  ...
'1/2*EPL_u_radius' '1/2*sqrt(3)*EPL_u_radius';  ...
'-1/2*EPL_u_radius' '1/2*sqrt(3)*EPL_u_radius';  ...
'-EPL_u_radius' '0';  ...
'-1/2*EPL_u_radius' '-1/2*sqrt(3)*EPL_u_radius';  ...
'1/2*EPL_u_radius' '-1/2*sqrt(3)*EPL_u_radius'});
model.component('comp1').geom('geom1').create('ext2', 'Extrude');
model.component('comp1').geom('geom1').feature('ext2').setIndex('distance', 'EPL_u_thick*2', 0);
model.component('comp1').geom('geom1').feature('ext2').selection('input').set({'wp2'});
tmptag_elp=['elp',num2str(count_elp)];count_elp=count_elp+1;
model.component('comp1').geom('geom1').create(tmptag_elp, 'Ellipsoid');
model.component('comp1').geom('geom1').feature(tmptag_elp).set('semiaxes', {'EPL_u_thick/2/(1-EPL_u_flat)' 'EPL_u_thick/2/(1-EPL_u_flat)' 'EPL_u_thick/2'});
tmptag_int=['int',num2str(count_int)];count_int=count_int+1;
model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set([cellstr(tmptag_elp), cellstr('ext2')]);
TMP_tag_EPL_u=tmptag_int;

for loop_i=1:N_EPL_u_rows
    for loop_j=1:N_EPL_u_cols
        tmptag_epl=['copy', num2str(count_copy)];count_copy=count_copy+1;
        model.component('comp1').geom('geom1').create(tmptag_epl, 'Copy');
        eval(['model.component(''comp1'').geom(''geom1'').feature(tmptag_epl).set(''displx'', ''1.5*',num2str(loop_j-1),'*EPL_u_radius-EPL_u_rand_displacex'');'])
        eval(['model.component(''comp1'').geom(''geom1'').feature(tmptag_epl).set(''disply'', ''sqrt(3)*(',num2str((loop_i-1)+0.5*mod(loop_j-1,2)),')*EPL_u_radius-EPL_u_rand_displacey'');'])
        model.component('comp1').geom('geom1').feature(tmptag_epl).set('displz', 'EPL_u_z');
        model.component('comp1').geom('geom1').feature(tmptag_epl).selection('input').set(TMP_tag_EPL_u);
        
        tmptag_int=['int',num2str(count_int)];count_int=count_int+1;
        model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
        model.component('comp1').geom('geom1').feature(tmptag_int).set('keep', true);
        model.component('comp1').geom('geom1').feature(tmptag_int).set('intbnd', false);
        model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set([cellstr(GLB_tag_EPL_u_blk), cellstr(tmptag_epl)]);
        %check void domain
        model.component('comp1').geom('geom1').run(tmptag_int);
        if model.geom('geom1').obj(tmptag_int).getNDomains==1
            GLB_tag_EPL_u_set(tmp_count)=cellstr(tmptag_int);
            tmp_count=tmp_count+1;
            TMP_tag_EPL_u_set_final_del(tmp_count4TMP_tag)=cellstr(tmptag_epl);
            tmp_count4TMP_tag=tmp_count4TMP_tag+1;
        else
            model.component('comp1').geom('geom1').feature.remove(tmptag_int);
            model.component('comp1').geom('geom1').feature.remove(tmptag_epl);
        end
    end
end
tmptag_del=['del',num2str(count_del)];count_del=count_del+1;
model.component('comp1').geom('geom1').create(tmptag_del, 'Delete');
model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').init;
model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set([cellstr(TMP_tag_EPL_u),TMP_tag_EPL_u_set_final_del,cellstr(GLB_tag_EPL_u_blk)]);

%% palisade
load('tmp_PAL_MS_distribute.mat', 'tmp_rand_pos_mat', 'delta_cylinder_d_mat')

tmptag_blk=['blk',num2str(count_blk)];count_blk=count_blk+1;
model.component('comp1').geom('geom1').create(tmptag_blk, 'Block');
model.component('comp1').geom('geom1').feature(tmptag_blk).set('pos', {'0' '0' 'SPO_box_thick'});
model.component('comp1').geom('geom1').feature(tmptag_blk).set('size', {'PAL_box_length' 'PAL_box_length' 'PAL_box_thick'});
GLB_tag_PAL_blk=tmptag_blk;
GLB_tag_PAL_MS={};
GLB_tag_PAL_CHL_mem={};
GLB_tag_PAL_VAC={};
count_PAL=1;
GLB_mat_contact_PAL_MSVAC_idx=[];

for loop_i=1:N_PAL_rows
    for loop_j=1:N_PAL_cols
        % PAL_MS_surface
        tmptag_pal_1=['cyl', num2str(count_cyl)];
        model.component('comp1').geom('geom1').create(tmptag_pal_1, 'Cylinder');
        count_cyl=count_cyl+1;
        tmpstr_displx=['(PAL_MS_radius*2+MODEL_ddis)*',num2str(loop_j-1)];
        tmpstr_disply=['(PAL_MS_radius*2+MODEL_ddis)*',num2str(loop_i-1)];
        tmpstr_displz=['SPO_box_thick+MODEL_ddis+(PAL_MS_radius)*(1-PAL_MS_lcap_flat)+',num2str(tmp_rand_pos_mat(loop_i,loop_j))];
        eval(['model.component(''comp1'').geom(''geom1'').feature(tmptag_pal_1).set(''pos'', {''',tmpstr_displx,''' ''',tmpstr_disply,''' ''',tmpstr_displz,'''});'])
        model.component('comp1').geom('geom1').feature(tmptag_pal_1).set('r', 'PAL_MS_radius');
        eval(['model.component(''comp1'').geom(''geom1'').feature(tmptag_pal_1).set(''h'', ''PAL_MS_height-',num2str(delta_cylinder_d_mat(loop_i,loop_j)),''');'])
        
        tmptag_pal_2=['elp', num2str(count_elp)];
        model.component('comp1').geom('geom1').create(tmptag_pal_2, 'Ellipsoid');
        count_elp=count_elp+1;
        eval(['model.component(''comp1'').geom(''geom1'').feature(tmptag_pal_2).set(''pos'', {''',tmpstr_displx,''' ''',tmpstr_disply,''' ''',tmpstr_displz,'''});'])
        model.component('comp1').geom('geom1').feature(tmptag_pal_2).set('semiaxes', {'PAL_MS_radius' 'PAL_MS_radius' '(PAL_MS_radius)*(1-PAL_MS_lcap_flat)'});
        
        tmptag_pal_3=['elp', num2str(count_elp)];
        count_elp=count_elp+1;
        tmpstr_displz=['SPO_box_thick+MODEL_ddis+(PAL_MS_radius)*(1-PAL_MS_lcap_flat)+PAL_MS_height+',num2str(tmp_rand_pos_mat(loop_i,loop_j)),'-',num2str(delta_cylinder_d_mat(loop_i,loop_j))];
        model.component('comp1').geom('geom1').create(tmptag_pal_3, 'Ellipsoid');
        eval(['model.component(''comp1'').geom(''geom1'').feature(tmptag_pal_3).set(''pos'', {''',tmpstr_displx,''' ''',tmpstr_disply,''' ''',tmpstr_displz,'''});'])
        model.component('comp1').geom('geom1').feature(tmptag_pal_3).set('semiaxes', {'PAL_MS_radius' 'PAL_MS_radius' '(PAL_MS_radius)*(1-PAL_MS_ucap_flat)'});
        
        tmptag_uni=['uni', num2str(count_uni)];
        model.component('comp1').geom('geom1').create(tmptag_uni, 'Union');
        count_uni=count_uni+1;
        model.component('comp1').geom('geom1').feature(tmptag_uni).set('intbnd', false);
        tmptag_set={tmptag_pal_1,tmptag_pal_2,tmptag_pal_3};
        model.component('comp1').geom('geom1').feature(tmptag_uni).selection('input').set(tmptag_set);
        
        GLB_tag_PAL_MS{count_PAL}=tmptag_uni;
        %count_PAL=count_PAL+1;
        
        % PAL CHL membrane outer surface
        tmptag_pal_1=['cyl', num2str(count_cyl)];
        model.component('comp1').geom('geom1').create(tmptag_pal_1, 'Cylinder');
        count_cyl=count_cyl+1;
        %tmpstr_displx=['(PAL_MS_radius*2+MODEL_ddis)*',num2str(loop_j-1)];
        %tmpstr_disply=['(PAL_MS_radius*2+MODEL_ddis)*',num2str(loop_i-1)];
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
        tmptag_PAL_mem_o=tmptag_uni;
        
        % PAL CHL membrane inner surface
        tmptag_pal_1=['cyl', num2str(count_cyl)];
        model.component('comp1').geom('geom1').create(tmptag_pal_1, 'Cylinder');
        count_cyl=count_cyl+1;
        %tmpstr_displx=['(PAL_MS_radius*2+MODEL_ddis)*',num2str(loop_j-1)];
        %tmpstr_disply=['(PAL_MS_radius*2+MODEL_ddis)*',num2str(loop_i-1)];
        tmpstr_displz=['SPO_box_thick+MODEL_ddis+(PAL_MS_radius)*(1-PAL_MS_lcap_flat)+PAL_dis_chl2wall+PAL_thick_chl_layer+',num2str(tmp_rand_pos_mat(loop_i,loop_j))];
        eval(['model.component(''comp1'').geom(''geom1'').feature(tmptag_pal_1).set(''pos'', {''',tmpstr_displx,''' ''',tmpstr_disply,''' ''',tmpstr_displz,'''});'])
        model.component('comp1').geom('geom1').feature(tmptag_pal_1).set('r', 'PAL_MS_radius-PAL_dis_chl2wall-PAL_thick_chl_layer');
        eval(['model.component(''comp1'').geom(''geom1'').feature(tmptag_pal_1).set(''h'', ''PAL_MS_height-2*PAL_dis_chl2wall-2*PAL_thick_chl_layer-',num2str(delta_cylinder_d_mat(loop_i,loop_j)),''');'])
        
        tmptag_pal_2=['elp', num2str(count_elp)];
        model.component('comp1').geom('geom1').create(tmptag_pal_2, 'Ellipsoid');
        count_elp=count_elp+1;
        eval(['model.component(''comp1'').geom(''geom1'').feature(tmptag_pal_2).set(''pos'', {''',tmpstr_displx,''' ''',tmpstr_disply,''' ''',tmpstr_displz,'''});'])
        model.component('comp1').geom('geom1').feature(tmptag_pal_2).set('semiaxes', {'PAL_MS_radius-PAL_dis_chl2wall-PAL_thick_chl_layer' 'PAL_MS_radius-PAL_dis_chl2wall-PAL_thick_chl_layer' '(PAL_MS_radius-PAL_dis_chl2wall-PAL_thick_chl_layer)*(1-PAL_MS_lcap_flat)'});
        
        tmptag_pal_3=['elp', num2str(count_elp)];
        count_elp=count_elp+1;
        tmpstr_displz=['SPO_box_thick+MODEL_ddis+(PAL_MS_radius)*(1-PAL_MS_lcap_flat)+PAL_MS_height+',num2str(tmp_rand_pos_mat(loop_i,loop_j)),'-PAL_dis_chl2wall-PAL_thick_chl_layer-',num2str(delta_cylinder_d_mat(loop_i,loop_j))];
        model.component('comp1').geom('geom1').create(tmptag_pal_3, 'Ellipsoid');
        eval(['model.component(''comp1'').geom(''geom1'').feature(tmptag_pal_3).set(''pos'', {''',tmpstr_displx,''' ''',tmpstr_disply,''' ''',tmpstr_displz,'''});'])
        model.component('comp1').geom('geom1').feature(tmptag_pal_3).set('semiaxes', {'PAL_MS_radius-PAL_dis_chl2wall-PAL_thick_chl_layer' 'PAL_MS_radius-PAL_dis_chl2wall-PAL_thick_chl_layer' '(PAL_MS_radius-PAL_dis_chl2wall-PAL_thick_chl_layer)*(1-PAL_MS_ucap_flat)'});
        
        tmptag_uni=['uni', num2str(count_uni)];
        model.component('comp1').geom('geom1').create(tmptag_uni, 'Union');
        count_uni=count_uni+1;
        model.component('comp1').geom('geom1').feature(tmptag_uni).set('intbnd', false);
        tmptag_set={tmptag_pal_1,tmptag_pal_2,tmptag_pal_3};
        model.component('comp1').geom('geom1').feature(tmptag_uni).selection('input').set(tmptag_set);
        tmptag_PAL_mem_i=tmptag_uni;
        
        % create PAL CHL membrane
        tmptag_dif=['dif',num2str(count_dif)];count_dif=count_dif+1;
        model.component('comp1').geom('geom1').create(tmptag_dif, 'Difference');
        model.component('comp1').geom('geom1').feature(tmptag_dif).set('intbnd', false);
        model.component('comp1').geom('geom1').feature(tmptag_dif).selection('input').set(tmptag_PAL_mem_o);%set({'blk1'})
        model.component('comp1').geom('geom1').feature(tmptag_dif).selection('input2').set(tmptag_PAL_mem_i);
        
		GLB_tag_PAL_CHL_mem{count_PAL}=tmptag_dif;
		
		% PAL VAC
		tmptag_pal_1=['cyl', num2str(count_cyl)];
        model.component('comp1').geom('geom1').create(tmptag_pal_1, 'Cylinder');
        count_cyl=count_cyl+1;
        %tmpstr_displx=['(PAL_MS_radius*2+MODEL_ddis)*',num2str(loop_j-1)];
        %tmpstr_disply=['(PAL_MS_radius*2+MODEL_ddis)*',num2str(loop_i-1)];
        tmpstr_displz=['SPO_box_thick+MODEL_ddis+(PAL_MS_radius)*(1-PAL_MS_lcap_flat)+PAL_dis_chl2wall+PAL_thick_chl_layer+PAL_MS_radius*(PAL_dis_mit2chl_rt+PAL_dis_vac2mit_rt)+',num2str(tmp_rand_pos_mat(loop_i,loop_j))];
        eval(['model.component(''comp1'').geom(''geom1'').feature(tmptag_pal_1).set(''pos'', {''',tmpstr_displx,''' ''',tmpstr_disply,''' ''',tmpstr_displz,'''});'])
        model.component('comp1').geom('geom1').feature(tmptag_pal_1).set('r', 'PAL_MS_radius-PAL_dis_chl2wall-PAL_thick_chl_layer-PAL_MS_radius*(PAL_dis_mit2chl_rt+PAL_dis_vac2mit_rt)');
        eval(['model.component(''comp1'').geom(''geom1'').feature(tmptag_pal_1).set(''h'', ''PAL_MS_height-2*PAL_dis_chl2wall-2*PAL_thick_chl_layer-2*PAL_MS_radius*(PAL_dis_mit2chl_rt+PAL_dis_vac2mit_rt)-',num2str(delta_cylinder_d_mat(loop_i,loop_j)),''');'])
        
        tmptag_pal_2=['elp', num2str(count_elp)];
        model.component('comp1').geom('geom1').create(tmptag_pal_2, 'Ellipsoid');
        count_elp=count_elp+1;
        eval(['model.component(''comp1'').geom(''geom1'').feature(tmptag_pal_2).set(''pos'', {''',tmpstr_displx,''' ''',tmpstr_disply,''' ''',tmpstr_displz,'''});'])
        model.component('comp1').geom('geom1').feature(tmptag_pal_2).set('semiaxes', {'PAL_MS_radius-PAL_dis_chl2wall-PAL_thick_chl_layer-PAL_MS_radius*(PAL_dis_mit2chl_rt+PAL_dis_vac2mit_rt)' 'PAL_MS_radius-PAL_dis_chl2wall-PAL_thick_chl_layer-PAL_MS_radius*(PAL_dis_mit2chl_rt+PAL_dis_vac2mit_rt)' '(PAL_MS_radius-PAL_dis_chl2wall-PAL_thick_chl_layer-PAL_MS_radius*(PAL_dis_mit2chl_rt+PAL_dis_vac2mit_rt))*(1-PAL_MS_lcap_flat)'});
        
        tmptag_pal_3=['elp', num2str(count_elp)];
        count_elp=count_elp+1;
        tmpstr_displz=['SPO_box_thick+MODEL_ddis+(PAL_MS_radius)*(1-PAL_MS_lcap_flat)+PAL_MS_height+',num2str(tmp_rand_pos_mat(loop_i,loop_j)),'-PAL_dis_chl2wall-PAL_thick_chl_layer-PAL_MS_radius*(PAL_dis_mit2chl_rt+PAL_dis_vac2mit_rt)-',num2str(delta_cylinder_d_mat(loop_i,loop_j))];
        model.component('comp1').geom('geom1').create(tmptag_pal_3, 'Ellipsoid');
        eval(['model.component(''comp1'').geom(''geom1'').feature(tmptag_pal_3).set(''pos'', {''',tmpstr_displx,''' ''',tmpstr_disply,''' ''',tmpstr_displz,'''});'])
        model.component('comp1').geom('geom1').feature(tmptag_pal_3).set('semiaxes', {'PAL_MS_radius-PAL_dis_chl2wall-PAL_thick_chl_layer-PAL_MS_radius*(PAL_dis_mit2chl_rt+PAL_dis_vac2mit_rt)' 'PAL_MS_radius-PAL_dis_chl2wall-PAL_thick_chl_layer-PAL_MS_radius*(PAL_dis_mit2chl_rt+PAL_dis_vac2mit_rt)' '(PAL_MS_radius-PAL_dis_chl2wall-PAL_thick_chl_layer-PAL_MS_radius*(PAL_dis_mit2chl_rt+PAL_dis_vac2mit_rt))*(1-PAL_MS_ucap_flat)'});
        
        tmptag_uni=['uni', num2str(count_uni)];
        model.component('comp1').geom('geom1').create(tmptag_uni, 'Union');
        count_uni=count_uni+1;
        model.component('comp1').geom('geom1').feature(tmptag_uni).set('intbnd', false);
        tmptag_set={tmptag_pal_1,tmptag_pal_2,tmptag_pal_3};
        model.component('comp1').geom('geom1').feature(tmptag_uni).selection('input').set(tmptag_set);
        
        GLB_tag_PAL_VAC{count_PAL}=tmptag_uni;
        
        %%% intersect with box?
        %tmpdbl_displx=(PAL_MS_radius*2+MODEL_ddis)*(loop_j-1);
        %tmpdbl_disply=(PAL_MS_radius*2+MODEL_ddis)*(loop_i-1);
        if (loop_i==1) ...,
                || (loop_i==N_PAL_rows) ...,
                || (loop_j==1) ...,
                || (loop_j==N_PAL_cols)
            GLB_mat_contact_PAL_MSVAC_idx(count_PAL)=1;
        else
            GLB_mat_contact_PAL_MSVAC_idx(count_PAL)=0;
        end
        
        count_PAL=count_PAL+1;
    end
end

%% PAL CHL
load('tmp_PAL_CHL_distribute.mat', 'PAL_CHL_center_final_pts','PAL_CHL_sph_selected_idx','PAL_CHL_sph_notselected_idx')
%PAL_CHL_center_final_pts_selected=PAL_CHL_center_final_pts(PAL_CHL_sph_selected_idx,:);
GLB_tag_PAL_CHL={};%TYPE=cell array of cell array
GLB_mat_contact_PAL_CHL_idx={};%TYPE=cell array of tmp_mat_contact_PAL_CHL_idx
Num_pal=numel(GLB_tag_PAL_MS);
disp('[Msg] start process PAL CHL')
for loop_i=1:Num_pal
    tmp_pts=PAL_CHL_center_final_pts{loop_i}(PAL_CHL_sph_selected_idx{loop_i},:);
    %%%%create CHL sphere
    tmptag_set_PAL_CHL={};
    for loop_j=1:size(tmp_pts,1)
        tmptag_sph=['sph',num2str(count_sph)];count_sph=count_sph+1;
        model.component('comp1').geom('geom1').create(tmptag_sph, 'Sphere');
        model.component('comp1').geom('geom1').feature(tmptag_sph).set('r', 'PAL_radius_chl_sphere');
        model.component('comp1').geom('geom1').feature(tmptag_sph).set('pos', string(tmp_pts(loop_j,:)));%{'5e-6' '0' '150e-6'}
        tmptag_set_PAL_CHL{loop_j}=tmptag_sph;
    end
    model.component('comp1').geom('geom1').run(tmptag_sph);
%     %merge CHL sphere %% used in distribute_MS_CHL script for calculating coverage
%     tmptag_uni=['uni', num2str(count_uni)];count_uni=count_uni+1;
%     model.component('comp1').geom('geom1').create(tmptag_uni, 'Union');
%     model.component('comp1').geom('geom1').feature(tmptag_uni).set('intbnd', false);
%     model.component('comp1').geom('geom1').feature(tmptag_uni).selection('input').set(tmptag_set);
%     %calculate surface area
%     model.component('comp1').geom('geom1').measure.selection.init(2);
%     tmp_NFaces=model.geom('geom1').obj(PAL_CH_tag_set{loop_i}).getNFaces;
%     model.component('comp1').geom('geom1').measure.selection.set(PAL_CH_tag_set{loop_i},1:tmp_NFaces);%set('int1',[1,2,3]);%[1,2,3] is corresponding boundary
%     full_area(loop_i)=model.geom('geom1').measure().getVolume();%Remark: getVolume here actually getArea
%     %int operation
%     tmptag_int=['int',num2str(count_int)];count_int=count_int+1;
%     model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
%     model.component('comp1').geom('geom1').feature(tmptag_int).set('keep', true);
%     model.component('comp1').geom('geom1').feature(tmptag_int).set('intbnd', false);
%     tmpcellstr=[PAL_CH_tag_set(loop_i),cellstr(tmptag_uni)];
%     model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set(tmpcellstr);
%     %calculate projected area
%     model.component('comp1').geom('geom1').run(tmptag_int);
%     model.component('comp1').geom('geom1').measure.selection.init(2);
%     tmp_NFaces=model.geom('geom1').obj(tmptag_int).getNFaces;%%alternative method: meshstats - getNEntities()
%     model.component('comp1').geom('geom1').measure.selection.set(tmptag_int,1:tmp_NFaces);%set('int1',[1,2,3]);%[1,2,3] is corresponding boundary
%     project_area(loop_i)=model.geom('geom1').measure().getVolume();
%     coverage(loop_i)=project_area(loop_i)/full_area(loop_i);

    %%%% chop chloroplast sphere
    %%%%%% store contact information into a matrix mat_contact_cell_idx;
    %%%%%% if ith cell overlaps jth cell, then mat_contact_cell_idx(i,j)=1
    %%%%%% if ith cell overlaps with box, then mat_contact_cell_idx(i,i)=1
    %%%%%% mat_contact_cell_idx is symetric
    tmp_mat_contact_PAL_CHL_idx=zeros(size(tmp_pts,1));% initial
    tmp_radius_2=(PAL_radius_chl_sphere)^2;
    for loop_i_CHL=1:size(tmp_pts,1)
        for loop_j_CHL=loop_i_CHL:size(tmp_pts,1)
            if loop_j_CHL==loop_i_CHL
                % contact with box?
                if (tmp_pts(loop_i_CHL,1)<PAL_radius_chl_sphere) ...,
                        || (tmp_pts(loop_i_CHL,1)>SPO_box_length-(PAL_radius_chl_sphere)) ...,
                        || (tmp_pts(loop_i_CHL,2)<PAL_radius_chl_sphere) ...,
                        || (tmp_pts(loop_i_CHL,2)>SPO_box_length-(PAL_radius_chl_sphere)) ...,
                        || (tmp_pts(loop_i_CHL,3)<PAL_radius_chl_sphere) ...,
                        || (tmp_pts(loop_i_CHL,3)>SPO_box_thick-(PAL_radius_chl_sphere))
                    tmp_mat_contact_PAL_CHL_idx(loop_i_CHL,loop_i_CHL)=1;
                end
            else
                if sum((tmp_pts(loop_i_CHL,:)-tmp_pts(loop_j_CHL,:)).^2) < tmp_radius_2*4 % (2*radius)^2
                    tmp_mat_contact_PAL_CHL_idx(loop_i_CHL,loop_j_CHL)=1;
                    tmp_mat_contact_PAL_CHL_idx(loop_j_CHL,loop_i_CHL)=1;
                end
            end
        end
    end
    %save box contact information for later process (Line???)
    GLB_mat_contact_PAL_CHL_idx{loop_i}=tmp_mat_contact_PAL_CHL_idx;
    
    %%%%%% process each CHL sphere with information from tmp_mat_contact_PAL_CHL_idx
    for loop_i_CHL=1:size(tmp_pts,1)
        tmp_idx_contact=find(tmp_mat_contact_PAL_CHL_idx(loop_i_CHL,:)==1);
        for tmp_loop_j_CHL=1:numel(tmp_idx_contact)
            loop_j_CHL=tmp_idx_contact(tmp_loop_j_CHL);
            if loop_j_CHL > loop_i_CHL
                tmp_center_a=tmp_pts(loop_i_CHL,:);
                tmp_center_b=tmp_pts(loop_j_CHL,:);
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
                
                tmptag_wp=['wp', num2str(count_wp)]; count_wp=count_wp+1;
                model.component('comp1').geom('geom1').create(tmptag_wp, 'WorkPlane');
                model.component('comp1').geom('geom1').feature(tmptag_wp).set('planetype', 'coordinates');
                model.component('comp1').geom('geom1').feature(tmptag_wp).set('genpoints', [tmp_pt1; tmp_pt2; tmp_pt3]);
                model.component('comp1').geom('geom1').feature(tmptag_wp).set('unite', true);
                %partition both obj with workplane
                %target = obj loop_i
                tmptag_obj=tmptag_set_PAL_CHL{loop_i_CHL};%{} since it's a cell array
                tmptag_par=['par', num2str(count_par)];count_par=count_par+1;
                model.component('comp1').geom('geom1').create(tmptag_par, 'Partition');
                model.component('comp1').geom('geom1').feature(tmptag_par).set('partitionwith', 'workplane');
                model.component('comp1').geom('geom1').feature(tmptag_par).set('workplane', tmptag_wp);
                model.component('comp1').geom('geom1').feature(tmptag_par).selection('input').set(tmptag_obj);
                %split par
                tmptag_spl=['spl', num2str(count_spl)];count_spl=count_spl+1;
                model.component('comp1').geom('geom1').create(tmptag_spl, 'Split');
                model.component('comp1').geom('geom1').feature(tmptag_spl).selection('input').set({tmptag_par});
                %[check]% rare case split still result 1 domain
                model.component('comp1').geom('geom1').run(tmptag_spl);
                if numel(model.geom('geom1').feature(tmptag_spl).objectNames())==1
                    tmptag_set_PAL_CHL(loop_i_CHL)=cellstr(tmptag_spl);
                else
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
                        tmptag_set_PAL_CHL(loop_i_CHL)=cellstr(tmptag_spl_dom2);
                    else
                        model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(tmptag_spl_dom2);
                        tmptag_set_PAL_CHL(loop_i_CHL)=cellstr(tmptag_spl_dom1);
                    end
                end
                
                %target = obj loop_j
                tmptag_obj=tmptag_set_PAL_CHL{loop_j_CHL};%{} since it's a cell array
                tmptag_par=['par', num2str(count_par)];count_par=count_par+1;
                model.component('comp1').geom('geom1').create(tmptag_par, 'Partition');
                model.component('comp1').geom('geom1').feature(tmptag_par).set('partitionwith', 'workplane');
                model.component('comp1').geom('geom1').feature(tmptag_par).set('workplane', tmptag_wp);
                model.component('comp1').geom('geom1').feature(tmptag_par).selection('input').set(tmptag_obj);
                %split par
                tmptag_spl=['spl', num2str(count_spl)];count_spl=count_spl+1;
                model.component('comp1').geom('geom1').create(tmptag_spl, 'Split');
                model.component('comp1').geom('geom1').feature(tmptag_spl).selection('input').set({tmptag_par});
                %[check]% rare case split still result 1 domain
                model.component('comp1').geom('geom1').run(tmptag_spl);
                if numel(model.geom('geom1').feature(tmptag_spl).objectNames())==1
                    tmptag_set_PAL_CHL(loop_j_CHL)=cellstr(tmptag_spl);
                else
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
                        tmptag_set_PAL_CHL(loop_j_CHL)=cellstr(tmptag_spl_dom2);
                    else
                        model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(tmptag_spl_dom2);
                        tmptag_set_PAL_CHL(loop_j_CHL)=cellstr(tmptag_spl_dom1);
                    end
                end
            end
        end

        %contact with box?
        %if tmp_mat_contact_PAL_CHL_idx(loop_i_CHL,loop_i_CHL)~=0   
        %end
    end
    
    
    %%%% intersect each chl sphere with chl membrane
    TMP_tag_set_final_del={};%tmp_count=1;
    tmp_CHL_sph_selected_idx=PAL_CHL_sph_selected_idx{loop_i};
    for loop_i_CHL=1:numel(tmp_CHL_sph_selected_idx)
        tmptag_obj=tmptag_set_PAL_CHL{loop_i_CHL};%{} since it's a cell array
        tmptag_int=['int',num2str(count_int)];count_int=count_int+1;
        model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
        model.component('comp1').geom('geom1').feature(tmptag_int).set('keep', true);
        model.component('comp1').geom('geom1').feature(tmptag_int).set('intbnd', false);
        model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set([GLB_tag_PAL_CHL_mem(loop_i), cellstr(tmptag_obj)]);
        tmptag_set_PAL_CHL{loop_i_CHL}=tmptag_int;
        TMP_tag_set_final_del{loop_i_CHL}=tmptag_obj;
    end
    GLB_tag_PAL_CHL{loop_i}=tmptag_set_PAL_CHL;
    
    tmptag_del=['del',num2str(count_del)];count_del=count_del+1;
    model.component('comp1').geom('geom1').create(tmptag_del, 'Delete');
    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').init;
    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set([TMP_tag_set_final_del,GLB_tag_PAL_CHL_mem(loop_i)]);
    
    disp(['[Msg] finish process PAL CHL No.',num2str(loop_i),'. In total ',num2str(Num_pal)])
end

%% spongy
load('tmp_SPO_MS_distribute.mat', 'SPO_MS_center_set_X', 'SPO_MS_center_set_Y', 'SPO_MS_center_set_Z', 'dis_SPO_MS_center')

tmptag_blk=['blk',num2str(count_blk)];count_blk=count_blk+1;
model.component('comp1').geom('geom1').create(tmptag_blk, 'Block');
model.component('comp1').geom('geom1').feature(tmptag_blk).set('pos', [0 0 0]);
model.component('comp1').geom('geom1').feature(tmptag_blk).set('size', {'SPO_box_length' 'SPO_box_length' 'SPO_box_thick'});
GLB_tag_SPO_blk=tmptag_blk;

Num_spo=numel(SPO_MS_center_set_X);
count_SPO=1;
GLB_tag_SPO_MS={};
GLB_tag_SPO_CHL_os={};
GLB_tag_SPO_CHL_is={};
GLB_tag_SPO_CHL_mem={};
GLB_tag_SPO_VAC={};
GLB_mat_contact_SPO_MSVAC_idx=[];

for loop=1:Num_spo
    tmp_center_X=SPO_MS_center_set_X(loop);
    tmp_center_Y=SPO_MS_center_set_Y(loop);
    tmp_center_Z=SPO_MS_center_set_Z(loop);
    tmptag_spo=['sph',num2str(count_sph)];
    model.component('comp1').geom('geom1').create(tmptag_spo, 'Sphere');
    count_sph=count_sph+1;
    model.component('comp1').geom('geom1').feature(tmptag_spo).set('r', 'SPO_MS_radius');
    tmpstr_pos={num2str(tmp_center_X), num2str(tmp_center_Y), num2str(tmp_center_Z)};
    model.component('comp1').geom('geom1').feature(tmptag_spo).set('pos', tmpstr_pos);    
    GLB_tag_SPO_MS{count_SPO}=tmptag_spo;
    
    tmptag_spo_os=['sph',num2str(count_sph)];
    model.component('comp1').geom('geom1').create(tmptag_spo_os, 'Sphere');
    count_sph=count_sph+1;
    model.component('comp1').geom('geom1').feature(tmptag_spo_os).set('r', 'SPO_MS_radius-SPO_dis_chl2wall');
    tmpstr_pos={num2str(tmp_center_X), num2str(tmp_center_Y), num2str(tmp_center_Z)};
    model.component('comp1').geom('geom1').feature(tmptag_spo_os).set('pos', tmpstr_pos);
    GLB_tag_SPO_CHL_os{count_SPO}=tmptag_spo_os;
    
    tmptag_spo_is=['sph',num2str(count_sph)];
    model.component('comp1').geom('geom1').create(tmptag_spo_is, 'Sphere');
    count_sph=count_sph+1;
    model.component('comp1').geom('geom1').feature(tmptag_spo_is).set('r', 'SPO_MS_radius-SPO_dis_chl2wall-SPO_thick_chl_layer');
    tmpstr_pos={num2str(tmp_center_X), num2str(tmp_center_Y), num2str(tmp_center_Z)};
    model.component('comp1').geom('geom1').feature(tmptag_spo_is).set('pos', tmpstr_pos);
    GLB_tag_SPO_CHL_is{count_SPO}=tmptag_spo_is;
    
    tmptag_spo_vac=['sph',num2str(count_sph)];
    model.component('comp1').geom('geom1').create(tmptag_spo_vac, 'Sphere');
    count_sph=count_sph+1;
    model.component('comp1').geom('geom1').feature(tmptag_spo_vac).set('r', 'SPO_MS_radius-SPO_dis_chl2wall-SPO_thick_chl_layer-SPO_MS_radius*(SPO_dis_mit2chl_rt+SPO_dis_vac2mit_rt)');
    tmpstr_pos={num2str(tmp_center_X), num2str(tmp_center_Y), num2str(tmp_center_Z)};
    model.component('comp1').geom('geom1').feature(tmptag_spo_vac).set('pos', tmpstr_pos);
    GLB_tag_SPO_VAC{count_SPO}=tmptag_spo_vac;
    
    count_SPO=count_SPO+1;
end

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

%%%process each cell with information from mat_contact_cell_idx
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
            %%%outer surface
            %%%shift wp toward/backward by distance SPO_dis_chl2wall
            tmp_pt1_a_o=tmp_pt1+tmp_vec1/norm(tmp_vec1)*SPO_dis_chl2wall;
            tmp_pt1_b_o=tmp_pt1-tmp_vec1/norm(tmp_vec1)*SPO_dis_chl2wall;
            tmp_pt2_a_o=tmp_pt1_a_o+tmp_vec2;
            tmp_pt2_b_o=tmp_pt1_b_o+tmp_vec2;
            tmp_pt3_a_o=tmp_pt1_a_o+tmp_vec3;
            tmp_pt3_b_o=tmp_pt1_b_o+tmp_vec3;
            %%%inner surface
            %%%shift wp toward/backward by distance SPO_dis_chl2wall
            tmp_pt1_a_i=tmp_pt1+tmp_vec1/norm(tmp_vec1)*(SPO_dis_chl2wall+SPO_thick_chl_layer);
            tmp_pt1_b_i=tmp_pt1-tmp_vec1/norm(tmp_vec1)*(SPO_dis_chl2wall+SPO_thick_chl_layer);
            tmp_pt2_a_i=tmp_pt1_a_i+tmp_vec2;
            tmp_pt2_b_i=tmp_pt1_b_i+tmp_vec2;
            tmp_pt3_a_i=tmp_pt1_a_i+tmp_vec3;
            tmp_pt3_b_i=tmp_pt1_b_i+tmp_vec3;
            %%%vac surface
            %%%shift wp toward/backward by distance SPO_dis_chl2wall
            tmp_pt1_a_vac=tmp_pt1+tmp_vec1/norm(tmp_vec1)*(SPO_dis_chl2wall+SPO_thick_chl_layer+SPO_MS_radius*(SPO_dis_mit2chl_rt+SPO_dis_vac2mit_rt));
            tmp_pt1_b_vac=tmp_pt1-tmp_vec1/norm(tmp_vec1)*(SPO_dis_chl2wall+SPO_thick_chl_layer+SPO_MS_radius*(SPO_dis_mit2chl_rt+SPO_dis_vac2mit_rt));
            tmp_pt2_a_vac=tmp_pt1_a_vac+tmp_vec2;
            tmp_pt2_b_vac=tmp_pt1_b_vac+tmp_vec2;
            tmp_pt3_a_vac=tmp_pt1_a_vac+tmp_vec3;
            tmp_pt3_b_vac=tmp_pt1_b_vac+tmp_vec3;
            
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
            
            %%% (2/4)process CHL outer surface
            tmptag_wp_a=['wp', num2str(count_wp)]; count_wp=count_wp+1;
            model.component('comp1').geom('geom1').create(tmptag_wp_a, 'WorkPlane');
            model.component('comp1').geom('geom1').feature(tmptag_wp_a).set('planetype', 'coordinates');
            model.component('comp1').geom('geom1').feature(tmptag_wp_a).set('genpoints', [tmp_pt1_a_o; tmp_pt2_a_o; tmp_pt3_a_o]);
            model.component('comp1').geom('geom1').feature(tmptag_wp_a).set('unite', true);
            tmptag_wp_b=['wp', num2str(count_wp)]; count_wp=count_wp+1;
            model.component('comp1').geom('geom1').create(tmptag_wp_b, 'WorkPlane');
            model.component('comp1').geom('geom1').feature(tmptag_wp_b).set('planetype', 'coordinates');
            model.component('comp1').geom('geom1').feature(tmptag_wp_b).set('genpoints', [tmp_pt1_b_o; tmp_pt2_b_o; tmp_pt3_b_o]);
            model.component('comp1').geom('geom1').feature(tmptag_wp_b).set('unite', true);
            %partition both obj with workplane
            %target = obj loop_i
            tmptag_obj=GLB_tag_SPO_CHL_os{loop_i};%{} since it's a cell array
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
                    GLB_tag_SPO_CHL_os(loop_i)=cellstr(tmptag_spl_dom2);
                else
                    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(cellstr(tmptag_spl_dom2));
                    GLB_tag_SPO_CHL_os(loop_i)=cellstr(tmptag_spl_dom1);
                end
            else
                GLB_tag_SPO_CHL_os(loop_i)=cellstr(tmptag_spl);
            end
            %target = obj loop_j
            tmptag_obj=GLB_tag_SPO_CHL_os{loop_j};%{} since it's a cell array
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
                    GLB_tag_SPO_CHL_os(loop_j)=cellstr(tmptag_spl_dom2);
                else
                    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(tmptag_spl_dom2);
                    GLB_tag_SPO_CHL_os(loop_j)=cellstr(tmptag_spl_dom1);
                end
            else
                GLB_tag_SPO_CHL_os(loop_j)=cellstr(tmptag_spl);
            end
            
            %%% (3/4)process CHL inner surface
            tmptag_wp_a=['wp', num2str(count_wp)]; count_wp=count_wp+1;
            model.component('comp1').geom('geom1').create(tmptag_wp_a, 'WorkPlane');
            model.component('comp1').geom('geom1').feature(tmptag_wp_a).set('planetype', 'coordinates');
            model.component('comp1').geom('geom1').feature(tmptag_wp_a).set('genpoints', [tmp_pt1_a_i; tmp_pt2_a_i; tmp_pt3_a_i]);
            model.component('comp1').geom('geom1').feature(tmptag_wp_a).set('unite', true);
            tmptag_wp_b=['wp', num2str(count_wp)]; count_wp=count_wp+1;
            model.component('comp1').geom('geom1').create(tmptag_wp_b, 'WorkPlane');
            model.component('comp1').geom('geom1').feature(tmptag_wp_b).set('planetype', 'coordinates');
            model.component('comp1').geom('geom1').feature(tmptag_wp_b).set('genpoints', [tmp_pt1_b_i; tmp_pt2_b_i; tmp_pt3_b_i]);
            model.component('comp1').geom('geom1').feature(tmptag_wp_b).set('unite', true);
            %partition both obj with workplane
            %target = obj loop_i
            tmptag_obj=GLB_tag_SPO_CHL_is{loop_i};%{} since it's a cell array
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
                    GLB_tag_SPO_CHL_is(loop_i)=cellstr(tmptag_spl_dom2);
                else
                    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(cellstr(tmptag_spl_dom2));
                    GLB_tag_SPO_CHL_is(loop_i)=cellstr(tmptag_spl_dom1);
                end
            else
                GLB_tag_SPO_CHL_is(loop_i)=cellstr(tmptag_spl);
            end
            %target = obj loop_j
            tmptag_obj=GLB_tag_SPO_CHL_is{loop_j};%{} since it's a cell array
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
                    GLB_tag_SPO_CHL_is(loop_j)=cellstr(tmptag_spl_dom2);
                else
                    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(tmptag_spl_dom2);
                    GLB_tag_SPO_CHL_is(loop_j)=cellstr(tmptag_spl_dom1);
                end
            else
                GLB_tag_SPO_CHL_is(loop_j)=cellstr(tmptag_spl);
            end
            
            %%% (4/4)process VAC surface
            tmptag_wp_a=['wp', num2str(count_wp)]; count_wp=count_wp+1;
            model.component('comp1').geom('geom1').create(tmptag_wp_a, 'WorkPlane');
            model.component('comp1').geom('geom1').feature(tmptag_wp_a).set('planetype', 'coordinates');
            model.component('comp1').geom('geom1').feature(tmptag_wp_a).set('genpoints', [tmp_pt1_a_vac; tmp_pt2_a_vac; tmp_pt3_a_vac]);
            model.component('comp1').geom('geom1').feature(tmptag_wp_a).set('unite', true);
            tmptag_wp_b=['wp', num2str(count_wp)]; count_wp=count_wp+1;
            model.component('comp1').geom('geom1').create(tmptag_wp_b, 'WorkPlane');
            model.component('comp1').geom('geom1').feature(tmptag_wp_b).set('planetype', 'coordinates');
            model.component('comp1').geom('geom1').feature(tmptag_wp_b).set('genpoints', [tmp_pt1_b_vac; tmp_pt2_b_vac; tmp_pt3_b_vac]);
            model.component('comp1').geom('geom1').feature(tmptag_wp_b).set('unite', true);
            %partition both obj with workplane
            %target = obj loop_i
            tmptag_obj=GLB_tag_SPO_VAC{loop_i};%{} since it's a cell array
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
                    GLB_tag_SPO_VAC(loop_i)=cellstr(tmptag_spl_dom2);
                else
                    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(cellstr(tmptag_spl_dom2));
                    GLB_tag_SPO_VAC(loop_i)=cellstr(tmptag_spl_dom1);
                end
            else
                GLB_tag_SPO_VAC(loop_i)=cellstr(tmptag_spl);
            end
            %target = obj loop_j
            tmptag_obj=GLB_tag_SPO_VAC{loop_j};%{} since it's a cell array
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
                    GLB_tag_SPO_VAC(loop_j)=cellstr(tmptag_spl_dom2);
                else
                    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(tmptag_spl_dom2);
                    GLB_tag_SPO_VAC(loop_j)=cellstr(tmptag_spl_dom1);
                end
            else
                GLB_tag_SPO_VAC(loop_j)=cellstr(tmptag_spl);
            end
            
        end
    end
    
    %%%create SPO CHL membrane
    tmptag_dif=['dif',num2str(count_dif)];count_dif=count_dif+1;
    model.component('comp1').geom('geom1').create(tmptag_dif, 'Difference');
    model.component('comp1').geom('geom1').feature(tmptag_dif).set('intbnd', false);
    model.component('comp1').geom('geom1').feature(tmptag_dif).selection('input').set(GLB_tag_SPO_CHL_os(loop_i));
    model.component('comp1').geom('geom1').feature(tmptag_dif).selection('input2').set(GLB_tag_SPO_CHL_is(loop_i));
    GLB_tag_SPO_CHL_mem{loop_i}=tmptag_dif;
    
    %contact with box?
    %if mat_contact_cell_idx(loop_i,loop_i)~=0    
    %end
end

%% SPO CHL
load('tmp_SPO_CHL_distribute.mat', 'SPO_CHL_center_final_pts', 'SPO_CHL_sph_selected_idx', 'SPO_CHL_sph_notselected_idx')
%SPO_CHL_center_final_pts_selected=SPO_CHL_center_final_pts(SPO_CHL_sph_selected_idx,:);
GLB_tag_SPO_CHL={};%TYPE=cell array of cell array
GLB_mat_contact_SPO_CHL_idx={};%TYPE=cell array of tmp_mat_contact_SPO_CHL_idx
disp('[Msg] start process SPO CHL')
for loop_i=1:Num_spo
    tmp_pts=SPO_CHL_center_final_pts{loop_i}(SPO_CHL_sph_selected_idx{loop_i},:);
    %%%%create CHL sphere
    tmptag_set_SPO_CHL={};
    for loop_j=1:size(tmp_pts,1)
        tmptag_sph=['sph',num2str(count_sph)];count_sph=count_sph+1;
        model.component('comp1').geom('geom1').create(tmptag_sph, 'Sphere');
        model.component('comp1').geom('geom1').feature(tmptag_sph).set('r', 'SPO_radius_chl_sphere');
        model.component('comp1').geom('geom1').feature(tmptag_sph).set('pos', string(tmp_pts(loop_j,:)));%{'5e-6' '0' '150e-6'}
        tmptag_set_SPO_CHL{loop_j}=tmptag_sph;
    end
    model.component('comp1').geom('geom1').run(tmptag_sph);

    %%%% chop chloroplast sphere
    %%%%%% store contact information into a matrix mat_contact_cell_idx;
    %%%%%% if ith cell overlaps jth cell, then mat_contact_cell_idx(i,j)=1
    %%%%%% if ith cell overlaps with box, then mat_contact_cell_idx(i,i)=1
    %%%%%% mat_contact_cell_idx is symetric
    tmp_mat_contact_SPO_CHL_idx=zeros(size(tmp_pts,1));% initial
    tmp_radius_2=(SPO_radius_chl_sphere)^2;
    for loop_i_CHL=1:size(tmp_pts,1)
        for loop_j_CHL=loop_i_CHL:size(tmp_pts,1)
            if loop_j_CHL==loop_i_CHL
                % contact with box?
                if (tmp_pts(loop_i_CHL,1)<SPO_radius_chl_sphere) ...,
                        || (tmp_pts(loop_i_CHL,1)>SPO_box_length-(SPO_radius_chl_sphere)) ...,
                        || (tmp_pts(loop_i_CHL,2)<SPO_radius_chl_sphere) ...,
                        || (tmp_pts(loop_i_CHL,2)>SPO_box_length-(SPO_radius_chl_sphere)) ...,
                        || (tmp_pts(loop_i_CHL,3)<SPO_radius_chl_sphere) ...,
                        || (tmp_pts(loop_i_CHL,3)>SPO_box_thick-(SPO_radius_chl_sphere))
                    tmp_mat_contact_SPO_CHL_idx(loop_i_CHL,loop_i_CHL)=1;
                end
            else
                if sum((tmp_pts(loop_i_CHL,:)-tmp_pts(loop_j_CHL,:)).^2) < tmp_radius_2*4 % (2*radius)^2
                    tmp_mat_contact_SPO_CHL_idx(loop_i_CHL,loop_j_CHL)=1;
                    tmp_mat_contact_SPO_CHL_idx(loop_j_CHL,loop_i_CHL)=1;
                end
            end
        end
    end
    %save box contact information for later process (Line???)
    GLB_mat_contact_SPO_CHL_idx{loop_i}=tmp_mat_contact_SPO_CHL_idx;
    
    %%%%%% process each CHL sphere with information from tmp_mat_contact_SPO_CHL_idx
    for loop_i_CHL=1:size(tmp_pts,1)
        tmp_idx_contact=find(tmp_mat_contact_SPO_CHL_idx(loop_i_CHL,:)==1);
        for tmp_loop_j_CHL=1:numel(tmp_idx_contact)
            loop_j_CHL=tmp_idx_contact(tmp_loop_j_CHL);
            if loop_j_CHL > loop_i_CHL
                tmp_center_a=tmp_pts(loop_i_CHL,:);
                tmp_center_b=tmp_pts(loop_j_CHL,:);
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
                
                tmptag_wp=['wp', num2str(count_wp)]; count_wp=count_wp+1;
                model.component('comp1').geom('geom1').create(tmptag_wp, 'WorkPlane');
                model.component('comp1').geom('geom1').feature(tmptag_wp).set('planetype', 'coordinates');
                model.component('comp1').geom('geom1').feature(tmptag_wp).set('genpoints', [tmp_pt1; tmp_pt2; tmp_pt3]);
                model.component('comp1').geom('geom1').feature(tmptag_wp).set('unite', true);
                %partition both obj with workplane
                %target = obj loop_i
                tmptag_obj=tmptag_set_SPO_CHL{loop_i_CHL};%{} since it's a cell array
                tmptag_par=['par', num2str(count_par)];count_par=count_par+1;
                model.component('comp1').geom('geom1').create(tmptag_par, 'Partition');
                model.component('comp1').geom('geom1').feature(tmptag_par).set('partitionwith', 'workplane');
                model.component('comp1').geom('geom1').feature(tmptag_par).set('workplane', tmptag_wp);
                model.component('comp1').geom('geom1').feature(tmptag_par).selection('input').set(tmptag_obj);
                %split par
                tmptag_spl=['spl', num2str(count_spl)];count_spl=count_spl+1;
                model.component('comp1').geom('geom1').create(tmptag_spl, 'Split');
                model.component('comp1').geom('geom1').feature(tmptag_spl).selection('input').set({tmptag_par});
                %[check]% rare case split still result 1 domain
                model.component('comp1').geom('geom1').run(tmptag_spl);
                if numel(model.geom('geom1').feature(tmptag_spl).objectNames())==1
                    tmptag_set_SPO_CHL(loop_i_CHL)=cellstr(tmptag_spl);
                else
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
                        tmptag_set_SPO_CHL(loop_i_CHL)=cellstr(tmptag_spl_dom2);
                    else
                        model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(tmptag_spl_dom2);
                        tmptag_set_SPO_CHL(loop_i_CHL)=cellstr(tmptag_spl_dom1);
                    end
                end
                
                %target = obj loop_j
                tmptag_obj=tmptag_set_SPO_CHL{loop_j_CHL};%{} since it's a cell array
                tmptag_par=['par', num2str(count_par)];count_par=count_par+1;
                model.component('comp1').geom('geom1').create(tmptag_par, 'Partition');
                model.component('comp1').geom('geom1').feature(tmptag_par).set('partitionwith', 'workplane');
                model.component('comp1').geom('geom1').feature(tmptag_par).set('workplane', tmptag_wp);
                model.component('comp1').geom('geom1').feature(tmptag_par).selection('input').set(tmptag_obj);
                %split par
                tmptag_spl=['spl', num2str(count_spl)];count_spl=count_spl+1;
                model.component('comp1').geom('geom1').create(tmptag_spl, 'Split');
                model.component('comp1').geom('geom1').feature(tmptag_spl).selection('input').set({tmptag_par});
                %[check]% rare case split still result 1 domain
                model.component('comp1').geom('geom1').run(tmptag_spl);
                if numel(model.geom('geom1').feature(tmptag_spl).objectNames())==1
                    tmptag_set_SPO_CHL(loop_j_CHL)=cellstr(tmptag_spl);
                else
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
                        tmptag_set_SPO_CHL(loop_j_CHL)=cellstr(tmptag_spl_dom2);
                    else
                        model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(tmptag_spl_dom2);
                        tmptag_set_SPO_CHL(loop_j_CHL)=cellstr(tmptag_spl_dom1);
                    end
                end
            end
        end

        %contact with box?
        if tmp_mat_contact_SPO_CHL_idx(loop_i_CHL,loop_i_CHL)~=0
            
        end
    end
    %%%% intersect each chl sphere with chl membrane
    TMP_tag_set_final_del={};%tmp_count=1;
    tmp_CHL_sph_selected_idx=SPO_CHL_sph_selected_idx{loop_i};
    for loop_i_CHL=1:numel(tmp_CHL_sph_selected_idx)
        tmptag_obj=tmptag_set_SPO_CHL{loop_i_CHL};%{} since it's a cell array
        tmptag_int=['int',num2str(count_int)];count_int=count_int+1;
        model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
        model.component('comp1').geom('geom1').feature(tmptag_int).set('keep', true);
        model.component('comp1').geom('geom1').feature(tmptag_int).set('intbnd', false);
        model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set([GLB_tag_SPO_CHL_mem(loop_i), cellstr(tmptag_obj)]);
        tmptag_set_SPO_CHL{loop_i_CHL}=tmptag_int;
        TMP_tag_set_final_del{loop_i_CHL}=tmptag_obj;
    end
    GLB_tag_SPO_CHL{loop_i}=tmptag_set_SPO_CHL;
    
    tmptag_del=['del',num2str(count_del)];count_del=count_del+1;
    model.component('comp1').geom('geom1').create(tmptag_del, 'Delete');
    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').init;
    model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set([TMP_tag_set_final_del,GLB_tag_SPO_CHL_mem(loop_i)]);

    disp(['[Msg] finish process SPO CHL No.',num2str(loop_i),'. In total ',num2str(Num_spo)])
end

%%%% count chloroplast number
%%%% VAR output: GLB_count_PAL_chl GLB_count_SPO_chl
%%%% these two vars are used in e_RT_geo_export_v0_1.m and
%%%% e_geo_AVRcal_v0_1.m
GLB_count_PAL_chl=[];
GLB_count_SPO_chl=[];
for loop_i=1:Num_pal%numel(GLB_tag_PAL_MS)
    tmptag_set_PAL_CHL=GLB_tag_PAL_CHL{loop_i};
    GLB_count_PAL_chl(loop_i)=numel(tmptag_set_PAL_CHL);
end
for loop_i=1:Num_spo%numel(GLB_tag_SPO_MS)
    tmptag_set_SPO_CHL=GLB_tag_SPO_CHL{loop_i};
    GLB_count_SPO_chl(loop_i)=numel(tmptag_set_SPO_CHL);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE module 1/2
%% save mph and mat for Ray Tracing (e_RT_geo_export)
mphsave(model,'tmpmph4RT_geo_export.mph');
save tmpsave4RT_geo_export.mat -regexp '^(?!(model|ans)$).'
disp("[Msg] exporting mat and mph for e_RT_geo_export module")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% start from this point
%% GLB_tag_* -> GLB_RT_tag_* and GLB_tag_*
%% [NOT revised]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PAL IAS
disp("[Msg] start to dif objects to generate domains")
tmptag_dif=['dif',num2str(count_dif)];count_dif=count_dif+1;
model.component('comp1').geom('geom1').create(tmptag_dif, 'Difference');
model.component('comp1').geom('geom1').feature(tmptag_dif).set('intbnd', false);
model.component('comp1').geom('geom1').feature(tmptag_dif).set('keep', true);
model.component('comp1').geom('geom1').feature(tmptag_dif).selection('input').set(GLB_tag_PAL_blk);%set({'blk1'})
model.component('comp1').geom('geom1').feature(tmptag_dif).selection('input2').set(GLB_tag_PAL_MS);

model.component('comp1').geom('geom1').run(tmptag_dif);%% IMPORTANT: run before measure
model.component('comp1').geom('geom1').measure.selection.init(3);
model.component('comp1').geom('geom1').measure.selection.set(tmptag_dif, 1);
try
    tmp_vol_IAS_PAL = model.geom('geom1').measure().getVolume()/(PAL_box_thick*PAL_box_length*PAL_box_length)
catch
    disp('[Warning][e_geo_main]: Fail to measure volume of PAL IAS.');
end

%% SPO IAS
tmptag_dif=['dif',num2str(count_dif)];count_dif=count_dif+1;
model.component('comp1').geom('geom1').create(tmptag_dif, 'Difference');
%count_dif=count_dif+1;
model.component('comp1').geom('geom1').feature(tmptag_dif).selection('input').set(GLB_tag_SPO_blk);
model.component('comp1').geom('geom1').feature(tmptag_dif).selection('input2').set(GLB_tag_SPO_MS);
model.component('comp1').geom('geom1').feature(tmptag_dif).set('intbnd', false);
model.component('comp1').geom('geom1').feature(tmptag_dif).set('keep', true);
%%calculate tmp_vol_IAS_SPO
model.component('comp1').geom('geom1').run(tmptag_dif);
model.component('comp1').geom('geom1').measure.selection.init(3);
model.component('comp1').geom('geom1').measure.selection.set(tmptag_dif, 1);
try
    tmp_vol_IAS_SPO = model.geom('geom1').measure().getVolume()/(SPO_box_thick*SPO_box_length*SPO_box_length)
catch
    disp('[Warning][e_geo_main]: Fail to measure volume of SPO IAS.');
end

%% PAL mitochondria
load('tmp_PAL_MIT_distribute.mat', 'PAL_MIT_center_final_pts')
GLB_tag_PAL_MIT={};%TYPE=cell array of cell array -- each cell is a cell array recording MIT tags for each PAL
GLB_mat_contact_PAL_MIT_idx={};%TYPE=cell array of tmp_mat_contact_PAL_MIT_idx
for loop_i=1:Num_pal
    tmp_pts=PAL_MIT_center_final_pts{loop_i};
    %%%%create MIT sphere
    tmptag_set_PAL_MIT={};
    %%%%%% store contact information into a matrix mat_contact_cell_idx;
    tmp_mat_contact_PAL_MIT_idx=zeros(1,size(tmp_pts,1));% initial
    for loop_j=1:size(tmp_pts,1)
        tmptag_sph=['sph',num2str(count_sph)];count_sph=count_sph+1;
        model.component('comp1').geom('geom1').create(tmptag_sph, 'Sphere');
        model.component('comp1').geom('geom1').feature(tmptag_sph).set('r', 'PAL_MS_radius*PAL_radius_mit_rt');
        model.component('comp1').geom('geom1').feature(tmptag_sph).set('pos', string(tmp_pts(loop_j,:)));%{'5e-6' '0' '150e-6'}
        tmptag_set_PAL_MIT{loop_j}=tmptag_sph;
        
        if GLB_mat_contact_PAL_MSVAC_idx(loop_i)==1
            if (tmp_pts(loop_j,1)<PAL_MS_radius*PAL_radius_mit_rt) ...,
                    || (tmp_pts(loop_j,1)>SPO_box_length-(PAL_MS_radius*PAL_radius_mit_rt)) ...,
                    || (tmp_pts(loop_j,2)<PAL_MS_radius*PAL_radius_mit_rt) ...,
                    || (tmp_pts(loop_j,2)>SPO_box_length-(PAL_MS_radius*PAL_radius_mit_rt))
                tmp_mat_contact_PAL_MIT_idx(loop_j)=1;
            end
        end
    end
    GLB_tag_PAL_MIT{loop_i}=tmptag_set_PAL_MIT; 
    GLB_mat_contact_PAL_MIT_idx{loop_i}=tmp_mat_contact_PAL_MIT_idx;
end
disp("[Msg] finish process PAL MIT")

%% SPO mitochondria
load('tmp_SPO_MIT_distribute.mat', 'SPO_MIT_center_final_pts')
GLB_tag_SPO_MIT={};%TYPE=cell array of cell array -- each cell is a cell array recording MIT tags for each SPO
GLB_mat_contact_SPO_MIT_idx={};%TYPE=cell array of tmp_mat_contact_SPO_MIT_idx
for loop_i=1:numel(GLB_tag_SPO_MS)
    tmp_pts=SPO_MIT_center_final_pts{loop_i};
    %%%%create MIT sphere
    tmptag_set_SPO_MIT={};
    %%%%%% store contact information into a matrix mat_contact_cell_idx;
    tmp_mat_contact_SPO_MIT_idx=zeros(1,size(tmp_pts,1));% initial
    for loop_j=1:size(tmp_pts,1)
        tmptag_sph=['sph',num2str(count_sph)];count_sph=count_sph+1;
        model.component('comp1').geom('geom1').create(tmptag_sph, 'Sphere');
        model.component('comp1').geom('geom1').feature(tmptag_sph).set('r', 'SPO_MS_radius*SPO_radius_mit_rt');
        model.component('comp1').geom('geom1').feature(tmptag_sph).set('pos', string(tmp_pts(loop_j,:)));%{'5e-6' '0' '150e-6'}
        tmptag_set_SPO_MIT{loop_j}=tmptag_sph;
        
        if GLB_mat_contact_SPO_MSVAC_idx(loop_i)==1
            if (tmp_pts(loop_j,1)<SPO_MS_radius*SPO_radius_mit_rt) ...,
                    || (tmp_pts(loop_j,1)>SPO_box_length-(SPO_MS_radius*SPO_radius_mit_rt)) ...,
                    || (tmp_pts(loop_j,2)<SPO_MS_radius*SPO_radius_mit_rt) ...,
                    || (tmp_pts(loop_j,2)>SPO_box_length-(SPO_MS_radius*SPO_radius_mit_rt))
                tmp_mat_contact_SPO_MIT_idx(loop_j)=1;
            end
        end
    end
    GLB_tag_SPO_MIT{loop_i}=tmptag_set_SPO_MIT; 
    GLB_mat_contact_SPO_MIT_idx{loop_i}=tmp_mat_contact_SPO_MIT_idx;
end
disp("[Msg] finish process SPO MIT")

%% PAL -- Trim domain; contribute to selection
GLB_tag_PAL_CYT={};
TMP_tag_set_final_del={};tmp_count=1;
for loop_i=1:numel(GLB_tag_PAL_MS)
    %%% diff to obtain CYT
    tmptag_dif=['dif',num2str(count_dif)];count_dif=count_dif+1;
    model.component('comp1').geom('geom1').create(tmptag_dif, 'Difference');
    model.component('comp1').geom('geom1').feature(tmptag_dif).set('intbnd', false);
    model.component('comp1').geom('geom1').feature(tmptag_dif).set('keep', true);
    model.component('comp1').geom('geom1').feature(tmptag_dif).selection('input').set(GLB_tag_PAL_MS(loop_i));%set({'blk1'})
    model.component('comp1').geom('geom1').feature(tmptag_dif).selection('input2').set([GLB_tag_PAL_CHL{loop_i}, GLB_tag_PAL_MIT{loop_i}, GLB_tag_PAL_VAC(loop_i)]);
    GLB_tag_PAL_CYT{loop_i}=tmptag_dif;
    TMP_tag_set_final_del{tmp_count}=GLB_tag_PAL_MS{loop_i};tmp_count=tmp_count+1;
    
    %%% Trim CYT CHL MIT VAC with box
    %%% based on box contact information
    %%% GLB_mat_contact_PAL_MSVAC_idx -- array
    %%%     -- TAG -- GLB_tag_PAL_CYT
    %%%     -- TAG -- GLB_tag_PAL_VAC
    %%% GLB_mat_contact_PAL_CHL_idx -- cell array of matrix
    %%%     -- TAG -- GLB_tag_PAL_CHL
    %%% GLB_mat_contact_PAL_MIT_idx -- cell array of matrix
    %%%     -- TAG -- GLB_tag_PAL_MIT
    
    if GLB_mat_contact_PAL_MSVAC_idx(loop_i)==1
        %%%%% trim CYT
        tmptag_int=['int',num2str(count_int)];count_int=count_int+1;
        model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
        model.component('comp1').geom('geom1').feature(tmptag_int).set('keep', true);
        model.component('comp1').geom('geom1').feature(tmptag_int).set('intbnd', false);
        model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set([cellstr(GLB_tag_PAL_blk), cellstr(GLB_tag_PAL_CYT{loop_i})]);
        TMP_tag_set_final_del{tmp_count}=GLB_tag_PAL_CYT{loop_i};tmp_count=tmp_count+1;
        GLB_tag_PAL_CYT{loop_i}=tmptag_int;
        %%%%% trim VAC
        tmptag_int=['int',num2str(count_int)];count_int=count_int+1;
        model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
        model.component('comp1').geom('geom1').feature(tmptag_int).set('keep', true);
        model.component('comp1').geom('geom1').feature(tmptag_int).set('intbnd', false);
        model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set([cellstr(GLB_tag_PAL_blk), cellstr(GLB_tag_PAL_VAC{loop_i})]);
        TMP_tag_set_final_del{tmp_count}=GLB_tag_PAL_VAC{loop_i};tmp_count=tmp_count+1;
        GLB_tag_PAL_VAC{loop_i}=tmptag_int;
        
        %%%%% trim CHL; for loop
        tmptag_set_PAL_CHL=GLB_tag_PAL_CHL{loop_i};
        tmp_mat_contact_PAL_CHL_idx=GLB_mat_contact_PAL_CHL_idx{loop_i};
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
        
        %%%%% trim MIT; for loop
        tmptag_set_PAL_MIT=GLB_tag_PAL_MIT{loop_i};
        tmp_mat_contact_PAL_MIT_idx=GLB_mat_contact_PAL_MIT_idx{loop_i};
        for loop_i_MIT=1:numel(tmptag_set_PAL_MIT)
            if tmp_mat_contact_PAL_MIT_idx(loop_i_MIT)==1
                tmptag_int=['int',num2str(count_int)];count_int=count_int+1;
                model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
                model.component('comp1').geom('geom1').feature(tmptag_int).set('keep', true);
                model.component('comp1').geom('geom1').feature(tmptag_int).set('intbnd', false);
                model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set([cellstr(GLB_tag_PAL_blk), cellstr(tmptag_set_PAL_MIT{loop_i_MIT})]);
                TMP_tag_set_final_del{tmp_count}=tmptag_set_PAL_MIT{loop_i_MIT};tmp_count=tmp_count+1;
                tmptag_set_PAL_MIT{loop_i_MIT}=tmptag_int;
            end
        end
        GLB_tag_PAL_MIT{loop_i}=tmptag_set_PAL_MIT;%update tmp to GLB
        
        %%%%% end trim
    end
    
end
%%% delete TMP_tag_set_final_del in each loop
tmptag_del=['del',num2str(count_del)];count_del=count_del+1;
model.component('comp1').geom('geom1').create(tmptag_del, 'Delete');
model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').init;
model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(TMP_tag_set_final_del);

%% SPO -- Trim domain; contribute to selection
GLB_tag_SPO_CYT={};
TMP_tag_set_final_del={};tmp_count=1;
for loop_i=1:numel(GLB_tag_SPO_MS)
    %%% diff to obtain CYT
    tmptag_dif=['dif',num2str(count_dif)];count_dif=count_dif+1;
    model.component('comp1').geom('geom1').create(tmptag_dif, 'Difference');
    model.component('comp1').geom('geom1').feature(tmptag_dif).set('intbnd', false);
    model.component('comp1').geom('geom1').feature(tmptag_dif).set('keep', true);
    model.component('comp1').geom('geom1').feature(tmptag_dif).selection('input').set(GLB_tag_SPO_MS(loop_i));%set({'blk1'})
    model.component('comp1').geom('geom1').feature(tmptag_dif).selection('input2').set([GLB_tag_SPO_CHL{loop_i}, GLB_tag_SPO_MIT{loop_i}, GLB_tag_SPO_VAC(loop_i)]);
    GLB_tag_SPO_CYT{loop_i}=tmptag_dif;
    TMP_tag_set_final_del{tmp_count}=GLB_tag_SPO_MS{loop_i};tmp_count=tmp_count+1;
    
    %%% Trim CYT CHL MIT VAC with box
    %%% based on box contact information
    %%% GLB_mat_contact_SPO_MSVAC_idx -- array
    %%%     -- TAG -- GLB_tag_SPO_CYT
    %%%     -- TAG -- GLB_tag_SPO_VAC
    %%% GLB_mat_contact_SPO_CHL_idx -- cell array of matrix
    %%%     -- TAG -- GLB_tag_SPO_CHL
    %%% GLB_mat_contact_SPO_MIT_idx -- cell array of matrix
    %%%     -- TAG -- GLB_tag_SPO_MIT
    
    if GLB_mat_contact_SPO_MSVAC_idx(loop_i)==1
        %%%%% trim CYT
        tmptag_int=['int',num2str(count_int)];count_int=count_int+1;
        model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
        model.component('comp1').geom('geom1').feature(tmptag_int).set('keep', true);
        model.component('comp1').geom('geom1').feature(tmptag_int).set('intbnd', false);
        model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set([cellstr(GLB_tag_SPO_blk), cellstr(GLB_tag_SPO_CYT{loop_i})]);
        TMP_tag_set_final_del{tmp_count}=GLB_tag_SPO_CYT{loop_i};tmp_count=tmp_count+1;
        GLB_tag_SPO_CYT{loop_i}=tmptag_int;
        %%%%% trim VAC
        tmptag_int=['int',num2str(count_int)];count_int=count_int+1;
        model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
        model.component('comp1').geom('geom1').feature(tmptag_int).set('keep', true);
        model.component('comp1').geom('geom1').feature(tmptag_int).set('intbnd', false);
        model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set([cellstr(GLB_tag_SPO_blk), cellstr(GLB_tag_SPO_VAC{loop_i})]);
        TMP_tag_set_final_del{tmp_count}=GLB_tag_SPO_VAC{loop_i};tmp_count=tmp_count+1;
        GLB_tag_SPO_VAC{loop_i}=tmptag_int;
        
        %%%%% trim CHL; for loop
        tmptag_set_SPO_CHL=GLB_tag_SPO_CHL{loop_i};
        tmp_mat_contact_SPO_CHL_idx=GLB_mat_contact_SPO_CHL_idx{loop_i};
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
        
        %%%%% trim MIT; for loop
        tmptag_set_SPO_MIT=GLB_tag_SPO_MIT{loop_i};
        tmp_mat_contact_SPO_MIT_idx=GLB_mat_contact_SPO_MIT_idx{loop_i};
        for loop_i_MIT=1:numel(tmptag_set_SPO_MIT)
            if tmp_mat_contact_SPO_MIT_idx(loop_i_MIT)==1
                tmptag_int=['int',num2str(count_int)];count_int=count_int+1;
                model.component('comp1').geom('geom1').create(tmptag_int, 'Intersection');
                model.component('comp1').geom('geom1').feature(tmptag_int).set('keep', true);
                model.component('comp1').geom('geom1').feature(tmptag_int).set('intbnd', false);
                model.component('comp1').geom('geom1').feature(tmptag_int).selection('input').set([cellstr(GLB_tag_SPO_blk), cellstr(tmptag_set_SPO_MIT{loop_i_MIT})]);
                TMP_tag_set_final_del{tmp_count}=tmptag_set_SPO_MIT{loop_i_MIT};tmp_count=tmp_count+1;
                tmptag_set_SPO_MIT{loop_i_MIT}=tmptag_int;
            end
        end
        GLB_tag_SPO_MIT{loop_i}=tmptag_set_SPO_MIT;%update tmp to GLB
        
        %%%%% end trim
    end
    
end
%%% delete TMP_tag_set_final_del in each loop
tmptag_del=['del',num2str(count_del)];count_del=count_del+1;
model.component('comp1').geom('geom1').create(tmptag_del, 'Delete');
model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').init;
model.component('comp1').geom('geom1').feature(tmptag_del).selection('input').set(TMP_tag_set_final_del);
model.geom('geom1').run(tmptag_del);

%% form assembly
% model.component('comp1').geom('geom1').feature('fin').label('Form Assembly');
% model.component('comp1').geom('geom1').feature('fin').set('action', 'assembly');
% model.component('comp1').geom('geom1').feature('fin').set('createpairs', false);
% model.component('comp1').geom('geom1').run;
% model.component('comp1').geom('geom1').run('fin');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SAVE module 2/2
%% save mph and mat for Reaction-Diffusion (RD)
save tmp_all.mat -regexp '^(?!(model|ans)$).'
mphsave(model, 'tmp_geo_t1.mph')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
