% eLeaf: 3D model of dicot leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 0.1

%% calculate AVR
%% prepare leaf dim and [chl] for ray tracing
%% modified from e_physics
clear all;
load parainput.mat

import com.comsol.model.*
import com.comsol.model.util.*

model=mphload('tmp_geomIP_mesh_nocresel.mph');
load SAVE_e_geom.mat

%% volume profile of each chloroplasts
tmp_count_chl=0;
for loop_cell=1:size(GLB_tag_PAL_CHL,2)
    for loop_chlo=1:size(GLB_tag_PAL_CHL{loop_cell},2)
        tmptag_dom1=GLB_tag_PAL_CHL{loop_cell}{loop_chlo};
        model.component('comp1').geom('geom1').measure.selection.init(3);
        model.component('comp1').geom('geom1').measure.selection.set(tmptag_dom1, 1);
        tmp_vol_dom1 = model.geom('geom1').measure().getVolume();
        
        tmp_count_chl=tmp_count_chl+1;
        chlo_vol_profile(tmp_count_chl,1)=tmp_vol_dom1;
    end
end
for loop_cell=1:size(GLB_tag_SPO_CHL,2)
    for loop_chlo=1:size(GLB_tag_SPO_CHL{loop_cell},2)
        tmptag_dom1=GLB_tag_SPO_CHL{loop_cell}{loop_chlo};
        model.component('comp1').geom('geom1').measure.selection.init(3);
        model.component('comp1').geom('geom1').measure.selection.set(tmptag_dom1, 1);
        tmp_vol_dom1 = model.geom('geom1').measure().getVolume();
        
        tmp_count_chl=tmp_count_chl+1;
        chlo_vol_profile(tmp_count_chl,1)=tmp_vol_dom1;
    end
end

% %%volume of all chloroplast
% tmp_sel=mphgetselection(model.selection('uni1'));
% tmp_bnd_num=tmp_sel.entities;
% model.geom('geom1').measureFinal.selection.geom('geom1', 3);
% model.geom('geom1').measureFinal.selection.set(tmp_bnd_num);
% all_chlo_vol=model.geom('geom1').measureFinal().getVolume();
xmax=PAL_box_length;
ymax=PAL_box_length;
zmax=SPO_box_thick+PAL_box_thick+EPL_u_thick+0.1e-6;
xmin=0;ymin=0;zmin=-EPL_l_thick-0.1e-6;
leaf_surface_area=(xmax-xmin)*(ymax-ymin);
chlo_vol_all=sum(chlo_vol_profile);
AVR=chlo_vol_all/leaf_surface_area;
%AVR=3e-5;%30mL/m2 in Zhu et al., 2007, i.e. 3e-5; for eLeaf_rice, ca 1.0~1.5e-5
chlo_vol_PAL=sum(chlo_vol_profile(1:sum(GLB_count_PAL_chl)));
chlo_vol_SPO=sum(chlo_vol_profile((sum(GLB_count_PAL_chl)+1):end));

%convert [chl] from mg m-2 to g/m3
% ([chl]*1e-3)/(AVR)
%chl_con=463.4425;%mg/m2
%chl_con_ratio_PvS=1.1;

% eqn1:
% (chl_PAL*chlo_vol_PAL+chl_SPO*chlo_vol_SPO)/leaf_surface_area=chl_con*1e3
chl_con_SPO_4RT=chl_con*1e-3*leaf_surface_area/(chl_con_ratio_PvS*chlo_vol_PAL+chlo_vol_SPO);
chl_con_PAL_4RT=chl_con_SPO_4RT*chl_con_ratio_PvS;

fid=fopen('count_chl4RT','w');
fprintf(fid,'%d %d\n', Num_pal, Num_spo);
fprintf(fid,'%d %e\n',[GLB_count_PAL_chl;chl_con_PAL_4RT*ones(1,Num_pal)]);
fprintf(fid,'%d %e\n',[GLB_count_SPO_chl;chl_con_SPO_4RT*ones(1,Num_spo)]);
fclose(fid);

save SAVE_e_geom4RT.mat -regexp '^(?!(model|ans)$).'
save chlo_vol_profile.mat chlo_vol_profile