% eLeaf: 3D model of dicot leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 0.1

function e_geo_addcresel_v0_1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update 2021 April
% - mph geometry pass the test of free tetrahedron mesh
% before export object for ray tracing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%path('C:\Program Files\COMSOL\COMSOL52a\Multiphysics\mli',path)
%mphstart(2036);
%import com.comsol.model.*
%import com.comsol.model.util.*
model=mphload('tmp_geomIP_mesh_nocresel.mph');
load list_cresel
for count_list_cresel=1:size(list_cresel,2)
    tmptag_cresel=list_cresel(count_list_cresel);
    %model.geom('geom1').feature(tmptag_cresel).set('createselection',true);%COMSOL5.3
    model.geom('geom1').feature(tmptag_cresel).set('selresult',true);%COMSOL5.3
end
model.geom('geom1').run;
mphsave(model,'tmp_geomIP_mesh_cresel.mph');
%ModelUtil.disconnect;
%exit;
end