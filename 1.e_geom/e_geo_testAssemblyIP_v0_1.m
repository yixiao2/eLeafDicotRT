% eLeaf: 3D model of dicot leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 0.1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update 2020/12/14
% - mph geometry pass the test of free tetrahedron mesh
% before export object for ray tracing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model.geom('geom1').feature('fin').name('Form Assembly');
model.geom('geom1').feature('fin').set('action', 'assembly');
model.geom('geom1').feature('fin').set('createpairs', 'on');
model.geom('geom1').feature('fin').set('imprint', 'on');
%model.geom('geom1').feature('fin').set('repairtol', '1e-5');
%model.geom('geom1').feature(tmptag_IAS).set('repairtoltype', 'auto');
model.geom('geom1').run;
mphsave(model,'tmp_geomIP_nocresel.mph')

%% test free tetrahedron
try
    %model.mesh.create('mesh1', 'geom1');
    %model.mesh('mesh1').feature.create('ftet1', 'FreeTet');
    %model.mesh('mesh1').run;
    
    % output tags for create selection
%     if loop_idx==-1&&loop_idy==-1
        save SAVE_e_geom.mat -regexp '^(?!(model|ans)$).'
%     end
    %save list_cresel.mat list_cresel
    %toc
catch
    display('FAIL free tetrahedron mesh');
    %pause;
    %ModelUtil.remove('Model_PAL_MS_distribute');
    %ModelUtil.remove('Model_PAL_CHL_distribute');
    %ModelUtil.remove('Model_PAL_MIT_distribute');
    %ModelUtil.remove('Model_SPO_MS_distribute');
    %ModelUtil.remove('Model_SPO_CHL_distribute');
    %ModelUtil.remove('Model_SPO_MIT_distribute');
    ModelUtil.remove('Model_3Dleaf_dicot')
%     e_geo_PAL_MS_distribute_v0_1;
%     e_geo_PAL_CHL_distribute_v0_1;
%     e_geo_PAL_MIT_distribute_v0_1;
%     e_geo_SPO_MS_distribute_v0_1;
%     e_geo_SPO_CHL_distribute_v0_1;
%     e_geo_SPO_MIT_distribute_v0_1;
%     e_geo_main_v0_1;
%     e_geo_test_v0_1;
end
mphsave(model,'tmp_geomIP_mesh_nocresel.mph')
%% end test