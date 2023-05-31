function e_geo_run_eleaf_dicot_v0_1(GLB_RUN_MODE, CFG_PARA_COM, CFG_VARARGIN)
%run_eleaf_dicot_v0_1
% demo: e_geo_run_eleaf_dicot_v0_1(1,1,ones(1,13))
%
%GLB_RUN_MODE=2;% [1 2 22 23 3 4 5]
%
%GLB_RT_MODE
% 1: local mode, parpool, and then parfor
% 2: cluster mode. write all bash commands to sh files
% 3: debug mode. test part of the ray tracing bash command with for-loop.
%
%CFG_PARA_COM=zeros(1,13);

save('TMP_CFG.mat','GLB_RUN_MODE','CFG_*');

rng shuffle

%GLB_RT_MODE=1; % 1=local parfor mode, 2=cluster mode, 3=debug mode
%if GLB_RUN_MODE==4
%    GLB_RT_MODE=CFG_VARARGIN(3);
%end

FLAG_RT_TRAIN=1; % 1=yes, 0=no
rep_num=1;

switch GLB_RUN_MODE
    case 1
        %% [1/?] 3D reconstruction
        FILE_NAME_DIARY=['RUNLOG_MODE1_',datestr(now, 'yymmddHHMM'),'.log'];
        diary(FILE_NAME_DIARY)
        tic
        e_geo_parainput_v0_1_b4fit(CFG_PARA_COM);toc;
        e_geo_PAL_MS_distribute_v0_1;toc;
        e_geo_PAL_CHL_distribute_v0_1;toc;
        e_geo_PAL_MIT_distribute_v0_1;toc;
        e_geo_SPO_MS_distribute_v0_1;toc;
        e_geo_SPO_CHL_distribute_v0_1_Sc;toc;
        e_geo_SPO_MIT_distribute_v0_1;toc;
        tic;e_geo_main_v0_1b;toc;%%SAVE_e_geom.mat
        %%in the middle, save tmpsave4RT_geo_export.mat
        e_geo_testAssemblyIP_v0_1;%%test form assembly with imprint
        %e_geo_addcresel_v0_1;
        e_geo_AVRcal_v0_1;%%SAVE_e_geom4RT.mat
        toc
        diary off
    case 2
        %% [2/?] Ray tracing
        tic
        %%%% [2-1] RT prepare
        if(exist('../2.e_raytracing','dir')~=7)
            disp('no dir 2.e_raytracing');
        else
            copyfile('tmpsave4RT_geo_export.mat','../2.e_raytracing/geo_export');
            copyfile('tmpmph4RT_geo_export.mph','../2.e_raytracing/geo_export');
            %copyfile('tmp_MS3D.mat','../2.e_raytracing/geo_export');
            %copyfile('tmpCK_geomIP_cad_mesh_cresel.mph','../2.e_raytracing/geo_export');
            %copyfile('tmpCK_geomIP_cad_mesh_cresel.mph','../2.5.eleaf_fvcb_fit');
            copyfile('SAVE_e_geom4RT.mat','../2.e_raytracing/geo_export');
            %copyfile('save_e_geom.mat','../2.5.eleaf_fvcb_fit');

            copyfile('count_chl4RT','../2.e_raytracing');%output from e_geo_AVR_cal_v0_1//e_RT_geo_export_v0_1.m
        end

        %generate leaf box (boundary of ray tracing)
        fileID = fopen('leaf','w');
        %fprintf(fileID,'8\n12\n');
        fprintf(fileID,'ply\nformat ascii 1.0\n');
        fprintf(fileID,'element vertex 8\n');
        fprintf(fileID,'property float x\n');
        fprintf(fileID,'property float y\n');
        fprintf(fileID,'property float z\n');
        fprintf(fileID,'element face 12\n');
        fprintf(fileID,'property list uchar int vertex_indices\nend_header\n');
        load('SAVE_e_geom4RT.mat','xmin','xmax','ymin','ymax','zmin','zmax');
        fprintf(fileID,'%e %e %e\n',xmin,ymax,zmin);
        fprintf(fileID,'%e %e %e\n',xmax,ymax,zmin);
        fprintf(fileID,'%e %e %e\n',xmax,ymin,zmin);
        fprintf(fileID,'%e %e %e\n',xmin,ymin,zmin);
        fprintf(fileID,'%e %e %e\n',xmin,ymax,zmax);
        fprintf(fileID,'%e %e %e\n',xmax,ymax,zmax);
        fprintf(fileID,'%e %e %e\n',xmax,ymin,zmax);
        fprintf(fileID,'%e %e %e\n',xmin,ymin,zmax);
        fprintf(fileID,'3 0 1 3\n3 1 2 3\n3 0 3 4\n3 3 4 7\n3 0 1 4\n3 1 4 5\n3 1 2 5\n3 2 5 6\n3 2 3 7\n3 2 6 7\n3 4 5 7\n3 5 6 7\n');
        fclose(fileID);
        copyfile('leaf','../2.e_raytracing');

        %%%% [2-2] RT geometry export and manually check
        %%%% run geo_export to generate triangle meshed surface
        cd ../2.e_raytracing/geo_export
        %%%% [2-2][1/7] Try export surface mesh with repair/delete_small_edges
        e_RT_geo_export_v0_1
        %%%% [2-2][2/7] Check exported ply files manually
        %e_RT_geo_export_vis_geo_manualcheck_v0_1
        %%%% [2-2][3/7] Save mph files for objs failed in 1&2
        %%%% GLB_savefailmph_mode=1
        e_RT_geo_export_savefailmph_v0_1

        cd ../../1.e_geom/
        toc

    case 22
        %% manual fix if automatic fix failed
        tic
        clearvars
        cd ../2.e_raytracing/geo_export
        %%pause();%%%% manually check mph files and try to export ply files
        %%%% [2-2][4/7][optional] GLB_nodeleterepair_mode=2 for some surfaces which cannot by exported
        %%%% use TMP_SAVE_RT_geo_export_aftervisualcheck_full.mat
        SAVE_FAIL_LIST_ext_EPL_l=[];
        SAVE_FAIL_LIST_ext_EPL_u=[];
        SAVE_FAIL_LIST_ext_PAL_MS=[];
        SAVE_FAIL_LIST_ext_PAL_VAC=[];
        SAVE_FAIL_LIST_ext_PAL_CHL=[4,10;4,30;4,35;4,77];
        SAVE_FAIL_LIST_ext_SPO_MS=[];
        SAVE_FAIL_LIST_ext_SPO_VAC=[];
        SAVE_FAIL_LIST_ext_SPO_CHL=[22,1];
        save TMP_SAVE_RT_geo_export_aftervisualcheck_full.mat SAVE_FAIL_LIST_ext_*
        %%%%%%e_RT_geo_export_savefailmph_v0_1

        %%%% [2-2][5/7] Check all ply files again, record any failed objs
        %%e_RT_geo_export_vis_geo_manualcheck_v0_1
        %%%% [2-2][6/7] Debug missing failed obj;
        %%%% GLB_savefailmph_mode=2
        e_RT_geo_export_savefailmph_v0_1
        %%%% [2-2][7/7] GLB_nodeleterepair_mode=2 for some cases
        %%%%%%e_RT_geo_export_savefailmph_v0_1

        cd ../../1.e_geom/
        toc
    case 3
        tic
        cd ../2.e_raytracing/geo_export
        % disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        % disp('eLeaf_dicot Phase I finished.')
        % disp('Please check mesh exported for RayTracing,')
        % disp('and continue Phase II.')
        % disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        % pause()

        %generate Defs.h in geo_export
        load SAVE_RT_geo_export.mat Num_pal Num_spo GLB_count_PAL_chl GLB_count_SPO_chl GLB_tag_*
        copyfile('Defs_template.h','Defs.h');
        fileID=fopen('Defs.h','a');
        load SAVE_e_geom4RT.mat xmax xmin ymax ymin zmax zmin GLB*
        fprintf(fileID,'#define xmax %e\n#define xmin %e\n#define zmax %e\n#define zmin %e\n#define ymax %e\n#define ymin %e\n\n',xmax,xmin,zmax,zmin,ymax,ymin);
        fprintf(fileID,'#define pal_num %d\n#define spo_num %d\n#define ms_num (pal_num+spo_num)\nint count_ms;\n',Num_pal,Num_spo);
        fprintf(fileID,'#define ms_max_chl_num %d\nint ms_chl_num[ms_num];\n',max([GLB_count_PAL_chl,GLB_count_SPO_chl]));
        fprintf(fileID,'double ms_chl_con[ms_num];\ndouble ms_chl_con_v2[ms_num][ms_max_chl_num];\n');
        fprintf(fileID,'int flag_format_count_chl;\n');
        fprintf(fileID,'#define epl_l_num %d\n#define epl_u_num %d\n#define nonms_num %d\n',numel(GLB_tag_EPL_l_set),numel(GLB_tag_EPL_u_set),numel(GLB_tag_EPL_l_set)+numel(GLB_tag_EPL_u_set));
        fprintf(fileID,'int count_nonms;\n\n');
        fprintf(fileID,'Object *p_cell_ms[ms_num];\nObject *p_chl_ms[ms_num][ms_max_chl_num];\nObject *p_vac_ms[ms_num];\nObject *p_leaf;\nObject *p_cell_ns[nonms_num];\n\n');
        fprintf(fileID,'#endif\n');
        fclose(fileID);
        copyfile('Defs.h','../');

        %copy ms*.ply to ../MS/
        copyfile('PAL*.ply','../MS/');
        copyfile('SPO*.ply','../MS/');
        %copy ns*.ply to ../nonMS/
        copyfile('EPL*.ply','../nonMS/');

        %% check inside-belong-mother_node relationship
        cd ..
        system('make');
        system('./inner_check count_chl4RT');%pause();%%check the output msg. If no problem, then continue.

        cd ../1.e_geom/
        toc
    case 4
        GLB_RT_MODE=CFG_VARARGIN(1);
        cd ../2.e_raytracing/
        if(exist('TMP','dir')~=7)
            mkdir('TMP');
        end
        %% ray tracing
        %%%% run 500*500 rays simulation by splitting into 100*100 tasks.
        %%%% use minimum of SAC_water & SAC_chl under RGB, then recalculate
        %%%% light absorption under RGB ligth
        %% [water, chl] SAC default = [7.5e-05*100, 4.94e+04*1e-4] blue 445nm
        %%                          = [0.000619*100, 1.22e+04*1e-4] green 560nm
        %%                          = [0.003108*100, 2.53e+04*1e-4] red 640nm
        %%                          = [0.000114*100, 4.26e+04*1e-4] blue 475nm
        %%                          = [0.002834*100, 2.34e+04*1e-4] red 625nm
        %% unit = [m-1, m2 g-1]
        %% exp for light ab profile measurement: blue=445nm; gree=561nm; red=638nm;
        %% exp for Licor red-blue measurement: blue=475nm; red=625nm;
        tic
        SAC_licor_blue475nm=[0.000114*100, 4.26e+04*1e-4];
        SAC_licor_red625nm=[0.002834*100, 2.34e+04*1e-4];

        if FLAG_RT_TRAIN==1
            % for profile training, tmp_SAC/20. Increase time complexity.
            tmp_SAC_water=min(SAC_licor_blue475nm(1),SAC_licor_red625nm(1))/20;
            tmp_SAC_chl=min(SAC_licor_blue475nm(2),SAC_licor_red625nm(2))/20;
        else
            tmp_SAC_water=min(SAC_licor_blue475nm(1),SAC_licor_red625nm(1));
            tmp_SAC_chl=min(SAC_licor_blue475nm(2),SAC_licor_red625nm(2));
        end

        %%%%%%%%% CASE5: cut x&y every 5 rays; batch number 100*100
        %rep_num=1;
        %CASE_NUM=5;
        RT_x_range=500;
        RT_y_range=500;
        RT_x_perthread=5;
        RT_y_perthread=5;
        num_loop_x=RT_x_range/RT_x_perthread;
        num_loop_y=RT_y_range/RT_y_perthread;
        count_batch=0;
        batch_cmd_all={};
        for loop_x=1:num_loop_x
            for loop_y=1:num_loop_y
                count_batch=count_batch+1;
                ray_x_start=RT_x_perthread*(loop_x-1);
                ray_x_end=RT_x_perthread*(loop_x);
                ray_y_start=RT_y_perthread*(loop_y-1);
                ray_y_end=RT_y_perthread*(loop_y);
                %%%% run RT using 1/20*chl_SAC, assuming 1/10 is the lower limit of
                %%%% [chl] and 1/2 is the lower limit of chl_SAC itself.
                %%%% 1.22e4*1e-4/20=0.0610
                tmp_bash_cmd=['./rt4abprofiletraining ',num2str(tmp_SAC_water,'%e'),' ',num2str(tmp_SAC_chl,'%e'),' ',...
                    'TMP/results_abevents_tmpnm_500x_rep',num2str(rep_num),'_',num2str(count_batch),' ',...
                    'count_chl4RT ',...
                    num2str(RT_x_range),' ',num2str(RT_y_range),' ',...
                    num2str(ray_x_start),' ',num2str(ray_x_end),' ',...
                    num2str(ray_y_start),' ',num2str(ray_y_end),' ',...
                    'TMP/results_srf_tmpnm_500x_rep',num2str(rep_num),'_',num2str(count_batch),' ',...
                    'TMP/results_sum_tmpnm_500x_rep',num2str(rep_num),'_',num2str(count_batch),...
                    ' >> TMP/results_RTlog_tmpnm_500x_',num2str(rep_num),'_',num2str(count_batch)];
                %system(['start ',tmp_batch_cmd])
                batch_cmd_all={batch_cmd_all{:},tmp_bash_cmd};
            end
        end
        
        switch GLB_RT_MODE
            case 1 %local parfor mode
                delete(gcp('nocreate'));
                parpool(8);
                parfor loop_bash=1:size(batch_cmd_all,2)
                    system(batch_cmd_all{loop_bash});
                end
                time_case(rep_num)=toc;
            case 2 %cluster mode - write all bash commands to sh files
                NUM_sh_files=200; %every 100 bash commands into ONE sh file
                count_num_sh_files=0;
                for loop_bash=1:size(batch_cmd_all,2)
                    if mod(loop_bash-1,NUM_sh_files)==0
                        count_num_sh_files=count_num_sh_files+1;
                        tmp_str_filename=['run_rtjobs_shfile_',num2str(count_num_sh_files),'.sh'];
                        fileID = fopen(tmp_str_filename,'w');

                        fprintf(fileID, '#! /bin/bash\n');
                        fprintf(fileID, '#SBATCH -p lowmem\n');
                        fprintf(fileID, '#SBATCH --mem=200m\n');
                        fprintf(fileID, '#SBATCH -N 1\n');
                        fprintf(fileID, '#SBATCH -n 1\n');
                    end
                    
                    fprintf(fileID,'%s\n',batch_cmd_all{loop_bash});

                    if mod(loop_bash,NUM_sh_files)==0 || loop_bash==size(batch_cmd_all,2)
                        fclose(fileID);
                    end
                end
                
            case 3 %test debug mode
                for loop_bash=1:100%size(batch_cmd_all,2)
                    system(batch_cmd_all{loop_bash});
                end
                time_case(rep_num)=toc;

            case 4 %check results_sum files
                count_batch=0;
                record_success=zeros(1,num_loop_x*num_loop_y);
                for loop_x=1:num_loop_x
                    for loop_y=1:num_loop_y
                        count_batch=count_batch+1;
                        tmp_file_name=['TMP/results_sum_tmpnm_500x_rep',num2str(rep_num),'_',num2str(count_batch)];
                        if(exist(tmp_file_name,'file')==2)
                            record_success(count_batch)=1;
                        end
                    end
                end
                failed_threads{rep_num}=find(record_success==0);
                disp(['Of all ',num2str(num_loop_x*num_loop_y),' parallel ray tracing tasks,']);
                disp([num2str(numel(find(record_success==0))),' tasks failed. The fail rate is ',num2str(numel(find(record_success==0))/(num_loop_x*num_loop_y)),'.'])
        end

        cd ../1.e_geom/
        toc

    case 5 %% 
        cd ../2.e_raytracing/
        %% trace_recal --> ab profiles under blue and red
        %rep_num=1;
        num_layer4lightprofile=10;
        car_chl_ratio=1/6.3;%% carotenoid/chl ratio=1/6.3
        PAR_SAC.blue445nm=[8.551e-05*100, (7.05277e+04+car_chl_ratio*1.68474e+05)*1e-4];
        PAR_SAC.green560nm=[0.000672*100, (1.1048e+04)*1e-4];
        PAR_SAC.red640nm=[0.003111*100, (3.12475e+04)*1e-4];
        %% recal blue
        RT_x_range=500;
        RT_y_range=500;
        RT_x_perthread=5;
        RT_y_perthread=5;
        num_loop_x=RT_x_range/RT_x_perthread;
        num_loop_y=RT_y_range/RT_y_perthread;
        tmp_bash_cmd=['./trace_recal ',num2str(PAR_SAC.blue445nm(1),'%e'),' ',num2str(PAR_SAC.blue445nm(2),'%e'),' count_chl4RT ',...
            num2str(RT_x_range),' ',num2str(RT_y_range),' ',...
            num2str(RT_x_perthread),' ',num2str(RT_y_perthread),' ',...
            'TMP/results_abevents_tmpnm_500x_rep',num2str(rep_num),'_ ',...
            'TMP/results_sum_tmpnm_500x_rep',num2str(rep_num),'_ ',...
            num2str(num_layer4lightprofile),' ',...
            'results_merged_abtri_475nm_500x_rep',num2str(rep_num),' ',...
            'results_merged_absrf_475nm_500x_rep',num2str(rep_num),' ',...
            'results_merged_abprofile_475nm_500x_rep',num2str(rep_num),'_layerN',num2str(num_layer4lightprofile),' ',...
            'results_merged_rtsum_475nm_500x_rep',num2str(rep_num)];
        system(tmp_bash_cmd);
        %% recal red
        RT_x_range=500;
        RT_y_range=25;
        RT_x_perthread=5;
        RT_y_perthread=1;
        num_loop_x=RT_x_range/RT_x_perthread;
        num_loop_y=RT_y_range/RT_y_perthread;
        tmp_bash_cmd=['./trace_recal ',num2str(PAR_SAC.red640nm(1),'%e'),' ',num2str(PAR_SAC.red640nm(2),'%e'),' count_chl4RT ',...
            num2str(RT_x_range),' ',num2str(RT_y_range),' ',...
            num2str(RT_x_perthread),' ',num2str(RT_y_perthread),' ',...
            'TMP/results_abevents_tmpnm_500x_rep',num2str(rep_num),'_ ',...
            'TMP/results_sum_tmpnm_500x_rep',num2str(rep_num),'_ ',...
            num2str(num_layer4lightprofile),' ',...
            'results_merged_abtri_625nm_500x_rep',num2str(rep_num),' ',...
            'results_merged_absrf_625nm_500x_rep',num2str(rep_num),' ',...
            'results_merged_abprofile_625nm_500x_rep',num2str(rep_num),'_layerN',num2str(num_layer4lightprofile),' ',...
            'results_merged_rtsum_625nm_500x_rep',num2str(rep_num)];
        system(tmp_bash_cmd);

        %% tar intermediate files from ray tracing
        %system('tar -zcf results_files_abevents.tar.gz results_abevents_* --remove-files');% tar -zxvf ***.tar.gz; -v = output log
        %system('tar -zcf results_files_srf.tar.gz results_srf_* --remove-files');
        %system('tar -zcf results_files_sum.tar.gz results_sum_* --remove-files');
        %system('tar -zcf results_files_RTlog.tar.gz results_RTlog_* --remove-files');

    case 6
        % %% [3/?] RT Training
        % %% Have measured light absorptance profile to incorporate?
        % %% Start eLeaf_dicot Ray Tracing training process
        % %% Includes
        % %% e_RT_training.m
        %
    case 7
        % %% [4/?] e_geo_physics: reaction-diffusion system
end

end
