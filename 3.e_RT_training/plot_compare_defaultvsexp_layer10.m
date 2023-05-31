%fun_RT_traing_v0_1()

%% RT training
%% parameters: 1) chlorophyll ratio; 2) SAC under RGB

%% blue-445nm; green-561nm; red-638nm
%%%% run 500*500 rays simulation by splitting into 100*100 tasks.
%%%% use minimum of SAC_water & SAC_chl under RGB, then recalculate
%%%% light absorption under RGB ligth
%% [water, chl] SAC default = [7.5e-05*100, 4.94e+04*1e-4] blue 445nm
%%                          = [0.000619*100, 1.22e+04*1e-4] green 560nm
%%                          = [0.003108*100, 2.53e+04*1e-4] red 640nm
%%                          = [0.000114*100, 4.26e+04*1e-4] blue 475nm
%%                          = [0.002834*100, 2.34e+04*1e-4] red 625nm
%% unit = [m-1, m2 g-1]
%% exp for light ab profile measurement: blue=445nm; green=561nm; red=638nm;
%% exp for Licor red-blue measurement: blue=475nm; red=625nm;
cd ../2.e_raytracing/
%% trace_recal --> ab profiles under 445nm, 560nm and 640nm
num_layer4lightprofile=10;
rep_num=1;
%%%% 445nm
SAC_default_blue445nm=[7.5e-05*100, 4.94e+04*1e-4];
RT_x_range=500;
RT_y_range=500;
RT_x_perthread=5;
RT_y_perthread=5;
num_loop_x=RT_x_range/RT_x_perthread;
num_loop_y=RT_y_range/RT_y_perthread;
tmp_bash_cmd=['./trace_recal ',num2str(SAC_default_blue445nm(1),'%e'),' ',num2str(SAC_default_blue445nm(2),'%e'),' count_chl4RT ',...
    num2str(RT_x_range),' ',num2str(RT_y_range),' ',...
    num2str(RT_x_perthread),' ',num2str(RT_y_perthread),' ',...
    'TMP/results_abevents_tmpnm_500x_rep',num2str(rep_num),'_ ',...
    'TMP/results_sum_tmpnm_500x_rep',num2str(rep_num),'_ ',...
    num2str(num_layer4lightprofile),' ',...
    'results_merged_abtri_445nm_500x_c1defaultSAC',' ',...
    'results_merged_absrf_445nm_500x_c1defaultSAC',' ',...
    'results_merged_abprofile_445nm_500x_c1defaultSAC','_layerN',num2str(num_layer4lightprofile),' ',...
    'results_merged_rtsum_445nm_500x_c1defaultSAC'];
system(tmp_bash_cmd);
tmp_file_name=['results_merged_rtsum_445nm_500x_c1defaultSAC'];
tmp_rtsum=importdata(tmp_file_name);
re_445nm(rep_num)=tmp_rtsum.data(2)/tmp_rtsum.data(5)
tr_445nm(rep_num)=tmp_rtsum.data(3)/tmp_rtsum.data(5)
ab_445nm(rep_num)=tmp_rtsum.data(1)/tmp_rtsum.data(5)
tmp_file_name=['results_merged_abprofile_445nm_500x_c1defaultSAC','_layerN',num2str(num_layer4lightprofile)];
tmp_abprofile=importdata(tmp_file_name);
abprofile_445nm(rep_num,:)=tmp_abprofile/tmp_rtsum.data(5);

%%%% 560nm
SAC_default_green560nm=[0.000619*100, 1.22e+04*1e-4];
RT_x_range=500;
RT_y_range=500;
RT_x_perthread=5;
RT_y_perthread=5;
num_loop_x=RT_x_range/RT_x_perthread;
num_loop_y=RT_y_range/RT_y_perthread;
tmp_bash_cmd=['./trace_recal ',num2str(SAC_default_green560nm(1),'%e'),' ',num2str(SAC_default_green560nm(2),'%e'),' count_chl4RT ',...
    num2str(RT_x_range),' ',num2str(RT_y_range),' ',...
    num2str(RT_x_perthread),' ',num2str(RT_y_perthread),' ',...
    'TMP/results_abevents_tmpnm_500x_rep',num2str(rep_num),'_ ',...
    'TMP/results_sum_tmpnm_500x_rep',num2str(rep_num),'_ ',...
    num2str(num_layer4lightprofile),' ',...
    'results_merged_abtri_560nm_500x_c1defaultSAC',' ',...
    'results_merged_absrf_560nm_500x_c1defaultSAC',' ',...
    'results_merged_abprofile_560nm_500x_c1defaultSAC','_layerN',num2str(num_layer4lightprofile),' ',...
    'results_merged_rtsum_560nm_500x_c1defaultSAC'];
system(tmp_bash_cmd);
tmp_file_name=['results_merged_rtsum_560nm_500x_c1defaultSAC'];
tmp_rtsum=importdata(tmp_file_name);
re_560nm(rep_num)=tmp_rtsum.data(2)/tmp_rtsum.data(5)
tr_560nm(rep_num)=tmp_rtsum.data(3)/tmp_rtsum.data(5)
ab_560nm(rep_num)=tmp_rtsum.data(1)/tmp_rtsum.data(5)
tmp_file_name=['results_merged_abprofile_560nm_500x_c1defaultSAC','_layerN',num2str(num_layer4lightprofile)];
tmp_abprofile=importdata(tmp_file_name);
abprofile_560nm(rep_num,:)=tmp_abprofile/tmp_rtsum.data(5);

%%%% 640nm
SAC_default_red640nm=[0.003108*100, 2.53e+04*1e-4];
RT_x_range=500;
RT_y_range=500;
RT_x_perthread=5;
RT_y_perthread=5;
num_loop_x=RT_x_range/RT_x_perthread;
num_loop_y=RT_y_range/RT_y_perthread;
tmp_bash_cmd=['./trace_recal ',num2str(SAC_default_red640nm(1),'%e'),' ',num2str(SAC_default_red640nm(2),'%e'),' count_chl4RT ',...
    num2str(RT_x_range),' ',num2str(RT_y_range),' ',...
    num2str(RT_x_perthread),' ',num2str(RT_y_perthread),' ',...
    'TMP/results_abevents_tmpnm_500x_rep',num2str(rep_num),'_ ',...
    'TMP/results_sum_tmpnm_500x_rep',num2str(rep_num),'_ ',...
    num2str(num_layer4lightprofile),' ',...
    'results_merged_abtri_640nm_500x_c1defaultSAC',' ',...
    'results_merged_absrf_640nm_500x_c1defaultSAC',' ',...
    'results_merged_abprofile_640nm_500x_c1defaultSAC','_layerN',num2str(num_layer4lightprofile),' ',...
    'results_merged_rtsum_640nm_500x_c1defaultSAC'];
system(tmp_bash_cmd);
tmp_file_name=['results_merged_rtsum_640nm_500x_c1defaultSAC'];
tmp_rtsum=importdata(tmp_file_name);
re_640nm(rep_num)=tmp_rtsum.data(2)/tmp_rtsum.data(5)
tr_640nm(rep_num)=tmp_rtsum.data(3)/tmp_rtsum.data(5)
ab_640nm(rep_num)=tmp_rtsum.data(1)/tmp_rtsum.data(5)
tmp_file_name=['results_merged_abprofile_640nm_500x_c1defaultSAC','_layerN',num2str(num_layer4lightprofile)];
tmp_abprofile=importdata(tmp_file_name);
abprofile_640nm(rep_num,:)=tmp_abprofile/tmp_rtsum.data(5);

cd ../3.e_RT_training/
load("ab_profile_exp_layer10.mat","re_*","tr_*","ab_layer_*");
figure
subplot(1,3,1)
plot(ab_layer_445nm_WT,'b-o');hold on;
plot(ab_layer_445nm_LCD,'r-*');hold on;
plot(abprofile_445nm(10:-1:1),'b--')
ylim([0,0.3])
subplot(1,3,2)
plot(ab_layer_561nm_WT,'b-o');hold on;
plot(ab_layer_561nm_LCD,'r-*');hold on;
plot(abprofile_560nm(10:-1:1),'b--')
ylim([0,0.3])
subplot(1,3,3)
plot(ab_layer_638nm_WT,'b-o');hold on;
plot(ab_layer_638nm_LCD,'r-*');hold on;
plot(abprofile_640nm(10:-1:1),'b--')
ylim([0,0.3])
set(gcf,'Color','w');

re_445nm_WT
tr_445nm_WT
1-re_445nm_WT-tr_445nm_WT

re_561nm_WT
tr_561nm_WT
1-re_561nm_WT-tr_561nm_WT

re_638nm_WT
tr_638nm_WT
1-re_638nm_WT-tr_638nm_WT

save compare_defaultvsexp_layer10.mat
%end