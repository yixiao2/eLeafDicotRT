% eLeaf: 3D model of dicot leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 0.1

function e_geo_parainput_v0_1_b4fit(CFG_PARA_COM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0 for wild-type; 1 for LCD1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%select_com=input('selection (e.g [1 1 0 0 0 0 0 0 0]):');
select_com=CFG_PARA_COM;
if size(select_com,2)~=13
    disp('input selection error');
    return;
end

%% other for model
MODEL_ddis=0.1e-6;%small distance to avoid point intersect
FLAG_debug_msg=0;

%% 1.upper epidermis
if select_com(1)==0
    EPL_u_thick  = 22e-6;
    EPL_u_e      = 0.77;%eccentricity: sqrt(1-b^2/a^2)
    EPL_u_flat   = sqrt(1-EPL_u_e^2);%flatness of elipsoid: b/a
    EPL_u_radius = EPL_u_thick/2/(1-EPL_u_flat)*0.75;%cut at 3/4 radius;hexagon cylinder
    EPL_u_rand_displacex=1/2*EPL_u_radius*rand;
    EPL_u_rand_displacey=sqrt(3)/2*EPL_u_radius*rand;
else
    EPL_u_thick  = 22e-6;
    EPL_u_e      = 0.77;%eccentricity: sqrt(1-b^2/a^2)
    EPL_u_flat   = sqrt(1-EPL_u_e^2);%flatness of elipsoid: b/a
    EPL_u_radius = EPL_u_thick/2/(1-EPL_u_flat)*0.75;%cut at 3/4 radius
    EPL_u_rand_displacex=1/2*EPL_u_radius*rand;
    EPL_u_rand_displacey=sqrt(3)/2*EPL_u_radius*rand;
end

%% 2. lower epidermis
if select_com(2)==0
    EPL_l_thick  = 16e-6;
    EPL_l_e      = 0.77;%eccentricity: sqrt(1-b^2/a^2)
    EPL_l_flat   = sqrt(1-EPL_l_e^2);%flatness of elipsoid: b/a
    EPL_l_radius = EPL_l_thick/2/(1-EPL_l_flat)*0.75;%cut at 3/4 radius
    EPL_l_rand_displacex=1/2*EPL_l_radius*rand;
    EPL_l_rand_displacey=sqrt(3)/2*EPL_l_radius*rand;
else
    EPL_l_thick  = 16e-6;
    EPL_l_e      = 0.77;%eccentricity: sqrt(1-b^2/a^2)
    EPL_l_flat   = sqrt(1-EPL_l_e^2);%flatness of elipsoid: b/a
    EPL_l_radius = EPL_l_thick/2/(1-EPL_l_flat)*0.75;%cut at 3/4 radius
    EPL_l_rand_displacex=1/2*EPL_l_radius*rand;
    EPL_l_rand_displacey=sqrt(3)/2*EPL_l_radius*rand;
end

%% 3. palisade tissue thickness
if select_com(3)==0
    PAL_box_thick= 112.1e-6;
else
    PAL_box_thick= 124.3e-6;
end
%% 4. palisade tissue porosity
if select_com(4)==0
    MS_porosity = 0.328;
else
    MS_porosity = 0.422;
end
%% 5. palisade cell radius
if select_com(5)==0
    PAL_MS_radius=22.3e-6/2;
else
    PAL_MS_radius=22.3e-6/2;
end
%% 6. palisade chloroplast coverage - ScperSm
if select_com(6)==0
    PAL_chl_ScperSm = 0.85; %Sc/Sm for chloroplasts in PAL
else
    PAL_chl_ScperSm = 0.85; %Sc/Sm for chloroplasts in PAL
end
%% 7. palisade chloroplast size
if select_com(7)==0
    tmp_PAL_radius_chl_sphere = 5.45e-6/2;% [to revise] can increase a bit?
    PAL_radius_chl_sphere = tmp_PAL_radius_chl_sphere; %radius of sphere
    PAL_diameter_chl_patch = tmp_PAL_radius_chl_sphere*2*0.8; %distance to search centroid to distribute chloroplasts
    PAL_thick_chl_layer = 2.87e-6; %thickness of chloroplast layer; later split by spheres to chloroplasts
else
    tmp_PAL_radius_chl_sphere = 5.45e-6/2;% [to revise] can increase a bit?
    PAL_radius_chl_sphere = tmp_PAL_radius_chl_sphere; %radius of sphere
    PAL_diameter_chl_patch = tmp_PAL_radius_chl_sphere*2*0.8; %distance to search centroid to distribute chloroplasts
    PAL_thick_chl_layer = 2.87e-6; %thickness of chloroplast layer; later split by spheres to chloroplasts
end
%% other for pal
N_PAL_cols=2;
N_PAL_rows=2;
if N_PAL_cols==N_PAL_rows
    PAL_box_length=PAL_MS_radius*(N_PAL_cols*2-2)+MODEL_ddis*(N_PAL_cols-1);
else
    error('INPUT not compatible. N_PAL_cols does not equal to N_PAL_rows.')
    %% to develop in future
end
PAL_MS_ucap_flat=0.7;
PAL_MS_lcap_flat=0.7;
N_EPL_u_cols=ceil((PAL_box_length+EPL_u_rand_displacex)/EPL_u_radius/1.5-0.5/1.5)+1;%initial displace 1/1.5
N_EPL_u_rows=ceil((PAL_box_length+EPL_u_rand_displacey)/EPL_u_radius/sqrt(3)+0.5);%initial displace 0.5 cell
N_EPL_l_cols=ceil((PAL_box_length+EPL_l_rand_displacex)/EPL_l_radius/1.5-0.5/1.5)+1;%initial displace 1/1.5
N_EPL_l_rows=ceil((PAL_box_length+EPL_l_rand_displacey)/EPL_l_radius/sqrt(3)+0.5);%initial displace 0.5 cell
PAL_dis_chl2wall = 0.3e-6; %thickness of cytosol between cell wall and chloroplast
PAL_dis_mit2chl_rt = 0.08; %distance from mitochondrion layer to inner boundary of CHL, suffix rt means ratio (relative to PAL_MS_radius)
PAL_radius_mit_rt = 0.04; %radius of mitochondria, rt (ratio) to PAL_MS_radius
PAL_dis_vac2mit_rt = 0.08; %PAL vac=CHL inner boundary - PAL_dis_mit2chl_rt - PAL_radius_mit_rt - PAL_dis_vac2mit_rt
PAL_diameter_mit_patch = 5e-6; %distance to search centroid to distribute mitochondria

%% 8. spongy tissue thickness
if select_com(8)==0
    SPO_box_thick= 140.6e-6;
else
    SPO_box_thick= 185.0e-6;
end
%% 9. spongy tissue porosity
if select_com(9)==0
    %MS_porosity = 0.328;
else
    %MS_porosity = 0.422;
end
%% 10. spongy cell radius
if select_com(10)==0
    SPO_MS_radius=19.3e-6/2;
else
    SPO_MS_radius=21.7e-6/2;
end
%% 11. spongy chloroplast coverage - ScperSm
if select_com(11)==0
    SPO_chl_ScperSm = 0.85; %Sc/Sm for chloroplasts in SPO
else
    SPO_chl_ScperSm = 0.85; %Sc/Sm for chloroplasts in SPO
end
%% 12. spongy chloroplast size
if select_com(12)==0
    tmp_SPO_radius_chl_sphere = 5.45e-6/2;
    SPO_radius_chl_sphere = tmp_SPO_radius_chl_sphere; %radius of sphere
    SPO_diameter_chl_patch = tmp_SPO_radius_chl_sphere*2*0.7; %distance to search centroid to distribute chloroplasts
    SPO_thick_chl_layer = 2.87e-6; %thickness of chloroplast layer; later split by spheres to chloroplasts
else
    tmp_SPO_radius_chl_sphere = 5.45e-6/2;
    SPO_radius_chl_sphere = tmp_SPO_radius_chl_sphere; %radius of sphere
    SPO_diameter_chl_patch = tmp_SPO_radius_chl_sphere*2*0.7; %distance to search centroid to distribute chloroplasts
    SPO_thick_chl_layer = 2.87e-6; %thickness of chloroplast layer; later split by spheres to chloroplasts
end
%% other for spo
SPO_box_length=PAL_box_length;
dis_SPO_MS_center = SPO_MS_radius*2*0.8;% minimum distance between two centers
SPO_dis_chl2wall = 0.3e-6; %thickness of cytosol between cell wall and chloroplast
SPO_dis_mit2chl_rt = 0.08; %distance from mitochondrion layer to inner boundary of CHL, suffix rt means ratio (relative to SPO_MS_radius)
SPO_radius_mit_rt = 0.04; %radius of mitochondria, rt (ratio) to PAL_MS_radius
SPO_dis_vac2mit_rt = 0.08; %SPO vac=CHL inner boundary - SPO_dis_mit2chl_rt - SPO_radius_mit_rt - SPO_dis_vac2mit_rt
SPO_diameter_mit_patch = 3e-6; %distance to search centroid to distribute mitochondria

%% 13 chlorophyll
if select_com(13)==0
    %convert [chl] from mg m-2 to g/m3
    % ([chl]*1e-3)/(AVR);% convertion in e_geo_AVRcal_v0_1.m when generating chl_vol_profile
    chl_con=463.4425;%mg/m2
    chl_con_ratio_PvS=1.0;
else
    %convert [chl] from mg m-2 to g/m3
    % ([chl]*1e-3)/(AVR);% convertion in e_geo_AVRcal_v0_1.m when generating chl_vol_profile
    chl_con=334.4463;%mg/m2
    chl_con_ratio_PvS=1.0;
end




save parainput.mat
disp('Input parameter loaded.')
end
