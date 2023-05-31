% eLeaf: 3D model of dicot leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 0.1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update 2020-Feb
% visualization meshed ply files

% celltype
%     0=leaf
%     1=ms, i.e. pal or spongy
%     2=ms_chl
%     3=non ms
%     4=ms_vac

clear;clc
%colormap([0.3,0.3,0.3;0.2,0.8,0.2;0.2,0.2,0.8]);
alpha=0.4%0.25;
linecolor=[0.5,0.5,0.5];%'none';%[0.8,0.8,0.8];
color_wall=[0.3,0.3,0.3];
color_chl=[0.2,0.8,0.2];
color_vac=[0.2,0.2,0.8];

load('SAVE_RT_geo_export.mat', 'GLB_tag_*')

%%%% 2. upper epidermis
for loop_i=1:numel(GLB_tag_EPL_u_set)
    tmp_name=['../nonMS/EPL_u_',num2str(loop_i),'.ply'];
    [tri,pts]=ply_read(tmp_name,'tri');
    trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'Facecolor',color_wall,'facealpha',alpha,'EdgeColor',linecolor);hold on;
end

%%%% 3-5. PAL
for loop_i=4:4%numel(GLB_tag_PAL_MS)
    %%%% 3. PAL MS
    tmp_name=['../MS/PAL_MS_',num2str(loop_i),'.ply'];
    [tri,pts]=ply_read(tmp_name,'tri');
    trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'Facecolor',color_wall,'facealpha',alpha,'EdgeColor',linecolor);hold on;
    
    %%%% 4. PAL VAC
    tmp_name=['../MS/PAL_VAC_',num2str(loop_i),'.ply'];
    [tri,pts]=ply_read(tmp_name,'tri');
    trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'Facecolor',color_vac,'facealpha',alpha,'EdgeColor',linecolor);hold on;
    
    %%%% 5. PAL CHL
    tmptag_set_PAL_CHL=GLB_tag_PAL_CHL{loop_i};
    for loop_i_CHL=1:numel(tmptag_set_PAL_CHL)
        tmp_name=['../MS/PAL_MS_',num2str(loop_i),'_CHL_',num2str(loop_i_CHL),'.ply'];
        [tri,pts]=ply_read(tmp_name,'tri');
        trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'Facecolor',color_chl,'facealpha',alpha,'EdgeColor',linecolor);hold on;
    end 
end

%%%% 6-8. SPO
for loop_i=11:11%numel(GLB_tag_SPO_MS)
    %%%% 6. SPO MS
    tmp_name=['../MS/SPO_MS_',num2str(loop_i),'.ply'];
    [tri,pts]=ply_read(tmp_name,'tri');
    %trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'Facecolor',color_wall,'facealpha',alpha,'EdgeColor',linecolor);hold on;
    
    %%%% 7. SPO VAC
    tmp_name=['../MS/SPO_VAC_',num2str(loop_i),'.ply'];
    [tri,pts]=ply_read(tmp_name,'tri');
    trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'Facecolor',color_vac,'facealpha',alpha,'EdgeColor',linecolor);hold on;
    
    %%%% 8. SPO CHL
    tmptag_set_SPO_CHL=GLB_tag_SPO_CHL{loop_i};
    for loop_i_CHL=1:numel(tmptag_set_SPO_CHL)
        tmp_name=['../MS/SPO_MS_',num2str(loop_i),'_CHL_',num2str(loop_i_CHL),'.ply'];
        [tri,pts]=ply_read(tmp_name,'tri');
        trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'Facecolor',color_chl,'facealpha',alpha,'EdgeColor',linecolor);hold on;
    end
end

%%%% 1. lower epidermis
for loop_i=1:numel(GLB_tag_EPL_l_set)
    tmp_name=['../nonMS/EPL_l_',num2str(loop_i),'.ply'];
    [tri,pts]=ply_read(tmp_name,'tri');
    trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'Facecolor',color_wall,'facealpha',alpha,'EdgeColor',linecolor);hold on;
end

axis equal;
%% plot ray
P=[1.152749e-05 1.213462e-05 2.459535e-04];
D=[2.387664e-01 4.533029e-01 -8.587823e-01];
dis=3e-06;%default 6e-6
plot3([P(1);P(1)+D(1)*dis],[P(2);P(2)+D(2)*dis],[P(3);P(3)+D(3)*dis],'r-o','linewidth',1);hold on;

%%%% just another line
P2=[1.152453e-05 1.212916e-05 2.459641e-04];
D2=[2.418012e-01 4.474303e-01 -8.610100e-01];
dis2=1.221707e-08;
plot3([P2(1);P2(1)+D2(1)*dis2],[P2(2);P2(2)+D2(2)*dis2],[P2(3);P2(3)+D2(3)*dis2],'b-o','linewidth',1);hold on;

zlim([P(3)-5e-6,P(3)+5e-6])

%%%% series of points
P3=[1.078820e-05 1.080735e-05 2.484794e-04;...,
    ];
plot3([P3(:,1);P2(1)],[P3(:,2);P2(2)],[P3(:,3);P2(3)],'b-o','linewidth',1);hold on;
axis equal
%%%% control region
zlim([P(3)-5e-6,P(3)+5e-6])

%% plot triangle
p1=[0, 1.6099999999999998e-05, 0.00023186388072289201];
p2=[1.1458979169813601e-06, 1.6099999999999998e-05, 0.00023259194047620799];
p3=[1.1600600540578899e-06, 1.6099999999999998e-05, 0.00023110205996436801];
plot3([p1(1);p2(1);p3(1);p1(1)],[p1(2);p2(2);p3(2);p1(2)],[p1(3);p2(3);p3(3);p1(3)],'b-','linewidth',1);hold on;
axis equal

tmpp1=[4.043143e-06 1.117360e-06 1.582415e-04];
tmpp2=[0.000000e+00 4.314533e-07 1.577182e-04];
tmpp3=[2.543243e-06 -5.293956e-23 1.573891e-04];
%tmpp4=[2.543243e-06 -3.761582e-37 1.573891e-04];
tmppall=[tmpp1;tmpp2;tmpp3];
plot3(tmppall(:,1)',tmppall(:,2)',tmppall(:,3)','r-o','linewidth',1);
hold on; axis equal;zlim([1.54e-4,1.6e-4])

%% plot RT_file4plot
tmp_paths_all=load('../tmp_debug_4633');
tmp_path4plot=tmp_paths_all(1:1035,3:5);
plot3(tmp_path4plot(:,1),tmp_path4plot(:,2),tmp_path4plot(:,3),'r-o','linewidth',1);hold on;
axis equal

%% for visualize the inner_check results
%% e.g. "[Warning] SPO_MS_1_CHL_1.ply PTS 37 is not inside SPO_MS_1.ply"
loop_i=10;
loop_i_CHL=25;
k_pts=8;
%%%% 6. SPO MS
tmp_name=['../MS/SPO_MS_',num2str(loop_i),'.ply'];
%[tri,pts]=ply_read_xy(tmp_name);
[tri,pts]=ply_read(tmp_name,'tri');
trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'Facecolor',color_wall,'facealpha',alpha,'EdgeColor',linecolor);hold on;

%%%% 7. SPO VAC
tmp_name=['../MS/SPO_VAC_',num2str(loop_i),'.ply'];
%[tri,pts]=ply_read_xy(tmp_name);
[tri,pts]=ply_read(tmp_name,'tri');
trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'Facecolor',color_vac,'facealpha',alpha,'EdgeColor',linecolor);hold on;

%%%% 8. SPO CHL
tmptag_set_SPO_CHL=GLB_tag_SPO_CHL{loop_i};
tmp_name=['../MS/SPO_MS_',num2str(loop_i),'_CHL_',num2str(loop_i_CHL),'.ply'];
%[tri,pts]=ply_read_xy(tmp_name);
[tri,pts]=ply_read(tmp_name,'tri');
trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'Facecolor',color_chl,'facealpha',alpha,'EdgeColor',linecolor);hold on;

P=pts(k_pts+1,:); %% since in C programme inner_check index of pts start from 0
D=[0,0,-1];
dis=6e-6;
plot3([P(1);P(1)+D(1)*dis],[P(2);P(2)+D(2)*dis],[P(3);P(3)+D(3)*dis],'r-o','linewidth',1);hold on;
axis equal

%%demo how to locate triangle with xyz from MATLAB figure
% p1=[6.923e-7,2.721e-6,7.389e-5];
% p2=[6.923e-7,3.327e-6,7.392e-5];
% p3=[1.119e-6,3.147e-6,7.392e-5];
% tmp=pts-p3;
% [value,idx]=min(sum(tmp.*tmp,2))