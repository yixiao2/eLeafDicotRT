% eLeaf: 3D model of dicot leaf photosynthesis
% @license: LGPL (GNU LESSER GENERAL PUBLIC LICENSE Version 3)
% @author: Yi Xiao <yixiao20@outlook.com>
% @version: 0.1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update 2020-April
% afte manually check all obj and save correct ply from Debug_*.mph
% visualize whole leaf again

clearvars
%colormap([0.3,0.3,0.3;0.2,0.8,0.2;0.2,0.2,0.8]);
tmpcolormap=[0.5,0.5,0.5;0.2,0.8,0.2;0.2,0.2,0.8];
alpha=0.3;%0.25;
linecolor=[0.3,0.3,0.3];%'none';%[0.8,0.8,0.8];
wall=0;
chl=1;
vac=2;

%input GLB_filedir_mode
prompt = 'Current folder for ply files:\n1=../MS & ../nonMS, 2=./geo_export\n';
GLB_filedir_mode = input(prompt);
%input GLB_MSCHLcheck_mode
prompt = 'Mode to check PAL&SPO CHL:\n1=pause at each CHL, 2=pause at each PAL&SPO\n';
GLB_MSCHLcheck_mode = input(prompt);
if GLB_MSCHLcheck_mode==1
    prompt = 'chl_check_mode=1, then which PAL&SPO to check?\n 1=load from MAT, 2=user input, 3=all\n';
    GLB_MSCHLcheck_select_mode = input(prompt);
    if GLB_MSCHLcheck_select_mode==1
        load SAVE_RT_aftervisualcheck_perMS.mat SAVE_FAIL_LIST_ext_PAL_VACCHL SAVE_FAIL_LIST_ext_SPO_VACCHL
        TMP_FAIL_LIST_PAL_VACCHL=SAVE_FAIL_LIST_ext_PAL_VACCHL;
        TMP_FAIL_LIST_SPO_VACCHL=SAVE_FAIL_LIST_ext_SPO_VACCHL;
    elseif GLB_MSCHLcheck_select_mode==2
        prompt = 'Input No. of PAL you would like to inspect\n';
        TMP_FAIL_LIST_PAL_VACCHL=input(prompt);
        prompt = 'Input No. of SPO you would like to inspect\n';
        TMP_FAIL_LIST_SPO_VACCHL=input(prompt);
        %%input format e.g.: [2,4] or 1:4 or 
    else
        %%do nothing
        TMP_FAIL_LIST_PAL_VACCHL=1:numel(GLB_tag_PAL_MS);
        TMP_FAIL_LIST_SPO_VACCHL=1:numel(GLB_tag_SPO_MS);
    end
end
%input GLB_manualcheck_save_mode
prompt = 'Mode to save result:\n1=SAVE_RT_aftervisualcheck_full.mat, 2=TMP_SAVE_RT_aftervisualcheck\n';
GLB_manualcheck_save_mode = input(prompt);

load('SAVE_RT_geo_export.mat', 'GLB_tag_*', 'SAVE_FAIL_LIST*')
ax1=subplot(1,2,1);
ax2=subplot(1,2,2);

if GLB_manualcheck_save_mode==1
    SAVE_FAIL_LIST_ext_EPL_u=SAVE_FAIL_LIST_EPL_u;
    SAVE_FAIL_LIST_ext_EPL_l=SAVE_FAIL_LIST_EPL_l;
    SAVE_FAIL_LIST_ext_PAL_MS=SAVE_FAIL_LIST_PAL_MS;
    SAVE_FAIL_LIST_ext_PAL_VAC=SAVE_FAIL_LIST_PAL_VAC;
    SAVE_FAIL_LIST_ext_PAL_CHL=SAVE_FAIL_LIST_PAL_CHL;
    SAVE_FAIL_LIST_ext_PAL_VACCHL=[];
    SAVE_FAIL_LIST_ext_SPO_MS=SAVE_FAIL_LIST_SPO_MS;
    SAVE_FAIL_LIST_ext_SPO_VACCHL=[];
    SAVE_FAIL_LIST_ext_SPO_VAC=SAVE_FAIL_LIST_SPO_VAC;%only used in 2nd batch
    SAVE_FAIL_LIST_ext_SPO_CHL=SAVE_FAIL_LIST_SPO_CHL;%only used in 2nd batch
else
    SAVE_FAIL_LIST_ext_EPL_u=[];
    SAVE_FAIL_LIST_ext_EPL_l=[];
    SAVE_FAIL_LIST_ext_PAL_MS=[];
    SAVE_FAIL_LIST_ext_PAL_VAC=[];
    SAVE_FAIL_LIST_ext_PAL_CHL=[];
    SAVE_FAIL_LIST_ext_PAL_VACCHL=[];
    SAVE_FAIL_LIST_ext_SPO_MS=[];
    SAVE_FAIL_LIST_ext_SPO_VACCHL=[];
    SAVE_FAIL_LIST_ext_SPO_VAC=[];
    SAVE_FAIL_LIST_ext_SPO_CHL=[];
end

%%%% 2. upper epidermis
%SAVE_FAIL_LIST_ext_EPL_u=SAVE_FAIL_LIST_EPL_u;
TMP_count_fail=size(SAVE_FAIL_LIST_EPL_u,1);
for loop_i=1:numel(GLB_tag_EPL_u_set)
    if GLB_filedir_mode==1
        tmp_name=['../nonMS/EPL_u_',num2str(loop_i),'.ply'];
    else
        tmp_name=['./EPL_u_',num2str(loop_i),'.ply'];
    end
    try
        %[tri,pts]=ply_read_xy(tmp_name);
        [tri,pts]=ply_read(tmp_name,'tri');
        axes(ax1);%subplot(1,2,1)%cell level
        trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facecolor',[0.5,0.5,0.5],'facealpha',alpha,'EdgeColor',linecolor);axis equal;
        axes(ax2);%subplot(1,2,2)%obj level
        trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facecolor',[0.5,0.5,0.5],'facealpha',alpha,'EdgeColor',linecolor);axis equal;
        %interactive manually check by visualization
        prompt = ['Mesh of ',tmp_name,' seems not OK? 1=true, 0/enter=false: '];
        tmp_str = input(prompt);
        if isempty(tmp_str)
            tmp_str = 0;
        end
        if tmp_str~=0
            TMP_count_fail=TMP_count_fail+1;
            SAVE_FAIL_LIST_ext_EPL_u(TMP_count_fail)=loop_i;
        end
    catch
        disp(['Doesn''t find this object: ''',tmp_name,'''']);
    end
end
%figure(f1);savefig('Fig_EPL_u.fig');close(gcf);

%%%% 3-5. PAL
%SAVE_FAIL_LIST_ext_PAL_MS=SAVE_FAIL_LIST_PAL_MS;
TMP_count_fail_1=size(SAVE_FAIL_LIST_PAL_MS,1);
if GLB_MSCHLcheck_mode==1
    %SAVE_FAIL_LIST_ext_PAL_VAC=SAVE_FAIL_LIST_PAL_VAC;
    TMP_count_fail_2=size(SAVE_FAIL_LIST_PAL_VAC,1);
    %SAVE_FAIL_LIST_ext_PAL_CHL=SAVE_FAIL_LIST_PAL_CHL;
    TMP_count_fail_3=size(SAVE_FAIL_LIST_PAL_CHL,1);
else
    %SAVE_FAIL_LIST_ext_PAL_VACCHL=[];
    TMP_count_fail_4=0;
end

if GLB_MSCHLcheck_mode==1
    TMP_loop_set=TMP_FAIL_LIST_PAL_VACCHL;
else
    TMP_loop_set=1:numel(GLB_tag_PAL_MS);
end
for loop_i=TMP_loop_set
    %%%% 3. PAL MS
    if GLB_filedir_mode==1
        tmp_name=['../MS/PAL_MS_',num2str(loop_i),'.ply'];
    else
        tmp_name=['./PAL_MS_',num2str(loop_i),'.ply'];
    end
    try
        %[tri,pts]=ply_read_xy(tmp_name);
        [tri,pts]=ply_read(tmp_name,'tri');
        axes(ax1);%subplot(1,2,1)
        cla(ax1)
        trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facecolor',[0.5,0.5,0.5],'facealpha',alpha,'EdgeColor',linecolor);hold on;axis equal;
        axes(ax2);%subplot(1,2,2)
        trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facecolor',[0.5,0.5,0.5],'facealpha',alpha,'EdgeColor',linecolor);axis equal;
        %interactive manually check by visualization
        prompt = ['Mesh of ',tmp_name,' seems not OK? 1=true, 0/enter=false: '];
        tmp_str = input(prompt);
        if isempty(tmp_str)
            tmp_str = 0;
        end
        if tmp_str~=0
            TMP_count_fail_1=TMP_count_fail_1+1;
            SAVE_FAIL_LIST_ext_PAL_MS(TMP_count_fail_1)=loop_i;
        end
    catch
        disp(['Doesn''t find this object: ''',tmp_name,'''']);
    end
    %close(gcf)
    %%%% 4. PAL VAC
    if GLB_filedir_mode==1
        tmp_name=['../MS/PAL_VAC_',num2str(loop_i),'.ply'];
    else
        tmp_name=['./PAL_VAC_',num2str(loop_i),'.ply'];
    end
    try
        %[tri,pts]=ply_read_xy(tmp_name);
        [tri,pts]=ply_read(tmp_name,'tri');
        axes(ax1);%subplot(1,2,1)
        cla(ax1)
        trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facecolor',[0.2,0.2,0.8],'facealpha',alpha,'EdgeColor',linecolor);hold on;axis equal;
        axes(ax2);%subplot(1,2,2)
        trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facecolor',[0.2,0.2,0.8],'facealpha',alpha,'EdgeColor',linecolor);axis equal;
        if GLB_MSCHLcheck_mode==1
            %interactive manually check by visualization
            prompt = ['Mesh of ',tmp_name,' seems not OK? 1=true, 0/enter=false: '];
            tmp_str = input(prompt);
            if isempty(tmp_str)
                tmp_str = 0;
            end
            if tmp_str~=0
                TMP_count_fail_2=TMP_count_fail_2+1;
                SAVE_FAIL_LIST_ext_PAL_VAC(TMP_count_fail_2)=loop_i;
            end
        end
    catch
        disp(['Doesn''t find this object: ''',tmp_name,'''']);
    end
    
    %%%% 5. PAL CHL
    tmptag_set_PAL_CHL=GLB_tag_PAL_CHL{loop_i};
    for loop_i_CHL=1:numel(tmptag_set_PAL_CHL)
        if GLB_filedir_mode==1
            tmp_name=['../MS/PAL_MS_',num2str(loop_i),'_CHL_',num2str(loop_i_CHL),'.ply'];
        else
            tmp_name=['./PAL_MS_',num2str(loop_i),'_CHL_',num2str(loop_i_CHL),'.ply'];
        end
        try
            %[tri,pts]=ply_read_xy(tmp_name);
            [tri,pts]=ply_read(tmp_name,'tri');
            axes(ax1);%subplot(1,2,1)
            trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facecolor',[0.2,0.8,0.2],'facealpha',alpha,'EdgeColor',linecolor);hold on;axis equal;
            axes(ax2);%subplot(1,2,2)
            trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facecolor',[0.2,0.8,0.2],'facealpha',alpha,'EdgeColor',linecolor);axis equal;
            if GLB_MSCHLcheck_mode==1
                %interactive manually check by visualization
                prompt = ['Mesh of ',tmp_name,' seems not OK? 1=true, 0/enter=false: '];
                tmp_str = input(prompt);
                if isempty(tmp_str)
                    tmp_str = 0;
                end
                if tmp_str~=0
                    TMP_count_fail_3=TMP_count_fail_3+1;
                    SAVE_FAIL_LIST_ext_PAL_CHL(TMP_count_fail_3,1)=loop_i;
                    SAVE_FAIL_LIST_ext_PAL_CHL(TMP_count_fail_3,2)=loop_i_CHL;
                end
            end
        catch
            disp(['Doesn''t find this object: ''',tmp_name,'''']);
        end
    end
    if GLB_MSCHLcheck_mode==2
        tmp_name=['./PAL_VACorCHL',num2str(loop_i),'.ply'];
        prompt = ['Mesh of ',tmp_name,' seems not OK? 1=true, 0/enter=false: '];
        tmp_str = input(prompt);
        if isempty(tmp_str)
            tmp_str = 0;
        end
        if tmp_str~=0
            TMP_count_fail_4=TMP_count_fail_4+1;
            SAVE_FAIL_LIST_ext_PAL_VACCHL(TMP_count_fail_4)=loop_i;
        end
    end
end

%%%% SPO first batch - only stop when finishing one cell
%%%% 6-8. SPO
%SAVE_FAIL_LIST_ext_SPO_MS=SAVE_FAIL_LIST_SPO_MS;
TMP_count_fail_1=size(SAVE_FAIL_LIST_SPO_MS,1);
%SAVE_FAIL_LIST_ext_SPO_VACCHL=[];
TMP_count_fail_4=0;
if GLB_MSCHLcheck_mode==1
    TMP_loop_set=TMP_FAIL_LIST_SPO_VACCHL;
else
    TMP_loop_set=1:numel(GLB_tag_SPO_MS);
end
for loop_i=TMP_loop_set%1:numel(GLB_tag_SPO_MS)
    %%%% 6. SPO MS
    if GLB_filedir_mode==1
        tmp_name=['../MS/SPO_MS_',num2str(loop_i),'.ply'];
    else
        tmp_name=['./SPO_MS_',num2str(loop_i),'.ply'];
    end
    try
        %[tri,pts]=ply_read_xy(tmp_name);
        [tri,pts]=ply_read(tmp_name,'tri');
        axes(ax1);%subplot(1,2,1)
        cla(ax1)
        trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facecolor',[0.5,0.5,0.5],'facealpha',alpha,'EdgeColor',linecolor);hold on;axis equal;
        axes(ax2);%subplot(1,2,1)
        trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facecolor',[0.5,0.5,0.5],'facealpha',alpha,'EdgeColor',linecolor);axis equal;
        %interactive manually check by visualization
        prompt = ['Mesh of ',tmp_name,' seems not OK? 1=true, 0/enter=false: '];
        tmp_str = input(prompt);
        if isempty(tmp_str)
            tmp_str = 0;
        end
        if tmp_str~=0
            TMP_count_fail_1=TMP_count_fail_1+1;
            SAVE_FAIL_LIST_ext_SPO_MS(TMP_count_fail_1)=loop_i;
        end
    catch
        disp(['Doesn''t find this object: ''',tmp_name,'''']);
    end
    %close(gcf)
    %%%% 7. SPO VAC
    if GLB_filedir_mode==1
        tmp_name=['../MS/SPO_VAC_',num2str(loop_i),'.ply'];
    else
        tmp_name=['./SPO_VAC_',num2str(loop_i),'.ply'];
    end
    try
        %[tri,pts]=ply_read_xy(tmp_name);
        [tri,pts]=ply_read(tmp_name,'tri');
        axes(ax1);%subplot(1,2,1)
        cla(ax1)
        trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facecolor',[0.2,0.2,0.8],'facealpha',alpha,'EdgeColor',linecolor);hold on;axis equal;
        axes(ax2);%subplot(1,2,2)
        trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facecolor',[0.2,0.2,0.8],'facealpha',alpha,'EdgeColor',linecolor);axis equal;
%         %interactive manually check by visualization
%         prompt = 'Mesh seems not OK?: 1=true, 0/enter=false ';
%         tmp_str = input(prompt);
%         if isempty(tmp_str)
%             tmp_str = 0;
%         end
%         if tmp_str~=0
%             TMP_count_fail_2=TMP_count_fail_2+1;
%             SAVE_FAIL_LIST_ext_SPO_VAC(TMP_count_fail_2)=loop_i;
%         end
    catch
        disp(['Doesn''t find this object: ''',tmp_name,'''']);
    end
    
    %%%% 8. SPO CHL
    tmptag_set_SPO_CHL=GLB_tag_SPO_CHL{loop_i};
    for loop_i_CHL=1:numel(tmptag_set_SPO_CHL)
        if GLB_filedir_mode==1
            tmp_name=['../MS/SPO_MS_',num2str(loop_i),'_CHL_',num2str(loop_i_CHL),'.ply'];
        else
            tmp_name=['./SPO_MS_',num2str(loop_i),'_CHL_',num2str(loop_i_CHL),'.ply'];
        end
        try
            %[tri,pts]=ply_read_xy(tmp_name);
            [tri,pts]=ply_read(tmp_name,'tri');
            axes(ax1);%subplot(1,2,1)
            trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facecolor',[0.2,0.8,0.2],'facealpha',alpha,'EdgeColor',linecolor);hold on;axis equal;
            axes(ax2);%subplot(1,2,2)
            trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facecolor',[0.2,0.8,0.2],'facealpha',alpha,'EdgeColor',linecolor);axis equal;
%             %interactive manually check by visualization
%             prompt = 'Mesh seems not OK?: 1=true, 0/enter=false ';
%             tmp_str = input(prompt);
%             if isempty(tmp_str)
%                 tmp_str = 0;
%             end
%             if tmp_str~=0
%                 TMP_count_fail_3=TMP_count_fail_3+1;
%                 SAVE_FAIL_LIST_ext_SPO_CHL(TMP_count_fail_3,1)=loop_i;
%                 SAVE_FAIL_LIST_ext_SPO_CHL(TMP_count_fail_3,2)=loop_i_CHL;
%             end
        catch
            disp(['Doesn''t find this object: ''',tmp_name,'''']);
        end
    end
    
    tmp_name=['./SPO_VACorCHL',num2str(loop_i),'.ply'];
    prompt = ['Mesh of ',tmp_name,' seems not OK? 1=true, 0/enter=false: '];
    tmp_str = input(prompt);
    if isempty(tmp_str)
        tmp_str = 0;
    end
    if tmp_str~=0
        TMP_count_fail_4=TMP_count_fail_4+1;
        SAVE_FAIL_LIST_ext_SPO_VACCHL(TMP_count_fail_4)=loop_i;
    end
    %savefig(['Fig_SPO_MS_',num2str(loop_i),'.fig']);close(gcf)
end

%%%% 1. lower epidermis
%SAVE_FAIL_LIST_ext_EPL_l=SAVE_FAIL_LIST_EPL_l;
TMP_count_fail=size(SAVE_FAIL_LIST_EPL_l,1);
for loop_i=1:numel(GLB_tag_EPL_l_set)
    if GLB_filedir_mode==1
        tmp_name=['../nonMS/EPL_l_',num2str(loop_i),'.ply'];
    else
        tmp_name=['./EPL_l_',num2str(loop_i),'.ply'];
    end
    try
        %[tri,pts]=ply_read_xy(tmp_name);
        [tri,pts]=ply_read(tmp_name,'tri');
        axes(ax1);%subplot(1,2,1)%cell level
        cla(ax1)
        trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facecolor',[0.5,0.5,0.5],'facealpha',alpha,'EdgeColor',linecolor);axis equal;
        axes(ax2);%subplot(1,2,2)%obj level
        trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facecolor',[0.5,0.5,0.5],'facealpha',alpha,'EdgeColor',linecolor);axis equal;
        %interactive manually check by visualization
        prompt = ['Mesh of ',tmp_name,' seems not OK? 1=true, 0/enter=false: '];
        tmp_str = input(prompt);
        if isempty(tmp_str)
            tmp_str = 0;
        end
        if tmp_str~=0
            TMP_count_fail=TMP_count_fail+1;
            SAVE_FAIL_LIST_ext_EPL_l(TMP_count_fail)=loop_i;
        end
    catch
        disp(['Doesn''t find this object: ''',tmp_name,'''']);
    end
end

if GLB_MSCHLcheck_mode==2 % otherwise no PAL second batch
%%%% PAL second batch
%%%% 3-5. PAL
%SAVE_FAIL_LIST_ext_PAL_VAC=SAVE_FAIL_LIST_PAL_VAC;
TMP_count_fail_2=size(SAVE_FAIL_LIST_PAL_VAC,1);
%SAVE_FAIL_LIST_ext_PAL_CHL=SAVE_FAIL_LIST_PAL_CHL;
TMP_count_fail_3=size(SAVE_FAIL_LIST_PAL_CHL,1);
for tmp_loop_i=1:numel(SAVE_FAIL_LIST_ext_PAL_VACCHL)
    loop_i=SAVE_FAIL_LIST_ext_PAL_VACCHL(tmp_loop_i);
%     %%%% 3. PAL MS
%     if GLB_filedir_mode==1
%         tmp_name=['../MS/PAL_MS_',num2str(loop_i),'.ply'];
%     else
%         tmp_name=['./PAL_MS_',num2str(loop_i),'.ply'];
%     end
%     try
%         %[tri,pts]=ply_read_xy(tmp_name);
%         [tri,pts]=ply_read(tmp_name,'tri');
%         axes(ax1);%subplot(1,2,1)
%         cla(ax1)
%         trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facecolor',[0.5,0.5,0.5],'facealpha',alpha,'EdgeColor',linecolor);hold on;axis equal;
%         axes(ax2);%subplot(1,2,2)
%         trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facecolor',[0.5,0.5,0.5],'facealpha',alpha,'EdgeColor',linecolor);axis equal;
%         %interactive manually check by visualization
%         prompt = ['Mesh of ',tmp_name,' seems not OK? 1=true, 0/enter=false: '];
%         tmp_str = input(prompt);
%         if isempty(tmp_str)
%             tmp_str = 0;
%         end
%         if tmp_str~=0
%             TMP_count_fail_1=TMP_count_fail_1+1;
%             SAVE_FAIL_LIST_ext_PAL_MS(TMP_count_fail_1)=loop_i;
%         end
%     catch
%         disp(['Doesn''t find this object: ''',tmp_name,'''']);
%     end
%     %close(gcf)
    %%%% 4. PAL VAC
    if GLB_filedir_mode==1
        tmp_name=['../MS/PAL_VAC_',num2str(loop_i),'.ply'];
    else
        tmp_name=['./PAL_VAC_',num2str(loop_i),'.ply'];
    end
    try
        %[tri,pts]=ply_read_xy(tmp_name);
        [tri,pts]=ply_read(tmp_name,'tri');
        axes(ax1);%subplot(1,2,1)
        cla(ax1)
        trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facecolor',[0.2,0.2,0.8],'facealpha',alpha,'EdgeColor',linecolor);hold on;axis equal;
        axes(ax2);%subplot(1,2,2)
        trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facecolor',[0.2,0.2,0.8],'facealpha',alpha,'EdgeColor',linecolor);axis equal;
        %interactive manually check by visualization
        prompt = ['Mesh of ',tmp_name,' seems not OK? 1=true, 0/enter=false: '];
        tmp_str = input(prompt);
        if isempty(tmp_str)
            tmp_str = 0;
        end
        if tmp_str~=0
            TMP_count_fail_2=TMP_count_fail_2+1;
            SAVE_FAIL_LIST_ext_PAL_VAC(TMP_count_fail_2)=loop_i;
        end
    catch
        disp(['Doesn''t find this object: ''',tmp_name,'''']);
    end
    
    %%%% 5. PAL CHL
    tmptag_set_PAL_CHL=GLB_tag_PAL_CHL{loop_i};
    for loop_i_CHL=1:numel(tmptag_set_PAL_CHL)
        if GLB_filedir_mode==1
            tmp_name=['../MS/PAL_MS_',num2str(loop_i),'_CHL_',num2str(loop_i_CHL),'.ply'];
        else
            tmp_name=['./PAL_MS_',num2str(loop_i),'_CHL_',num2str(loop_i_CHL),'.ply'];
        end
        try
            %[tri,pts]=ply_read_xy(tmp_name);
            [tri,pts]=ply_read(tmp_name,'tri');
            axes(ax1);%subplot(1,2,1)
            trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facecolor',[0.2,0.8,0.2],'facealpha',alpha,'EdgeColor',linecolor);hold on;axis equal;
            axes(ax2);%subplot(1,2,2)
            trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facecolor',[0.2,0.8,0.2],'facealpha',alpha,'EdgeColor',linecolor);axis equal;
            %interactive manually check by visualization
            prompt = ['Mesh of ',tmp_name,' seems not OK? 1=true, 0/enter=false: '];
            tmp_str = input(prompt);
            if isempty(tmp_str)
                tmp_str = 0;
            end
            if tmp_str~=0
                TMP_count_fail_3=TMP_count_fail_3+1;
                SAVE_FAIL_LIST_ext_PAL_CHL(TMP_count_fail_3,1)=loop_i;
                SAVE_FAIL_LIST_ext_PAL_CHL(TMP_count_fail_3,2)=loop_i_CHL;
            end
        catch
            disp(['Doesn''t find this object: ''',tmp_name,'''']);
        end
    end
    
    %savefig(['Fig_PAL_MS_',num2str(loop_i),'.fig']);close(gcf)
end
end

%%%% SPO second batch
%%%% 6-8. SPO
%SAVE_FAIL_LIST_ext_SPO_VAC=SAVE_FAIL_LIST_SPO_VAC;
TMP_count_fail_2=size(SAVE_FAIL_LIST_SPO_VAC,1);
%SAVE_FAIL_LIST_ext_SPO_CHL=SAVE_FAIL_LIST_SPO_CHL;
TMP_count_fail_3=size(SAVE_FAIL_LIST_SPO_CHL,1);
for tmp_loop_i=1:numel(SAVE_FAIL_LIST_ext_SPO_VACCHL)
    loop_i=SAVE_FAIL_LIST_ext_SPO_VACCHL(tmp_loop_i);
%     %%%% 6. SPO MS
%     if GLB_filedir_mode==1
%         tmp_name=['../MS/SPO_MS_',num2str(loop_i),'.ply'];
%     else
%         tmp_name=['./SPO_MS_',num2str(loop_i),'.ply'];
%     end
%     try
%         %[tri,pts]=ply_read_xy(tmp_name);
%         [tri,pts]=ply_read(tmp_name,'tri');
%         %axes(ax1);%subplot(1,3,1)
%         %trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facecolor',[0.5,0.5,0.5],'facealpha',alpha,'EdgeColor',linecolor);hold on;axis equal;
%         axes(ax2);%subplot(1,3,2)
%         cla(ax2)
%         trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facecolor',[0.5,0.5,0.5],'facealpha',alpha,'EdgeColor',linecolor);hold on;axis equal;
%         axes(ax3);%subplot(1,3,3)
%         trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facecolor',[0.5,0.5,0.5],'facealpha',alpha,'EdgeColor',linecolor);axis equal;
%         %interactive manually check by visualization
%         prompt = 'Mesh seems not OK?: 1=true, 0/enter=false ';
%         tmp_str = input(prompt);
%         if isempty(tmp_str)
%             tmp_str = 0;
%         end
%         if tmp_str~=0
%             TMP_count_fail_1=TMP_count_fail_1+1;
%             SAVE_FAIL_LIST_ext_SPO_MS(TMP_count_fail_1)=loop_i;
%         end
%     catch
%         disp(['Doesn''t find this object: ''',tmp_name,'''']);
%     end
%     %close(gcf)
    %%%% 7. SPO VAC
    if GLB_filedir_mode==1
        tmp_name=['../MS/SPO_VAC_',num2str(loop_i),'.ply'];
    else
        tmp_name=['./SPO_VAC_',num2str(loop_i),'.ply'];
    end
    try
        %[tri,pts]=ply_read_xy(tmp_name);
        [tri,pts]=ply_read(tmp_name,'tri');
        axes(ax1);%subplot(1,2,1)
        cla(ax1)
        trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facecolor',[0.2,0.2,0.8],'facealpha',alpha,'EdgeColor',linecolor);hold on;axis equal;
        axes(ax2);%subplot(1,2,2)
        trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facecolor',[0.2,0.2,0.8],'facealpha',alpha,'EdgeColor',linecolor);axis equal;
        %interactive manually check by visualization
        prompt = 'Mesh seems not OK?: 1=true, 0/enter=false ';
        tmp_str = input(prompt);
        if isempty(tmp_str)
            tmp_str = 0;
        end
        if tmp_str~=0
            TMP_count_fail_2=TMP_count_fail_2+1;
            SAVE_FAIL_LIST_ext_SPO_VAC(TMP_count_fail_2)=loop_i;
        end
    catch
        disp(['Doesn''t find this object: ''',tmp_name,'''']);
    end
    
    %%%% 8. SPO CHL
    tmptag_set_SPO_CHL=GLB_tag_SPO_CHL{loop_i};
    for loop_i_CHL=1:numel(tmptag_set_SPO_CHL)
        if GLB_filedir_mode==1
            tmp_name=['../MS/SPO_MS_',num2str(loop_i),'_CHL_',num2str(loop_i_CHL),'.ply'];
        else
            tmp_name=['./SPO_MS_',num2str(loop_i),'_CHL_',num2str(loop_i_CHL),'.ply'];
        end
        try
            %[tri,pts]=ply_read_xy(tmp_name);
            [tri,pts]=ply_read(tmp_name,'tri');
            axes(ax1);%subplot(1,2,1)
            trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facecolor',[0.2,0.8,0.2],'facealpha',alpha,'EdgeColor',linecolor);hold on;axis equal;
            axes(ax2);%subplot(1,2,2)
            trisurf(tri,pts(:,1),pts(:,2),pts(:,3),'facecolor',[0.2,0.8,0.2],'facealpha',alpha,'EdgeColor',linecolor);axis equal;
            %interactive manually check by visualization
            prompt = 'Mesh seems not OK?: 1=true, 0/enter=false ';
            tmp_str = input(prompt);
            if isempty(tmp_str)
                tmp_str = 0;
            end
            if tmp_str~=0
                TMP_count_fail_3=TMP_count_fail_3+1;
                SAVE_FAIL_LIST_ext_SPO_CHL(TMP_count_fail_3,1)=loop_i;
                SAVE_FAIL_LIST_ext_SPO_CHL(TMP_count_fail_3,2)=loop_i_CHL;
            end
        catch
            disp(['Doesn''t find this object: ''',tmp_name,'''']);
        end
    end
    
    %savefig(['Fig_SPO_MS_',num2str(loop_i),'.fig']);close(gcf)
end

if GLB_manualcheck_save_mode==1
    %if GLB_MSCHLcheck_mode==1
        save SAVE_RT_geo_export_aftervisualcheck_full.mat -regexp '^(?!(model|ans)$).'
    %elseif GLB_MSCHLcheck_mode==2
        %save SAVE_RT_geo_export_aftervisualcheck_perMS.mat -regexp '^(?!(model|ans)$).'
    %end
elseif GLB_manualcheck_save_mode==2
    %if GLB_MSCHLcheck_mode==1
        save TMP_SAVE_RT_geo_export_aftervisualcheck_full.mat -regexp '^(?!(model|ans)$).'
    %elseif GLB_MSCHLcheck_mode==2
        %save TMP_SAVE_RT_geo_export_aftervisualcheck_perMS.mat -regexp '^(?!(model|ans)$).'
    %end
end
