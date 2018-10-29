function varargout = Contrast_Def(varargin)
% CONTRAST_DEF MATLAB code for Contrast_Def.fig
%      CONTRAST_DEF, by itself, creates a new CONTRAST_DEF or raises the existing
%      singleton*.
%
%      H = CONTRAST_DEF returns the handle to a new CONTRAST_DEF or the handle to
%      the existing singleton*.
%
%      CONTRAST_DEF('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CONTRAST_DEF.M with the given input arguments.
%
%      CONTRAST_DEF('Property','Value',...) creates a new CONTRAST_DEF or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Contrast_Def_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Contrast_Def_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Contrast_Def

% Last Modified by GUIDE v2.5 10-Aug-2018 10:52:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Contrast_Def_OpeningFcn, ...
                   'gui_OutputFcn',  @Contrast_Def_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Contrast_Def is made visible.
function Contrast_Def_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Contrast_Def (see VARARGIN)

% Choose default command line output for Contrast_Def
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Contrast_Def wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Contrast_Def_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function con1_Callback(hObject, eventdata, handles)
% hObject    handle to con1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of con1 as text
%        str2double(get(hObject,'String')) returns contents of con1 as a double


% --- Executes during object creation, after setting all properties.
function con1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to con1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function con2_Callback(hObject, eventdata, handles)
% hObject    handle to con2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of con2 as text
%        str2double(get(hObject,'String')) returns contents of con2 as a double


% --- Executes during object creation, after setting all properties.
function con2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to con2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function con1_count_Callback(hObject, eventdata, handles)
% hObject    handle to con1_count (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of con1_count as text
%        str2double(get(hObject,'String')) returns contents of con1_count as a double


% --- Executes during object creation, after setting all properties.
function con1_count_CreateFcn(hObject, eventdata, handles)
% hObject    handle to con1_count (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function con2_count_Callback(hObject, eventdata, handles)
% hObject    handle to con2_count (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of con2_count as text
%        str2double(get(hObject,'String')) returns contents of con2_count as a double


% --- Executes during object creation, after setting all properties.
function con2_count_CreateFcn(hObject, eventdata, handles)
% hObject    handle to con2_count (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in two_cons.
function two_cons_Callback(hObject, eventdata, handles)
% hObject    handle to two_cons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of two_cons



function con_Callback(hObject, eventdata, handles)
% hObject    handle to con (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of con as text
%        str2double(get(hObject,'String')) returns contents of con as a double


% --- Executes during object creation, after setting all properties.
function con_CreateFcn(hObject, eventdata, handles)
% hObject    handle to con (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function con_count_Callback(hObject, eventdata, handles)
% hObject    handle to con_count (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of con_count as text
%        str2double(get(hObject,'String')) returns contents of con_count as a double


% --- Executes during object creation, after setting all properties.
function con_count_CreateFcn(hObject, eventdata, handles)
% hObject    handle to con_count (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in design.
function design_Callback(hObject, eventdata, handles)
% hObject    handle to design (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in run.
study_design = cellstr(spm_select(1,'mat','load study design'));
load(study_design{1});
assignin('base','study_design',study_design);

function prefix_Callback(hObject, eventdata, handles)
% hObject    handle to prefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of prefix as text
%        str2double(get(hObject,'String')) returns contents of prefix as a double


% --- Executes during object creation, after setting all properties.
function prefix_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prefix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('Starting creation of contrast definition...');

%% define file seperator 
f = filesep;
box_path=evalin('base','box_path');

%% get GUI input
    % definition of contrast definition name
    prefix = get(handles.prefix,'String'); 
    name_contrast = sprintf('%s_contrast',prefix);
    two_cons = get(handles.two_cons,'value'); % 1 = comparison of two contrasts out of one statistic
    cond_1 = get(handles.cond1,'value'); % 1 = only one condition
    cond_2 = get(handles.cond2,'value'); % 1 = several conditions
    
%% get study design info
    study_design=evalin('base','study_design');
    runs = study_design.number_sessions; % number runs
    dir_results=study_design.results_directory; % result directory
    name_design = study_design.name_design;
    if runs == 1
        single_run = str2double(study_design.identifier_session);
    end;
    
    exStats = study_design.exist_stats;
    if exStats == 1
        nr_subj=str2double(study_design.number_subjects); % number subjects
        load(study_design.subject_list); % subject list 
        stats=study_design.stats_directory; % name of stats directory
        path=study_design.stats_path; % path to stats directories
    end;
    
    ex4D = study_design.exist_4D; 

%%
    % initiate contrast_def structure
    contrast_def = struct();
    
if exStats == 1
    %get contrast_def information and extend study_design by number of
    %parametric modulators
    contrast_def.two_contrasts = two_cons;
    if runs == 1
        stats_filled = sprintf(stats,single_run);
    else
        stats_filled = sprintf(stats,1);
    end;

    if two_cons == 0 
        con=get(handles.con,'String');      
        file = [path f id{1} f stats_filled f con '*']; 
        con = dir(file);
        if length(con)>1
            con = con(2).name;
            contrast_def.contrast = con;
        else
            con = con(1).name;
            contrast_def.contrast = con;        
        end;
        
        spmmat = [path f id{1} f stats_filled f 'SPM.mat'];
        load(spmmat);        
        for i = 1:length(SPM.xCon)
            if strcmp({SPM.xCon(i).Vcon.fname}, {con})==1
                for ind_c = 1:length(SPM.xCon(i).c)
                    if SPM.xCon(i).c(ind_c) == 1
                        idx = ind_c;
                    end;
                end;
            end;
        end;
        contrast_def.number_regressor = idx;    
        [pathstr,name,ext] = fileparts(con);
        contrast_def.contrast_number = str2double(regexp(con,'[\d]+','match'));
        contrast_def.contrast_format = ext;
        
        % info parametric modulators
        if strcmp(SPM.Sess(1).U(idx).P(1).name,'none')
            study_design.number_parametric = 0;
            contrast_def.number_parametric = 0;

            nr_para = 0;
        else
            nr_para = size(SPM.Sess(1).U(idx).P,2);
            study_design.number_parametric = nr_para;
            contrast_def.number_parametric = nr_para;            
        end;
        cd(dir_results)
        save(name_design,'study_design');
    elseif two_cons == 1 
        con1=get(handles.con1,'String');
        contrast_def.contrast1 = con1;        
        con=con1;
        spmmat = [path f id{1} f stats_filled f 'SPM.mat'];
        load(spmmat);        
        for i = 1:length(SPM.xCon)
            if strcmp({SPM.xCon(i).Vcon.fname}, [con '.nii'])==1
                idx = find(SPM.xCon(i).c==1);
                ext = '.nii';
            elseif strcmp({SPM.xCon(i).Vcon.fname}, [con '.img'])==1
                idx = find(SPM.xCon(i).c==1);
                ext = '.img';
            end;
        end;
        contrast_def.number_regressor1 = idx;
        contrast_def.contrast_format = ext;
        contrast_def.contrast1_number = str2double(regexp(con1,'[\d]+','match'));

        % info parametric modulators 
        new_idx = 0;
        for count = 1:length(SPM.Sess.U)
            if new_idx == idx
                Uidx = count;
            end;
            if length(SPM.Sess.U(count).name)>1
                for count2 = 1:(length(SPM.Sess.U(count).name)-1)
                    if new_idx+count2 == idx
                        Uidx = count;
                    end;
                end; 
            end;    
            new_idx = new_idx + length(SPM.Sess.U(count).name);        
        end; 
        contrast_def.Uidx1 = Uidx;

        if strcmp(SPM.Sess(1).U(Uidx).P(1).name,'none')
            study_design.number_parametric1 = 0;
            nr_para1 = 0;
        else
            nr_para1 = size(SPM.Sess(1).U(Uidx).P,2);
            study_design.number_parametric1 = nr_para1;
        end;
        cd(dir_results);
        save(name_design,'study_design');

        con2=get(handles.con2,'String');
        contrast_def.contrast2 = con2;        
        con=con2;
        spmmat = [path f id{1} f stats_filled f 'SPM.mat'];
        load(spmmat);        
        for i = 1:length(SPM.xCon)
            if strcmp({SPM.xCon(i).Vcon.fname}, [con '.nii'])==1
                idx = find(SPM.xCon(i).c==1);
            elseif strcmp({SPM.xCon(i).Vcon.fname}, [con '.img'])==1
                idx = find(SPM.xCon(i).c==1);
            end;
        end;       
        
        contrast_def.number_regressor2 = idx;
        contrast_def.contrast2_number = str2double(regexp(con2,'[\d]+','match'));

        % info parametric modulators 
        new_idx = 0;
        for count = 1:length(SPM.Sess.U)
            if new_idx == idx
                Uidx = count;
            end;
            if length(SPM.Sess.U(count).name)>1
                for count2 = 1:(length(SPM.Sess.U(count).name)-1)
                    if new_idx+count2 == idx
                        Uidx = count;
                    end;
                end; 
            end;    
            new_idx = new_idx + length(SPM.Sess.U(count).name);        
        end;
        contrast_def.Uidx2 = Uidx;
        
        if strcmp(SPM.Sess(1).U(Uidx).P(1).name,'none')
            study_design.number_parametric2 = 0;
            nr_para2 = 0;
        else
            nr_para2 = size(SPM.Sess(1).U(Uidx).P,2);
            study_design.number_parametric2 = nr_para2;
        end;
        cd(dir_results);
        save(name_design,'study_design');    
    end;

% 3D to 4D conversion
cd(box_path);
load(['scripts_templates' filesep 'template_3D-4D.mat']); % template for 3D to 4D conversion
if two_cons == 0 
   
    % list contrast_def images
    disp('...load contrast images to create 4D images...');
    if runs == 1
        cd(path);
        con_list=cell(nr_subj,1);
        for j = 1:nr_subj
            cd (sprintf('%s',id{j}));
            cd (sprintf(stats,single_run));
            con_list{j,1}=[pwd f con];
            cd(path);
        end;
        eval(sprintf('con_img_%d = con_list;',single_run));
        
        %not yet implemented
%         if nr_para > 0
%             for i_par = 1:nr_para
%                 cd(path)
%                 con_list=cell(nr_subj,1);
%                 for j = 1:nr_subj
%                     cd (sprintf('%s',id{j}));
%                     cd (sprintf(stats,single_run));
%                     con_list{j,1}=[pwd f 'con_' num2str(contrast_def.contrast_number+i_par,'%04d') '.nii,1'];
%                     cd (path);
%                 end;
%                 eval(sprintf('con_img_%d_par%d = con_list;',single_run,i_par));
%             end;
%         end;
        
    else
        for i = 1:runs
            cd(path)
            con_list=cell(nr_subj,1);
            for j = 1:nr_subj
                cd(sprintf('%s',id{j}));
                cd(sprintf(stats,i));
                con_list{j,1}=[pwd f con];
                cd(path);
            end;
            eval(sprintf('con_img_%d = con_list;',i));
        end;
        
        if nr_para > 0
            for i_par = 1:nr_para
                for i = 1:runs
                    cd(path)
                    con_list=cell(nr_subj,1);
                    for j = 1:nr_subj
                        cd (sprintf('%s',id{j}));
                        cd (sprintf(stats,i));
                        con_list{j,1}=[pwd f 'con_' num2str(contrast_def.contrast_number+i_par,'%04d') ext ',1'];
                        cd(path);
                    end;
                    eval(sprintf('con_img_%d_par%d = con_list;',i,i_par));
                end;
            end;
        end;       
        
    end;

    % create batch
    disp('...creates batch...');
    if runs == 1
        img_name = sprintf('4D_%d.nii',single_run);
        matlabbatch{1}.spm.util.cat.name = img_name;
        eval(sprintf('matlabbatch{1}.spm.util.cat.vols = con_img_%d;',single_run));
        matlabbatch{1}.spm.util.cat.dtype = 4; % default for data type (INT16 - signed short) ; 0 : same data type as input images
         if nr_para > 0
            for i_par = 1:nr_para
                img_name = sprintf('4D_par%d_%d.nii',i_par,single_run);
                matlabbatch{1+i_par}.spm.util.cat.name = img_name;
                eval(sprintf('matlabbatch{1+i_par}.spm.util.cat.vols = con_img_%d_par%d;',single_run,i_par));
                matlabbatch{1+i_par}.spm.util.cat.dtype = 4; % default for data type (INT16 - signed short) ; 0 : same data type as input images
            end;
         end;
   
    else
        for i = 1:runs
            img_name = sprintf('4D_%d.nii',i);
            matlabbatch{i}.spm.util.cat.name = img_name;
            eval(sprintf('matlabbatch{i}.spm.util.cat.vols = con_img_%d;',i));
            matlabbatch{i}.spm.util.cat.dtype = 4; % default for data type (INT16 - signed short) ; 0 : same data type as input images
        end;
         if nr_para > 0
            for i_par = 1:nr_para
                    idx=runs+((i_par-1)*runs)+1;
                for j = 1:runs
                    img_name = sprintf('4D_par%d_%d.nii',i_par,j);
                    matlabbatch{idx}.spm.util.cat.name = img_name;
                    eval(sprintf('matlabbatch{idx}.spm.util.cat.vols = con_img_%d_par%d;',j,i_par));
                    matlabbatch{idx}.spm.util.cat.dtype = 4; % default for data type (INT16 - signed short) ; 0 : same data type as input images
                    idx=idx+1;
                end;
            end;
         end;
    end;

    % save batch
    cd(dir_results);
    save('batch_3D-4D.mat','matlabbatch');
    
    % run batch
    disp('...runs batch...');
    spm_jobman('serial','batch_3D-4D.mat');

    % go to results folder 
    cd(study_design.results_directory);

    
    %save study design in results folder
    save(name_design,'study_design');

elseif two_cons == 1 
    disp('...load contrast images to create 4D images...');
        for ind_con = 1:2
            eval(sprintf('con=contrast_def.contrast%d;',ind_con));
            eval(sprintf('con_count=contrast_def.contrast%d_number;',ind_con));
            eval(sprintf('nr_para = nr_para%d;',ind_con));
            if runs == 1
                cd(path)
                con_list=cell(nr_subj,1);
                for j = 1:nr_subj
                    cd(sprintf('%s',id{j}));
                    cd(sprintf(stats,single_run));
                    con_list{j,1}=[pwd f con ext];
                    cd(path);
                end;
                eval(sprintf('con_img_%d = con_list;',single_run));
                if nr_para > 0
                    for i_par = 1:nr_para
                        cd (path)
                        con_list=cell(nr_subj,1);
                        for j = 1:nr_subj
                            cd(sprintf('%s',id{j}));
                            cd(sprintf(stats,single_run));
                            con_list{j,1}=[pwd f 'con_' sprintf('%04d',con_count+i_par) ext ',1'];
                            cd(path);
                        end;
                        eval(sprintf('con_img_%d_par%d = con_list;',single_run,i_par));
                    end;
                end;            
            else
                for i = 1:runs
                    cd(path)
                    con_list=cell(nr_subj,1);
                    for j = 1:nr_subj
                        cd (sprintf('%s',id{j}));
                        cd (sprintf(stats,i));
                        con_list{j,1}=[pwd f con ext ',1'];
                        cd (path);
                     end;
                    eval(sprintf('con_img_%d = con_list;',i));
                end;
               if nr_para > 0
                    for i_par = 1:nr_para
                        for i = 1:runs
                            cd (path)
                            con_list=cell(nr_subj,1);
                            for j = 1:nr_subj
                                cd (sprintf('%s',id{j}));
                                cd (sprintf(stats,i));
                                con_list{j,1}=[pwd f 'con_' sprintf('%04d',con_count+i_par) ext ',1'];
                                cd (path);
                            end;
                            eval(sprintf('con_img_%d_par%d = con_list;',i,i_par));
                        end;
                    end;
                end;       
            
         end;
                
            % create batch
            disp('...creates batch...');
            if runs == 1
                    img_name = sprintf('4D_%s_%d.nii',con,single_run);
                    matlabbatch{1}.spm.util.cat.name = img_name;
                    eval(sprintf('matlabbatch{1}.spm.util.cat.vols = con_img_%d;',single_run));
                    matlabbatch{1}.spm.util.cat.dtype = 4; % default for data type (INT16 - signed short) ; 0 : same data type as input images
                     if nr_para > 0
                        for i_par = 1:nr_para
                            img_name = sprintf('4D_%s_par%d_%d.nii',con,i_par,single_run);
                            matlabbatch{1+i_par}.spm.util.cat.name = img_name;
                            eval(sprintf('matlabbatch{1+i_par}.spm.util.cat.vols = con_img_par%d_%d;',i_par,single_run));
                            matlabbatch{1+i_par}.spm.util.cat.dtype = 4; % default for data type (INT16 - signed short) ; 0 : same data type as input images
                        end;
                     end;             
            else
                for i = 1:runs
                    img_name = sprintf('4D_%s_%d.nii',con,i);
                    matlabbatch{i}.spm.util.cat.name = img_name;
                    eval(sprintf('matlabbatch{i}.spm.util.cat.vols = con_img_%d;',i));
                    matlabbatch{i}.spm.util.cat.dtype = 4; % default for data type (INT16 - signed short) ; 0 : same data type as input images
                end;
                  if nr_para > 0
                    for i_par = 1:nr_para
                        idx=runs+1;
                        for j = 1:runs
                            img_name = sprintf('4D_%s_par%d_%d.nii',con,i_par,j);
                            matlabbatch{idx}.spm.util.cat.name = img_name;
                            eval(sprintf('matlabbatch{idx}.spm.util.cat.vols = con_img_%d_par%d;',j,i_par));
                            matlabbatch{idx}.spm.util.cat.dtype = 4; % default for data type (INT16 - signed short) ; 0 : same data type as input images
                            idx=idx+1;
                        end;
                    end;
                 end;               
            end;
            
            % save batch
            cd (dir_results);
            
            save('batch_3D-4D.mat','matlabbatch');

            % run batch
            disp('...runs batch...');
            spm_jobman('serial','batch_3D-4D.mat');

            % go to results folder 
            cd (dir_results);

    clear con_list con_img_1 con_img_2 matlabbatch
    end;
  
end;
    disp('...moves 4D files to results folder...');
    for i_run = 1:runs
        directory=[study_design.stats_path f id{1} f sprintf(stats,i_run)];
        cd(directory)
        list = dir('4D*');
        if ~isempty(list)
            for i_list = 1:length(list)
                movefile (list(i_list).name,dir_results,'f');
            end;
        end;
    end;

elseif ex4D == 1
    cd(dir_results)
    list = dir('4D*');
    example = load_nii(list(1).name);
    study_design.number_subjects = length(example.img(1,1,1,:));
    if cond_1 == 1
        contrast_def.number_conditions = 1; 
    elseif cond_2 == 1
        conditions = {};
        for i = 1:length(list)
            name = list(i).name;
            [a,b] = strtok(name,'_');
            [c,d] = strtok(b,'_');
            if ~any(strcmp(conditions,c))    
                conditions{end+1,1} = c;
            end;
        end;
        contrast_def.number_conditions = length(conditions);
        contrast_def.conditions = conditions;
    end;  
    
end;

cd(dir_results)
%save contrast_def information
save(name_contrast,'contrast_def');
save(name_design,'study_design');
disp('...DONE');
cd(box_path);    


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1
axes(hObject);
imshow('logo.png');


% --- Executes on button press in exist4D.
function exist4D_Callback(hObject, eventdata, handles)
% hObject    handle to exist4D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of exist4D


% --- Executes on button press in exist4D2.
function exist4D2_Callback(hObject, eventdata, handles)
% hObject    handle to exist4D2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of exist4D2


% --- Executes on button press in cond1.
function cond1_Callback(hObject, eventdata, handles)
% hObject    handle to cond1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cond1


% --- Executes on button press in cond2.
function cond2_Callback(hObject, eventdata, handles)
% hObject    handle to cond2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cond2
