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

% Last Modified by GUIDE v2.5 03-Mar-2017 18:56:22

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
disp('starting creation of contrast definition');
box_path = pwd;
%get study design info
study_design=evalin('base','study_design');
runs = str2double(study_design.number_sessions);
nr_subj=str2double(study_design.number_subjects);
load(study_design.subject_list);
stats=study_design.stats_directory;
path=study_design.stats_path;
dir_results=study_design.results_directory;
name_design = study_design.name_design;

prefix = get(handles.prefix,'String');
name_contrast = sprintf('%s_contrast',prefix);

if runs == 1
    single_run = str2double(study_design.identifier_session);
end;

% initiate contrast_def structure
contrast_def = struct();

%get contrast_def information and extend study_design
two_cons = get(handles.two_cons,'value'); % 1 = comparison of two contrasts out of one statistic
contrast_def.two_contrasts = two_cons;
if runs == 1
    stats_filled = sprintf(stats,single_run);
else
    stats_filled = sprintf(stats,1);
end;


if two_cons == 0
    con=get(handles.con,'String');
    con_count = str2double(get(handles.con_count,'String'));
    file = sprintf('%s\\%s\\%s\\%s*',path,vp{1},stats_filled,con);
    con = dir(file);
    if isstruct(con)
        con = con(2).name;
        contrast_def.contrast = con;
    end;
    contrast_def.number_contrast = con_count;    
    [pathstr,name,ext] = fileparts(con);
    contrast_def.contrast_format = ext;

else
    con1=get(handles.con1,'String');
    con1_count = str2double(get(handles.con1_count,'String'));
    file = sprintf('%s\\%s\\%s\\%s*',path,vp{1},stats_filled,con1);
    con1 = dir(file);
    if isstruct(con1)
        con1 = con1(2).name;
        contrast_def.contrast1 = con1;
    end;
    contrast_def.number_contrast1 = con1_count;
    
    con2=get(handles.con2,'String');
    con2_count = str2double(get(handles.con2_count,'String'));
    file = sprintf('%s\\%s\\%s\\%s*',path,vp{1},stats_filled,con2);
    con2 = dir(file);
    if isstruct(con2)
        con2 = con2(2).name;
        contrast_def.contrast2 = con2;
    end;
    contrast_def.number_contrast2 = con2_count;
    [pathstr,name,ext] = fileparts(con2);
    contrast_def.contrast_format = ext;
end;

%% 3D to 4D conversion
load('template_3D-4D.mat'); % template for 3D to 4D conversion
name ='4D';


if two_cons == 0
    % list contrast_def images
    disp('...load contrast images to create 4D images...');
    if runs == 1
        cd (study_design.stats_path)
        con_list=cell(nr_subj,1);
        for j = 1:nr_subj
            cd (sprintf('%s',vp{j}));
            cd (sprintf(study_design.stats_directory,single_run));
            con_list{j,1}=sprintf('%s\\%s%s,1',pwd,con);
            cd (study_design.stats_path);
        end;
        eval(sprintf('con_img_%d = con_list;',single_run));
    else
        for i = 1:runs
            cd (study_design.stats_path)
            con_list=cell(nr_subj,1);
            for j = 1:nr_subj
                cd (sprintf('%s',vp{j}));
                cd (sprintf(study_design.stats_directory,i));
                con_list{j,1}=sprintf('%s\\%s,1',pwd,con);
                cd (study_design.stats_path);
            end;
            eval(sprintf('con_img_%d = con_list;',i));
        end;
    end;

    % create batch
    disp('...creates batch...');
    if runs == 1
        img_name = sprintf('%s_%d.nii',name,single_run);
        matlabbatch{i}.spm.util.cat.name = img_name;
        eval(sprintf('matlabbatch{1}.spm.util.cat.vols = con_img_%d;',single_run));
        matlabbatch{i}.spm.util.cat.dtype = 4; % default for data type (INT16 - signed short) ; 0 : same data type as input images
    else
        for i = 1:runs
            img_name = sprintf('%s_%d.nii',name,i);
            matlabbatch{i}.spm.util.cat.name = img_name;
            eval(sprintf('matlabbatch{i}.spm.util.cat.vols = con_img_%d;',i));
            matlabbatch{i}.spm.util.cat.dtype = 4; % default for data type (INT16 - signed short) ; 0 : same data type as input images
        end;
    end;

    % save batch
    cd(dir_results);
    save('batch_3D-4D.mat','matlabbatch');
    
    % run batch
    disp('...runs batch...');
    spm_jobman('serial','batch_3D-4D.mat');

    % go to results folder 
    cd (study_design.results_directory);

    disp('...moves 4D files to results folder...');
    % move files
    if runs == 1
        stats=sprintf(study_design.stats_directory,single_run);
        first=sprintf('%s\\%s\\%s',study_design.stats_path,vp{1},stats);
        file1 = sprintf('%s\\%s_%d.nii',first,name,single_run);
        file2 = sprintf('%s\\%s_%d.mat',first,name,single_run);
        movefile (file1,study_design.results_directory,'f');
        movefile (file2,study_design.results_directory,'f');       
    else
        for k = 1:runs
            stats=sprintf(study_design.stats_directory,k);
            first=sprintf('%s\\%s\\%s',study_design.stats_path,vp{1},stats);
            file1 = sprintf('%s\\%s_%d.nii',first,name,k);
            file2 = sprintf('%s\\%s_%d.mat',first,name,k);
            movefile (file1,study_design.results_directory,'f');
            movefile (file2,study_design.results_directory,'f');
        end;
    end;
    %save study design in results folder
    save(name_design,'study_design');

else
    disp('...load contrast images to create 4D images...');
    for ind_con = 1:2
    eval(sprintf('con=contrast_def.contrast%d;',ind_con));
        if runs == 1
            cd (study_design.stats_path)
            con_list=cell(nr_subj,1);
            for j = 1:nr_subj
                cd (sprintf('%s',vp{j}));
                cd (sprintf(study_design.stats_directory,single_run));
                con_list{j,1}=sprintf('%s\\%s,1',pwd,con);
                cd (study_design.stats_path);
            end;
            eval(sprintf('con_img_%d = con_list;',single_run));
        else
            for i = 1:runs
                cd (study_design.stats_path)
                con_list=cell(nr_subj,1);
                for j = 1:nr_subj
                    cd (sprintf('%s',vp{j}));
                    cd (sprintf(study_design.stats_directory,i));
                    con_list{j,1}=sprintf('%s\\%s,1',pwd,con);
                    cd (study_design.stats_path);
                 end;
            eval(sprintf('con_img_%d = con_list;',i));
            end;
        end;
                
            % create batch
            disp('...creates batch...');
            if runs == 1
                    img_name = sprintf('%s_%s_%d.nii',name,con,single_run);
                    matlabbatch{1}.spm.util.cat.name = img_name;
                    eval(sprintf('matlabbatch{1}.spm.util.cat.vols = con_img_%d;',single_run));
                    matlabbatch{1}.spm.util.cat.dtype = 4; % default for data type (INT16 - signed short) ; 0 : same data type as input images
            else
                for i = 1:runs
                    img_name = sprintf('%s_%s_%d.nii',name,con,i);
                    matlabbatch{i}.spm.util.cat.name = img_name;
                    eval(sprintf('matlabbatch{i}.spm.util.cat.vols = con_img_%d;',i));
                    matlabbatch{i}.spm.util.cat.dtype = 4; % default for data type (INT16 - signed short) ; 0 : same data type as input images
                end;
            end;
            
            % save batch
            save(batch_name,'matlabbatch');

            % run batch
            disp('...runs batch...');
            spm_jobman('serial',batch_name);

            % go to results folder 
            cd (study_design.results_directory);

            disp('...moves 4D files to results folder...');
            % move files
            if runs == 1
                stats=sprintf(study_design.stats_directory,single_run);
                first=sprintf('%s\\%s\\%s',study_design.stats_path,vp{1},stats);
                file1 = sprintf('%s\\%s_%s_%d.nii',first,name,con,single_run);
                file2 = sprintf('%s\\%s_%s_%d.mat',first,name,con,single_run);
                movefile (file1,study_design.results_directory,'f');
                movefile (file2,study_design.results_directory,'f');
            else
                for k = 1:runs
                    stats=sprintf(study_design.stats_directory,k);
                    first=sprintf('%s\\%s\\%s',study_design.stats_path,vp{1},stats);
                    file1 = sprintf('%s\\%s_%s_%d.nii',first,name,con,k);
                    file2 = sprintf('%s\\%s_%s_%d.mat',first,name,con,k);
                    movefile (file1,study_design.results_directory,'f');
                    movefile (file2,study_design.results_directory,'f');
                end;
            end;
   
        end;


end;
cd(study_design.results_directory)
%save contrast_def information
save(name_contrast,'contrast_def');
disp('...DONE');
cd(box_path);    
