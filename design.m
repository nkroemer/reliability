function varargout = design(varargin)
% DESIGN MATLAB code for design.fig
%      DESIGN, by itself, creates a new DESIGN or raises the existing
%      singleton*.
%
%      H = DESIGN returns the handle to a new DESIGN or the handle to
%      the existing singleton*.
%
%      DESIGN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DESIGN.M with the given input arguments.
%
%      DESIGN('Property','Value',...) creates a new DESIGN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before design_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to design_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help design

% Last Modified by GUIDE v2.5 03-Feb-2017 09:20:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @design_OpeningFcn, ...
                   'gui_OutputFcn',  @design_OutputFcn, ...
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


% --- Executes just before design is made visible.
function design_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to design (see VARARGIN)

% Choose default command line output for design
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes design wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = design_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function subjects_Callback(hObject, eventdata, handles)
% hObject    handle to subjects (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of subjects as text
%        str2double(get(hObject,'String')) returns contents of subjects as a double


% --- Executes during object creation, after setting all properties.
function subjects_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subjects (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function runs_Callback(hObject, eventdata, handles)
% hObject    handle to runs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of runs as text
%        str2double(get(hObject,'String')) returns contents of runs as a double


% --- Executes during object creation, after setting all properties.
function runs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to runs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in list.
function list_Callback(hObject, eventdata, handles)
% hObject    handle to list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
list_subj = cellstr(spm_select(Inf,'mat','load list with subjects'));
assignin('base','list',list_subj);

% --- Executes on button press in dir.
function dir_Callback(hObject, eventdata, handles)
% hObject    handle to dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dir_results = cellstr(spm_select(Inf,'dir','choose results directory'));
assignin('base','dir',dir_results);


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



function stats_Callback(hObject, eventdata, handles)
% hObject    handle to stats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stats as text
%        str2double(get(hObject,'String')) returns contents of stats as a double


% --- Executes during object creation, after setting all properties.
function stats_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in stats_dir.
function stats_dir_Callback(hObject, eventdata, handles)
% hObject    handle to stats_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in save.
dir_stats = cellstr(spm_select(Inf,'dir','choose stats directory'));
assignin('base','dir_stats',dir_stats);



function prefix_design_Callback(hObject, eventdata, handles)
% hObject    handle to prefix_design (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of prefix_design as text
%        str2double(get(hObject,'String')) returns contents of prefix_design as a double


% --- Executes during object creation, after setting all properties.
function prefix_design_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prefix_design (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
oldpointer = get(handles.figure1, 'pointer'); 
set(handles.figure1, 'pointer', 'watch') 
drawnow;

subjects = get(handles.subjects,'String');
runs = get(handles.runs,'String');
list=evalin('base','list');
dir=evalin('base','dir');
con=get(handles.con,'String');
stats = get(handles.stats,'String');
dir_stats=evalin('base','dir_stats');
load('template_3D-4D.mat');
prefix_design = get(handles.prefix_design,'String');
name_design = sprintf('%s_study_design.mat',prefix_design);



study_design = struct('number_subjects',subjects,'number_sessions',runs,'subject_list',list,'results_directory',dir,'contrast',con,'stats_directory',stats,'stats_path',dir_stats);
clearvars dir list;
assignin('base','study_design',study_design);
% %save study_design in results folder
% save(name_design,'study_design');

name ='4D';
batch_name = sprintf('batch_3Dto4D_%s.mat',con);

runs=str2double(study_design.number_sessions); 
nr_subj=str2double(study_design.number_subjects);
load(study_design.subject_list);
con=study_design.contrast;

disp('...load contrast images to create 4D images...');
for i = 1:runs
    cd (study_design.stats_path)
    con_list=cell(nr_subj,1);
            for j = 1:nr_subj
                cd (sprintf('%d',vp(j)));
                cd (sprintf(study_design.stats_directory,i));
                con_list{j,1}=sprintf('%s\\%s.nii,1',pwd,con);
                cd (study_design.stats_path);
            end;
    estr=sprintf('con_img_%d = con_list;',i);
    eval(estr);
end;

    % create batch
    disp('...creates batch...');
    for i = 1:runs
        img_name = sprintf('%s_%d.nii',name,i);
        matlabbatch{i}.spm.util.cat.name = img_name;
        if i == 1
        matlabbatch{i}.spm.util.cat.vols = con_img_1;
        elseif i == 2
        matlabbatch{i}.spm.util.cat.vols = con_img_2;
        elseif i == 3
        matlabbatch{i}.spm.util.cat.vols = con_img_3;
        end;
        matlabbatch{i}.spm.util.cat.dtype = 4; % default for data type (INT16 - signed short) ; 0 : same data type as input images
    end;

    % save batch
    save(batch_name,'matlabbatch');
   
    disp('...runs batch...');
    spm_jobman('serial',batch_name);

    % go to results folder 
    cd (study_design.results_directory);
        
    disp('...moves 4D files to results folder...');
    % move files
    for k = 1:runs
        stats=sprintf(study_design.stats_directory,k);
        first=sprintf('%s\\%d\\%s',study_design.stats_path,vp(1),stats);
        file1 = sprintf('%s\\%s_%d.nii',first,name,k);
        file2 = sprintf('%s\\%s_%d.mat',first,name,k);
        movefile (file1,study_design.results_directory,'f');
        movefile (file2,study_design.results_directory,'f');
    end;
    %save study design in results folder
    save(name_design,'study_design');


    set(handles.figure1, 'pointer', oldpointer)
    disp('...DONE');
    
