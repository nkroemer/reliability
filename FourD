function varargout = FourD(varargin)
% FOURD MATLAB code for FourD.fig
%      FOURD, by itself, creates a new FOURD or raises the existing
%      singleton*.
%
%      H = FOURD returns the handle to a new FOURD or the handle to
%      the existing singleton*.
%
%      FOURD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FOURD.M with the given input arguments.
%
%      FOURD('Property','Value',...) creates a new FOURD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FourD_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FourD_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FourD

% Last Modified by GUIDE v2.5 06-Jan-2017 10:33:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FourD_OpeningFcn, ...
                   'gui_OutputFcn',  @FourD_OutputFcn, ...
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


% --- Executes just before FourD is made visible.
function FourD_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FourD (see VARARGIN)

% Choose default command line output for FourD
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FourD wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FourD_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function nr_subj_Callback(hObject, eventdata, handles)
% hObject    handle to nr_subj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nr_subj as text
%        str2double(get(hObject,'String')) returns contents of nr_subj as a double



% --- Executes during object creation, after setting all properties.
function nr_subj_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nr_subj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function nr_runs_Callback(hObject, eventdata, handles)
% hObject    handle to nr_runs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nr_runs as text
%        str2double(get(hObject,'String')) returns contents of nr_runs as a double


% --- Executes during object creation, after setting all properties.
function nr_runs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nr_runs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function batchname_Callback(hObject, eventdata, handles)
% hObject    handle to batchname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of batchname as text
%        str2double(get(hObject,'String')) returns contents of batchname as a double


% --- Executes during object creation, after setting all properties.
function batchname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to batchname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function filename_Callback(hObject, eventdata, handles)
% hObject    handle to filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filename as text
%        str2double(get(hObject,'String')) returns contents of filename as a double


% --- Executes during object creation, after setting all properties.
function filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load ('template_3D-4D.mat');
load('vp_gesamt.mat');
n=str2double(get(handles.nr_subj,'String')); 
runs=str2double(get(handles.nr_runs,'String')); 
name = get(handles.filename,'String');
batch_name = get(handles.batchname,'String');

disp('...please choose contrast images...');
for i = 1:runs
        if i == 1
        con_img_1 = cellstr(spm_select(n,'nii','contrast images run 1'));
            for j = 1:n
                con_img_1{j,1} = sprintf('%s,1',con_img_1{j,1});
            end;
        elseif i == 2
        con_img_2 = cellstr(spm_select(n,'nii','contrast images run 2'));
            for j = 1:n
                con_img_2{j,1} = sprintf('%s,1',con_img_2{j,1});
            end;
        elseif i == 3
        con_img_3= cellstr(spm_select(n,'nii','contrast images run 3'));
            for j = 1:n
                con_img_3{j,1} = sprintf('%s,1',con_img_3{j,1});
            end;
        end;
    end;

    % create batch
    disp('...creates batch...');
    for i = 1:runs
        img_name = sprintf('%d_%s.nii',i,name);
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
    disp('...please choose output directory...');
    dir_results=spm_select(1,'dir','choose results directory');
    cd (dir_results);
    
    
    disp('...moves 4D files to results folder...');
    % move files
    for k = 1:runs
        img_name = sprintf('%d_%s',k,name);
        file1 = sprintf('M:\\SeSyN\\019\\BMBF_itech\\Juliane\\niftii\\%06d\\stats_%d_maskfix\\%s.nii',vp(1),k,img_name);
        file2 = sprintf('M:\\SeSyN\\019\\BMBF_itech\\Juliane\\niftii\\%06d\\stats_%d_maskfix\\%s.mat',vp(1),k,img_name);
        movefile (file1,dir_results,'f');
        movefile (file2,dir_results,'f');
       
        if k == 1
            img_name = sprintf('%d_%s.nii',k,name);
            first_4d = load_nii (img_name);
        elseif k == 2
            img_name = sprintf('%d_%s.nii',k,name);
            second_4d = load_nii (img_name);
        elseif k == 3
            img_name = sprintf('%d_%s.nii',k,name);
            third_4d = load_nii (img_name);
        end;

    end;
    disp('...DONE');
    
