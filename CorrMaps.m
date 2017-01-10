function varargout = CorrMaps(varargin)
% CORRMAPS MATLAB code for CorrMaps.fig
%      CORRMAPS, by itself, creates a new CORRMAPS or raises the existing
%      singleton*.
%
%      H = CORRMAPS returns the handle to a new CORRMAPS or the handle to
%      the existing singleton*.
%
%      CORRMAPS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CORRMAPS.M with the given input arguments.
%
%      CORRMAPS('Property','Value',...) creates a new CORRMAPS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CorrMaps_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CorrMaps_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CorrMaps

% Last Modified by GUIDE v2.5 09-Jan-2017 10:14:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CorrMaps_OpeningFcn, ...
                   'gui_OutputFcn',  @CorrMaps_OutputFcn, ...
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


% --- Executes just before CorrMaps is made visible.
function CorrMaps_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CorrMaps (see VARARGIN)

% Choose default command line output for CorrMaps
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CorrMaps wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CorrMaps_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function p12_Callback(hObject, eventdata, handles)
% hObject    handle to p12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p12 as text
%        str2double(get(hObject,'String')) returns contents of p12 as a double


% --- Executes during object creation, after setting all properties.
function p12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p13_Callback(hObject, eventdata, handles)
% hObject    handle to p13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p13 as text
%        str2double(get(hObject,'String')) returns contents of p13 as a double


% --- Executes during object creation, after setting all properties.
function p13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p23_Callback(hObject, eventdata, handles)
% hObject    handle to p23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of p23 as text
%        str2double(get(hObject,'String')) returns contents of p23 as a double


% --- Executes during object creation, after setting all properties.
function p23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function s12_Callback(hObject, eventdata, handles)
% hObject    handle to s12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of s12 as text
%        str2double(get(hObject,'String')) returns contents of s12 as a double


% --- Executes during object creation, after setting all properties.
function s12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to s12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function s13_Callback(hObject, eventdata, handles)
% hObject    handle to s13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of s13 as text
%        str2double(get(hObject,'String')) returns contents of s13 as a double


% --- Executes during object creation, after setting all properties.
function s13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to s13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function s23_Callback(hObject, eventdata, handles)
% hObject    handle to s23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of s23 as text
%        str2double(get(hObject,'String')) returns contents of s23 as a double


% --- Executes during object creation, after setting all properties.
function s23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to s23 (see GCBO)
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

function d1_Callback(hObject, eventdata, handles)
% hObject    handle to d1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of d1 as text
%        str2double(get(hObject,'String')) returns contents of d1 as a double


% --- Executes during object creation, after setting all properties.
function d1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to d1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function d2_Callback(hObject, eventdata, handles)
% hObject    handle to d2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of d2 as text
%        str2double(get(hObject,'String')) returns contents of d2 as a double


% --- Executes during object creation, after setting all properties.
function d2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to d2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function d3_Callback(hObject, eventdata, handles)
% hObject    handle to d3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of d3 as text
%        str2double(get(hObject,'String')) returns contents of d3 as a double


% --- Executes during object creation, after setting all properties.
function d3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to d3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function template_Callback(hObject, eventdata, handles)
% hObject    handle to template (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of template as text
%        str2double(get(hObject,'String')) returns contents of template as a double


% --- Executes during object creation, after setting all properties.
function template_CreateFcn(hObject, eventdata, handles)
% hObject    handle to template (see GCBO)
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
runs=str2double(get(handles.nr_runs,'String')); 

disp('...loads image dimensions..');
temp_img = get(handles.template,'String');
temp_img=load_nii(temp_img);
dim = size(temp_img.img);
x = dim(1);
y = dim(2);
z = dim(3);

name1 = get(handles.p12,'String');
name2 = get(handles.p23,'String');
name3 = get(handles.p13,'String');

name4 = get(handles.s12,'String');
name5 = get(handles.s23,'String');
name6 = get(handles.s13,'String');

disp('...please choose output directory...');
dir_results=spm_select(1,'dir','choose results directory');
cd (dir_results); 

% 1 vs 2
disp('...creates correlation maps...');

 if runs>1
     
first = get(handles.d1,'String');
first = load_nii(first);
first = first.img;
first (~first) = nan;

second = get(handles.d2,'String');
second = load_nii(second);
second = second.img;
second (~second) = nan;

% Pearson
% create correlation vector 
for ind_x = 1:x
    for ind_y = 1:y
        for ind_z = 1:z
            first_voxel = first (ind_x, ind_y, ind_z, :);
            second_voxel = second (ind_x, ind_y, ind_z, :);
            r = corrcoef(first_voxel,second_voxel, 'rows', 'pairwise');
            if isnan (r(1,2))
                r_vec_pear_1_2(ind_x, ind_y, ind_z) = 0;
            else
                r_vec_pear_1_2(ind_x, ind_y, ind_z) = r(1,2);
            end;
        end;
    end;
end;
           
% save correlation map

target_img = temp_img;
target_img.fileprefix = name1;
target_img.img = r_vec_pear_1_2;

save_nii(target_img,target_img.fileprefix);

% Spearman
% create correlation vector 
for ind_x = 1:x
    for ind_y = 1:y
        for ind_z = 1:z
            first_voxel = squeeze(first (ind_x, ind_y, ind_z, :));
            second_voxel = squeeze(second (ind_x, ind_y, ind_z, :));
            r = corr(first_voxel,second_voxel, 'Type','Spearman', 'rows', 'pairwise');
            r_vec_spea_1_2(ind_x, ind_y, ind_z) = r;
        end;
    end;
end;
           
% save correlation map

target_img = temp_img;
target_img.fileprefix = name4;
target_img.img = r_vec_spea_1_2;

save_nii(target_img,target_img.fileprefix);

end;
 
%% 2 vs 3
 if runs>2
     
second = get(handles.d2,'String');
second = load_nii(second);
second = second.img;
second (~second) = nan;

third = get(handles.d3,'String');
third = load_nii(third);
third = third.img;
third (~third) = nan;

% Pearson
% create correlation vector 
for ind_x = 1:x
    for ind_y = 1:y
        for ind_z = 1:z
            second_voxel = second (ind_x, ind_y, ind_z, :);
            third_voxel = third (ind_x, ind_y, ind_z, :);
            r = corrcoef(second_voxel,third_voxel, 'rows', 'pairwise');
            if isnan (r(1,2))
                r_vec_pear_2_3(ind_x, ind_y, ind_z) = 0;
            else
                r_vec_pear_2_3(ind_x, ind_y, ind_z) = r(1,2);
            end;
        end;
    end;
end;
           
% save correlation map

target_img = temp_img;
target_img.fileprefix = name2;
target_img.img = r_vec_pear_2_3;

save_nii(target_img,target_img.fileprefix);

% Spearman
% create correlation vector 
for ind_x = 1:x
    for ind_y = 1:y
        for ind_z = 1:z
            second_voxel = squeeze(second (ind_x, ind_y, ind_z, :));
            third_voxel = squeeze(third (ind_x, ind_y, ind_z, :));
            r = corr(second_voxel,third_voxel, 'Type','Spearman', 'rows', 'pairwise');
            r_vec_spea_2_3(ind_x, ind_y, ind_z) = r;
            
        end;
    end;
end;
           
% save correlation map

target_img = temp_img;
target_img.fileprefix = name5;
target_img.img = r_vec_spea_2_3;

save_nii(target_img,target_img.fileprefix);

 end;
 
 %% 1 vs 3
 if runs>2
     
first = get(handles.d1,'String');
first = load_nii(first);
first = first.img;
first (~first) = nan;

third = get(handles.d3,'String');
third = load_nii(third);
third = third.img;
third (~third) = nan;

% Pearson
% create correlation vector 
for ind_x = 1:x
    for ind_y = 1:y
        for ind_z = 1:z
            first_voxel = first (ind_x, ind_y, ind_z, :);
            third_voxel = third (ind_x, ind_y, ind_z, :);
            r = corrcoef(first_voxel,third_voxel, 'rows', 'pairwise');
            if isnan (r(1,2))
                r_vec_pear_1_3(ind_x, ind_y, ind_z) = 0;
            else
                r_vec_pear_1_3(ind_x, ind_y, ind_z) = r(1,2);
            end;
        end;
    end;
end;
           
% save correlation map

target_img = temp_img;
target_img.fileprefix = name3;
target_img.img = r_vec_pear_1_3;

save_nii(target_img,target_img.fileprefix);

% Spearman
% create correlation vector 
for ind_x = 1:x
    for ind_y = 1:y
        for ind_z = 1:z
            first_voxel = squeeze(first (ind_x, ind_y, ind_z, :));
            third_voxel = squeeze(third (ind_x, ind_y, ind_z, :));
            r = corr(first_voxel,third_voxel, 'Type','Spearman', 'rows', 'pairwise');
            r_vec_spea_1_3(ind_x, ind_y, ind_z) = r;
         end;
    end;
end;
           
% save correlation map

target_img = temp_img;
target_img.fileprefix = name6;
target_img.img = r_vec_spea_1_3;

save_nii(target_img,target_img.fileprefix);

end;

disp('...finished correlation maps')


% --- Executes on button press in show.
function show_Callback(hObject, eventdata, handles)
% hObject    handle to show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
struct=cellstr(spm_select(1,'nii','structural image'));
corr=cellstr(spm_select(1,'nii','correlation map'));
spm_mask(corr{1,1},struct{1,1},0.001);
masked=cellstr(spm_select(1,'nii','masked image "m*"'));
masked = load_nii(masked{1,1});
view_nii(masked);
