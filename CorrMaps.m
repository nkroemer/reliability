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

% Last Modified by GUIDE v2.5 08-Feb-2017 15:28:31

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

% --- Executes on button press in design.
function design_Callback(hObject, eventdata, handles)
% hObject    handle to design (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
study_design = cellstr(spm_select(1,'mat','load study design'));
load(study_design{1});
assignin('base','study_design',study_design);

% --- Executes on button press in split.
function split_Callback(hObject, eventdata, handles)
% hObject    handle to split (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of split


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

% --- Executes on button press in pearson.
function pearson_Callback(hObject, eventdata, handles)
% hObject    handle to pearson (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pearson


% --- Executes on button press in spearman.
function spearman_Callback(hObject, eventdata, handles)
% hObject    handle to spearman (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of spearman


% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
oldpointer = get(handles.figure1, 'pointer'); 
set(handles.figure1, 'pointer', 'watch') 
drawnow;

study_design=evalin('base','study_design');
runs=str2double(study_design.number_sessions); 
nr_subj=str2double(study_design.number_subjects);
load(study_design.subject_list);
con=study_design.contrast;
stats=study_design.stats_directory;
path=study_design.stats_path;

pear = get(handles.pearson,'value');
spea = get(handles.spearman,'value');
split = get(handles.split,'value');

disp('...loads image dimensions..');
stats_temp =sprintf(stats,1);
temp_img = sprintf('%s\\%d\\%s\\%s.nii',path,vp(1),stats_temp,con);
temp_img=load_nii(temp_img);
dim = size(temp_img.img);
x = dim(1);
y = dim(2);
z = dim(3);

name = get(handles.prefix,'String');

dir_results=study_design.results_directory;
cd (dir_results); 


if split == 0
    % 1 vs 2
    disp('...creates correlation maps...');

     if runs>1

    one = load_nii('4D_1.nii');
    one = one.img;
    one (~one) = nan;

    second = load_nii('4D_2.nii');
    second = second.img;
    second (~second) = nan;

    % create correlation vectors 
    r_vec_pear_1_2 = zeros(x,y,z);
    r_vec_spea_1_2 = zeros(x,y,z);
    z_r_vec_pear_1_2 = zeros(x,y,z);
    z_r_vec_spea_1_2 = zeros(x,y,z);
    
    for ind_x = 1:x
        for ind_y = 1:y
            for ind_z = 1:z
                first_voxel = one (ind_x, ind_y, ind_z, :);
                second_voxel = second (ind_x, ind_y, ind_z, :);
                % Pearson
                if pear==1
                    r = corrcoef(first_voxel,second_voxel, 'rows', 'pairwise');
                    if isnan (r(1,2))
                        r_vec_pear_1_2(ind_x, ind_y, ind_z) = 0;
                        z_r_vec_pear_1_2(ind_x, ind_y, ind_z) = 0;
                    else
                        r_vec_pear_1_2(ind_x, ind_y, ind_z) = r(1,2);
                        z_r_vec_pear_1_2(ind_x, ind_y, ind_z) = atanh(r(1,2));
                    end;
                end;
                %Spearman
                if spea==1
                    first_voxel = squeeze(one (ind_x, ind_y, ind_z, :));
                    second_voxel = squeeze(second (ind_x, ind_y, ind_z, :));
                  r = corr(first_voxel,second_voxel, 'Type','Spearman', 'rows', 'pairwise');
                  r_vec_spea_1_2(ind_x, ind_y, ind_z) = r;  
                  z_r_vec_spea_1_2(ind_x, ind_y, ind_z) = atanh(r);                    
                end;
            end;
        end;
    end;
                    % save correlation map
                    target_img = temp_img;
                    file = sprintf('%s_pear_1_2.nii',name);
                    target_img.fileprefix = file;
                    target_img.img = r_vec_pear_1_2;
                    save_nii(target_img,target_img.fileprefix);  
                    
                    target_img = temp_img;
                    file = sprintf('z_%s_pear_1_2.nii',name);
                    target_img.fileprefix = file;
                    target_img.img = z_r_vec_pear_1_2;
                    save_nii(target_img,target_img.fileprefix);                     

                    target_img = temp_img;
                    file = sprintf('%s_spea_1_2.nii',name);
                    target_img.fileprefix = file;
                    target_img.img = r_vec_spea_1_2;
                    save_nii(target_img,target_img.fileprefix);
                    
                    target_img = temp_img;
                    file = sprintf('z_%s_spea_1_2.nii',name);
                    target_img.fileprefix = file;
                    target_img.img = z_r_vec_spea_1_2;
                    save_nii(target_img,target_img.fileprefix);                   
     end; 
    %% 1 vs 3
    if runs>2
    third = load_nii('4D_3.nii');
    third = third.img;
    third (~third) = nan;

    % create correlation vectors 
    r_vec_pear_1_3 = zeros(x,y,z);
    r_vec_spea_1_3 = zeros(x,y,z);
    z_r_vec_pear_1_3 = zeros(x,y,z);
    z_r_vec_spea_1_3 = zeros(x,y,z);
    
    for ind_x = 1:x
        for ind_y = 1:y
            for ind_z = 1:z
                first_voxel = one (ind_x, ind_y, ind_z, :);
                third_voxel = third (ind_x, ind_y, ind_z, :);
                % Pearson
                if pear==1
                    r = corrcoef(first_voxel,third_voxel, 'rows', 'pairwise');
                    if isnan (r(1,2))
                        r_vec_pear_1_3(ind_x, ind_y, ind_z) = 0;
                        z_r_vec_pear_1_3(ind_x, ind_y, ind_z) = 0;
                       
                    else
                        r_vec_pear_1_3(ind_x, ind_y, ind_z) = r(1,2);
                        z_r_vec_pear_1_3(ind_x, ind_y, ind_z) = atanh(r(1,2));
                   
                    end;

                end;
                %Spearman
                if spea==1
                  first_voxel = squeeze(one (ind_x, ind_y, ind_z, :));
                  third_voxel = squeeze(third (ind_x, ind_y, ind_z, :));
                  r = corr(first_voxel,third_voxel, 'Type','Spearman', 'rows', 'pairwise');
                  r_vec_spea_1_3(ind_x, ind_y, ind_z) = r;  
                  z_r_vec_spea_1_3(ind_x, ind_y, ind_z) = atanh(r);  
                  

                end;
            end;
        end;

    end;
                    % save correlation map
                    target_img = temp_img;
                    file = sprintf('%s_pear_1_3.nii',name);
                    target_img.fileprefix = file;
                    target_img.img = r_vec_pear_1_3;
                    save_nii(target_img,target_img.fileprefix);
                    target_img = temp_img;
                    file = sprintf('%s_spea_1_3.nii',name);
                    target_img.fileprefix = file;
                    target_img.img = r_vec_spea_1_3;
                    save_nii(target_img,target_img.fileprefix);
                    
                    target_img = temp_img;
                    file = sprintf('z_%s_pear_1_3.nii',name);
                    target_img.fileprefix = file;
                    target_img.img = z_r_vec_pear_1_3;
                    save_nii(target_img,target_img.fileprefix);
                    target_img = temp_img;
                    file = sprintf('z_%s_spea_1_3.nii',name);
                    target_img.fileprefix = file;
                    target_img.img = r_vec_spea_1_3;
                    save_nii(target_img,target_img.fileprefix);

     % 2 vs 3
    % create correlation vectors 
    r_vec_pear_2_3 = zeros(x,y,z);
    r_vec_spea_2_3 = zeros(x,y,z);
    z_r_vec_pear_2_3 = zeros(x,y,z);
    z_r_vec_spea_2_3 = zeros(x,y,z);

    for ind_x = 1:x
        for ind_y = 1:y
            for ind_z = 1:z
                second_voxel = second (ind_x, ind_y, ind_z, :);
                third_voxel = third (ind_x, ind_y, ind_z, :);
                % Pearson
                if pear==1
                    r = corrcoef(second_voxel,third_voxel, 'rows', 'pairwise');
                    if isnan (r(1,2))
                        r_vec_pear_2_3(ind_x, ind_y, ind_z) = 0;
                        z_r_vec_pear_2_3(ind_x, ind_y, ind_z) = 0;

                    else
                        r_vec_pear_2_3(ind_x, ind_y, ind_z) = r(1,2);
                        z_r_vec_pear_2_3(ind_x, ind_y, ind_z) = atanh(r(1,2));

                    end;

                end;
                %Spearman
                if spea==1
                    second_voxel = squeeze(second (ind_x, ind_y, ind_z, :));
                    third_voxel = squeeze(third (ind_x, ind_y, ind_z, :));
                  r = corr(second_voxel,third_voxel, 'Type','Spearman', 'rows', 'pairwise');
                  r_vec_spea_2_3(ind_x, ind_y, ind_z) = r;  
                  z_r_vec_spea_2_3(ind_x, ind_y, ind_z) = atanh(r);  

                end;
            end;
        end;
    end;     
                    % save correlation map
                    target_img = temp_img;
                    file = sprintf('%s_pear_2_3.nii',name);
                    target_img.fileprefix = file;
                    target_img.img = r_vec_pear_2_3;
                    save_nii(target_img,target_img.fileprefix); 
                    target_img = temp_img;
                    file = sprintf('%s_spea_2_3.nii',name);
                    target_img.fileprefix = file;
                    target_img.img = r_vec_spea_2_3;
                    save_nii(target_img,target_img.fileprefix);
                    
                    target_img = temp_img;
                    file = sprintf('z_%s_pear_2_3.nii',name);
                    target_img.fileprefix = file;
                    target_img.img = z_r_vec_pear_2_3;
                    save_nii(target_img,target_img.fileprefix); 
                    target_img = temp_img;
                    file = sprintf('z_%s_spea_2_3.nii',name);
                    target_img.fileprefix = file;
                    target_img.img = z_r_vec_spea_2_3;
                    save_nii(target_img,target_img.fileprefix);
    end;

    if runs>3
    four = load_nii('4D_4.nii');
    four = four.img;
    four (~four) = nan;
    % 1 vs 4
    % create correlation vectors 
    r_vec_pear_1_4 = zeros(x,y,z);
    r_vec_spea_1_4 = zeros(x,y,z);
    z_r_vec_pear_1_4 = zeros(x,y,z);
    z_r_vec_spea_1_4 = zeros(x,y,z);
    
    for ind_x = 1:x
        for ind_y = 1:y
            for ind_z = 1:z
                first_voxel = one (ind_x, ind_y, ind_z, :);
                four_voxel = four (ind_x, ind_y, ind_z, :);
                % Pearson
                if pear==1
                    r = corrcoef(first_voxel,four_voxel, 'rows', 'pairwise');
                    if isnan (r(1,2))
                        r_vec_pear_1_4(ind_x, ind_y, ind_z) = 0;
                        z_r_vec_pear_1_4(ind_x, ind_y, ind_z) = 0;
                        
                    else
                        r_vec_pear_1_4(ind_x, ind_y, ind_z) = r(1,2);
                        z_r_vec_pear_1_4(ind_x, ind_y, ind_z) = atanh(r(1,2));
                       
                    end;

                end;
                %Spearman
                if spea==1
                first_voxel = squeeze(one (ind_x, ind_y, ind_z, :));
                four_voxel = squeeze(four (ind_x, ind_y, ind_z, :));
                  r = corr(first_voxel,four_voxel, 'Type','Spearman', 'rows', 'pairwise');
                  r_vec_spea_1_4(ind_x, ind_y, ind_z) = r;  
                  z_r_vec_spea_1_4(ind_x, ind_y, ind_z) = atanh(r);  
                 

                end;
            end;
        end;
    end;
                    % save correlation map
                    target_img = temp_img;
                    file = sprintf('%s_pear_1_4.nii',name);
                    target_img.fileprefix = file;
                    target_img.img = r_vec_pear_1_4;
                    save_nii(target_img,target_img.fileprefix);  
                    target_img = temp_img;
                    file = sprintf('%s_spea_1_4.nii',name);
                    target_img.fileprefix = file;
                    target_img.img = r_vec_spea_1_4;
                    save_nii(target_img,target_img.fileprefix);
                    
                    target_img = temp_img;
                    file = sprintf('z_%s_pear_1_4.nii',name);
                    target_img.fileprefix = file;
                    target_img.img = z_r_vec_pear_1_4;
                    save_nii(target_img,target_img.fileprefix);  
                    target_img = temp_img;
                    file = sprintf('z_%s_spea_1_4.nii',name);
                    target_img.fileprefix = file;
                    target_img.img = z_r_vec_spea_1_4;
                    save_nii(target_img,target_img.fileprefix);
    % 2 vs 4
    % create correlation vectors 
    r_vec_pear_2_4 = zeros(x,y,z);
    r_vec_spea_2_4 = zeros(x,y,z);
    z_r_vec_pear_2_4 = zeros(x,y,z);
    z_r_vec_spea_2_4 = zeros(x,y,z);
    
    for ind_x = 1:x
        for ind_y = 1:y
            for ind_z = 1:z
                second_voxel = second (ind_x, ind_y, ind_z, :);
                four_voxel = four (ind_x, ind_y, ind_z, :);
                % Pearson
                if pear==1
                    r = corrcoef(second_voxel,four_voxel, 'rows', 'pairwise');
                    if isnan (r(1,2))
                        r_vec_pear_2_4(ind_x, ind_y, ind_z) = 0;
                        z_r_vec_pear_2_4(ind_x, ind_y, ind_z) = 0;

                    else
                        r_vec_pear_2_4(ind_x, ind_y, ind_z) = r(1,2);
                        z_r_vec_pear_2_4(ind_x, ind_y, ind_z) = atanh(r(1,2));

                    end;

                end;
                %Spearman
                if spea==1
                second_voxel = squeeze(second (ind_x, ind_y, ind_z, :));
                four_voxel = squeeze(four (ind_x, ind_y, ind_z, :));
                  r = corr(second_voxel,four_voxel, 'Type','Spearman', 'rows', 'pairwise');
                  r_vec_spea_2_4(ind_x, ind_y, ind_z) = r;  
                  z_r_vec_spea_2_4(ind_x, ind_y, ind_z) = atanh(r);  
                  

                end;
            end;
        end;
    end;   
                    % save correlation map
                    target_img = temp_img;
                    file = sprintf('%s_pear_2_4.nii',name);
                    target_img.fileprefix = file;
                    target_img.img = r_vec_pear_2_4;
                    save_nii(target_img,target_img.fileprefix);  
                    target_img = temp_img;
                    file = sprintf('%s_spea_2_4.nii',name);
                    target_img.fileprefix = file;
                    target_img.img = r_vec_spea_2_4;
                    save_nii(target_img,target_img.fileprefix);
                    
                    target_img = temp_img;
                    file = sprintf('z_%s_pear_2_4.nii',name);
                    target_img.fileprefix = file;
                    target_img.img = z_r_vec_pear_2_4;
                    save_nii(target_img,target_img.fileprefix);  
                    target_img = temp_img;
                    file = sprintf('z_%s_spea_2_4.nii',name);
                    target_img.fileprefix = file;
                    target_img.img = z_r_vec_spea_2_4;
                    save_nii(target_img,target_img.fileprefix);

    % 3 vs 4
    % create correlation vectors 
    r_vec_pear_3_4 = zeros(x,y,z);
    r_vec_spea_3_4 = zeros(x,y,z);
    z_r_vec_pear_3_4 = zeros(x,y,z);
    z_r_vec_spea_3_4 = zeros(x,y,z);
    
    for ind_x = 1:x
        for ind_y = 1:y
            for ind_z = 1:z
                third_voxel = third (ind_x, ind_y, ind_z, :);
                four_voxel = four (ind_x, ind_y, ind_z, :);
                % Pearson
                if pear==1
                    r = corrcoef(third_voxel,four_voxel, 'rows', 'pairwise');
                    if isnan (r(1,2))
                        r_vec_pear_3_4(ind_x, ind_y, ind_z) = 0;
                        z_r_vec_pear_3_4(ind_x, ind_y, ind_z) = 0;

                    else
                        r_vec_pear_3_4(ind_x, ind_y, ind_z) = r(1,2);
                        z_r_vec_pear_3_4(ind_x, ind_y, ind_z) = atanh(r(1,2));

                    end;

                end;
                %Spearman
                if spea==1
                third_voxel = squeeze(third (ind_x, ind_y, ind_z, :));
                four_voxel = squeeze(four (ind_x, ind_y, ind_z, :));
                  r = corr(third_voxel,four_voxel, 'Type','Spearman', 'rows', 'pairwise');
                  r_vec_spea_3_4(ind_x, ind_y, ind_z) = r;  
                  z_r_vec_spea_3_4(ind_x, ind_y, ind_z) = atanh(r);  

                end;
            end;
        end;
    end;  
                    % save correlation map
                    target_img = temp_img;
                    file = sprintf('%s_pear_3_4.nii',name);
                    target_img.fileprefix = file;
                    target_img.img = r_vec_pear_3_4;
                    save_nii(target_img,target_img.fileprefix);   
                    target_img = temp_img;
                    file = sprintf('%s_spea_3_4.nii',name);
                    target_img.fileprefix = file;
                    target_img.img = r_vec_spea_3_4;
                    save_nii(target_img,target_img.fileprefix);
                    
                    target_img = temp_img;
                    file = sprintf('z_%s_pear_3_4.nii',name);
                    target_img.fileprefix = file;
                    target_img.img = z_r_vec_pear_3_4;
                    save_nii(target_img,target_img.fileprefix);   
                    target_img = temp_img;
                    file = sprintf('z_%s_spea_3_4.nii',name);
                    target_img.fileprefix = file;
                    target_img.img = z_r_vec_spea_3_4;
                    save_nii(target_img,target_img.fileprefix);
    end;

    if runs>4
    msg='Toolbox only allows 4 runs. Please adapt the script or contact the developers.';
    error(msg);
    end;
    
else
    for i = 1:runs

    img1=sprintf('4D_%d_split1.nii',i);
    one = load_nii(img1);
    one = one.img;
    one (~one) = nan;
    
    img2=sprintf('4D_%d_split2.nii',i);
    second = load_nii(img2);
    second = second.img;
    second (~second) = nan;

    % create correlation vectors 
    r_vec_pear_1_2 = zeros(x,y,z);
    r_vec_spea_1_2 = zeros(x,y,z);
    z_r_vec_pear_1_2 = zeros(x,y,z);
    z_r_vec_spea_1_2 = zeros(x,y,z);
    
    for ind_x = 1:x
        for ind_y = 1:y
            for ind_z = 1:z
                first_voxel = one (ind_x, ind_y, ind_z, :);
                second_voxel = second (ind_x, ind_y, ind_z, :);
                % Pearson
                if pear==1
                    r = corrcoef(first_voxel,second_voxel, 'rows', 'pairwise');
                    if isnan (r(1,2))
                        r_vec_pear_1_2(ind_x, ind_y, ind_z) = 0;
                        z_r_vec_pear_1_2(ind_x, ind_y, ind_z) = 0;
                        
                    else
                        r_vec_pear_1_2(ind_x, ind_y, ind_z) = r(1,2);
                        z_r_vec_pear_1_2(ind_x, ind_y, ind_z) = atanh(r(1,2));
                       
                    end;

                end;
                %Spearman
                if spea==1
                    first_voxel = squeeze(one (ind_x, ind_y, ind_z, :));
                    second_voxel = squeeze(second (ind_x, ind_y, ind_z, :));
                  r = corr(first_voxel,second_voxel, 'Type','Spearman', 'rows', 'pairwise');
                  r_vec_spea_1_2(ind_x, ind_y, ind_z) = r;  
                  z_r_vec_spea_1_2(ind_x, ind_y, ind_z) = atanh(r);  
                 
                end;
            end;
        end;
    end;
                    % save correlation map
                    target_img = temp_img;
                    file = sprintf('%s_pear_%d_split.nii',name,i);
                    target_img.fileprefix = file;
                    target_img.img = r_vec_pear_1_2;
                    save_nii(target_img,target_img.fileprefix);   

                    target_img = temp_img;
                    file = sprintf('%s_spea_%d_split.nii',name,i);
                    target_img.fileprefix = file;
                    target_img.img = r_vec_spea_1_2;
                    save_nii(target_img,target_img.fileprefix);
                    
                    target_img = temp_img;
                    file = sprintf('z_%s_pear_%d_split.nii',name,i);
                    target_img.fileprefix = file;
                    target_img.img = z_r_vec_pear_1_2;
                    save_nii(target_img,target_img.fileprefix);   

                    target_img = temp_img;
                    file = sprintf('z_%s_spea_%d_split.nii',name,i);
                    target_img.fileprefix = file;
                    target_img.img = z_r_vec_spea_1_2;
                    save_nii(target_img,target_img.fileprefix);
     end;  
disp('...finished correlation maps')
set(handles.figure1, 'pointer', oldpointer)
end;




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
