function varargout = ICC_vox(varargin)
% ICC_VOX MATLAB code for ICC_vox.fig
%      ICC_VOX, by itself, creates a new ICC_VOX or raises the existing
%      singleton*.
%
%      H = ICC_VOX returns the handle to a new ICC_VOX or the handle to
%      the existing singleton*.
%
%      ICC_VOX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ICC_VOX.M with the given input arguments.
%
%      ICC_VOX('Property','Value',...) creates a new ICC_VOX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ICC_vox_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ICC_vox_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ICC_vox

% Last Modified by GUIDE v2.5 10-Apr-2017 16:26:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ICC_vox_OpeningFcn, ...
                   'gui_OutputFcn',  @ICC_vox_OutputFcn, ...
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


% --- Executes just before ICC_vox is made visible.
function ICC_vox_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ICC_vox (see VARARGIN)

% Choose default command line output for ICC_vox
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ICC_vox wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ICC_vox_OutputFcn(hObject, eventdata, handles) 
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

% --- Executes on button press in con.
function con_Callback(hObject, eventdata, handles)
% hObject    handle to con (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of con


% --- Executes on button press in abs.
function abs_Callback(hObject, eventdata, handles)
% hObject    handle to abs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of abs

% --- Executes on button press in split.
function split_Callback(hObject, eventdata, handles)
% hObject    handle to split (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of split

% --- Executes on button press in con_def.
function con_def_Callback(hObject, eventdata, handles)
% hObject    handle to con_def (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
con_def = cellstr(spm_select(1,'mat','load contrast definition'));
load(con_def{1});
assignin('base','contrast_def',contrast_def);

% --- Executes on button press in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('starting calculation of ICC voxelwise');

boxpath = pwd;

%% set parameters
%get study design information
study_design=evalin('base','study_design');
contrast_def=evalin('base','contrast_def');

runs=str2double(study_design.number_sessions); 
nr_subj=str2double(study_design.number_subjects);
load(study_design.subject_list);
stats=study_design.stats_directory;
path=study_design.stats_path;
dir_results = study_design.results_directory;
if runs == 1
    single_run = str2double(study_design.identifier_session);
end;

%get GUI input
cons = get(handles.con,'value');
abs = get(handles.abs,'value');
split = get(handles.split,'value');

%load contrast information
two_cons = contrast_def.two_contrasts;
if two_cons == 0
    con=contrast_def.contrast;
else
    con1=contrast_def.contrast1;
    con2=contrast_def.contrast2;
    con1_count=contrast_def.number_contrast1;
    con2_count=contrast_def.number_contrast2;
end;

cd(dir_results);
%% load 4D images
if split == 0 && two_cons == 0
    for i = 1:runs
        img = load_nii(sprintf('4D_%d.nii',i));
        eval(sprintf('img_%d = img;',i));
        eval(sprintf('img_%d = img_%d.img;',i,i));
    end;
disp('...loads image dimensions..');
stats_temp =sprintf(stats,1);
temp_img = sprintf('%s\\%s\\%s\\%s',path,vp{1},stats_temp,con);
temp_img=load_nii(temp_img);
dim = size(temp_img.img);
x = dim(1);
y = dim(2);
z = dim(3);

elseif split == 1
    nr_para = study_design.number_parametric;
    if runs == 1
        img1 = load_nii(sprintf('4D_%d_split1.nii',single_run));
        eval(sprintf('img_%d_split1 = img1;',single_run));
        eval(sprintf('img_%d_split1 = img_%d_split1.img;',single_run,single_run));
        img2 = load_nii(sprintf('4D_%d_split2.nii',single_run));
        eval(sprintf('img_%d_split2 = img2;',single_run));
        eval(sprintf('img_%d_split2 = img_%d_split2.img;',single_run,single_run));
    else
        for i = 1:runs    
            img1 = load_nii(sprintf('4D_%d_split1.nii',i));
            eval(sprintf('img_%d_split1 = img1;',i));
            eval(sprintf('img_%d_split1 = img_%d_split1.img;',i,i));
            img2 = load_nii(sprintf('4D_%d_split2.nii',i));
            eval(sprintf('img_%d_split2 = img2;',i));
            eval(sprintf('img_%d_split2 = img_%d_split2.img;',i,i));
        end;
    end;
    if nr_para > 0
        for ind_para = 1:nr_para
            if runs == 1
                nii = sprintf('4D_split1_par%d_%d.nii',ind_para,single_run);
                img1 = load_nii(nii);
                eval(sprintf('img_%d_par%d_split1 = img1;',single_run,ind_para));
                eval(sprintf('img_%d_par%d_split1 = img_%d_par%d_split1.img;',single_run,ind_para,single_run,ind_para));
                nii2 = sprintf('4D_split2_par%d_%d.nii',ind_para,single_run);
                img2 = load_nii(nii2);
                eval(sprintf('img_%d_par%d_split2 = img2;',single_run,ind_para));
                eval(sprintf('img_%d_par%d_split2 = img_%d_par%d_split2.img;',single_run,ind_para,single_run,ind_para));
            else
               for i = 1:runs    
                    nii = sprintf('4D_split1_par%d_%d.nii',ind_para,i);
                    img1 = load_nii(nii);
                    eval(sprintf('img_%d_par%d_split1 = img1;',i,ind_para));
                    eval(sprintf('img_%d_par%d_split1 = img_%d_par%d_split1.img;',i,ind_para,i,ind_para));
                    nii2 = sprintf('4D_split2_par%d_%d.nii',ind_para,i);
                    img2 = load_nii(nii2);
                    eval(sprintf('img_%d_par%d_split2 = img2;',i,ind_para));
                    eval(sprintf('img_%d_par%d_split2 = img_%d_par%d_split2.img;',i,ind_para,i,ind_para));
               end; 
            end;
        end;
    end;
    disp('...loads image dimensions..');
    stats_temp =sprintf(stats,1);
    temp_img = sprintf('%s\\%s\\%s\\%s',path,vp{1},stats_temp,con);
    temp_img=load_nii(temp_img);
    dim = size(temp_img.img);
    x = dim(1);
    y = dim(2);
    z = dim(3);
    
elseif two_cons == 1
    
    if runs == 1
        img = load_nii(sprintf('4D_%s_%d.nii',con1,single_run));
        eval(sprintf('img_%d_%d = img;',single_run,con1_count));
        eval(sprintf('img_%d_%d = img_%d_%d.img;',single_run,con1_count,single_run,con1_count));
        img2 = load_nii(sprintf('4D_%s_%d.nii',con2,single_run));
        eval(sprintf('img_%d_%d = img2;',single_run,con2_count));
        eval(sprintf('img_%d_%d = img_%d_%d.img;',single_run,con2_count,single_run,con2_count));
    else
        for i = 1:runs
            img = load_nii(sprintf('4D_%s_%d.nii',con1,i));
            eval(sprintf('img_%d_%d = img;',i,con1_count));
            eval(sprintf('img_%d_%d = img_%d_%d.img;',i,con1_count,i,con1_count));
            img2 = load_nii(sprintf('4D_%s_%d.nii',con2,i));
            eval(sprintf('img_%d_%d = img2;',i,con2_count));
            eval(sprintf('img_%d_%d = img_%d_%d.img;',i,con2_count,i,con2_count));
        end;
    end;
    disp('...loads image dimensions..');
    stats_temp =sprintf(stats,1);
    temp_img = sprintf('%s\\%s\\%s\\%s',path,vp{1},stats_temp,con1);
    temp_img=load_nii(temp_img);
    dim = size(temp_img.img);
    x = dim(1);
    y = dim(2);
    z = dim(3);
end;


%% initialization 
if split == 0 && two_cons == 0
    data=zeros(nr_subj,runs);
else
    data=zeros(nr_subj,2);
end;

ICC_con=zeros(x,y,z);
ICC_abs=zeros(x,y,z);
z_ICC_con=zeros(x,y,z);
z_ICC_abs=zeros(x,y,z);

% initialize summary
summary = [];
cols={};

%% calculating ICCs
disp('calculating ICCs')
if split == 0 && two_cons == 0
    for ind_x = 1:x
        fprintf('scanning x = %d \n', ind_x);
        for ind_y = 1:y
           for ind_z = 1:z
                for ind_run = 1:runs
                    for ind_subj = 1:nr_subj
                       eval(sprintf('img%d_voxel(ind_subj,1)= img_%d (ind_x, ind_y, ind_z, ind_subj);',ind_run,ind_run));
                    end;
                    eval(sprintf('data(:,ind_run)=img%d_voxel;',ind_run));
                end;
                
                %ICC
                nsamples=nr_subj*runs;

                grandmean=0;
                for sub=1:nr_subj,     
                    for sess=1:runs,
                       grandmean= grandmean + data(sub,sess);
                    end;
                end;
                grandmean=grandmean./nsamples;

                sessionmean=zeros(runs,1);
                for sess=1:runs
                    for sub=1:nr_subj,  
                        sessionmean(sess) = sessionmean(sess) + data(sub,sess);
                    end;
                    sessionmean(sess)=sessionmean(sess)./nr_subj;
                end;

                subjmean=zeros(nr_subj,1);
                for sub=1:nr_subj
                    for sess=1:runs
                        subjmean(sub)=subjmean(sub) + data(sub,sess);
                    end;
                    subjmean(sub)=subjmean(sub)./runs;
                end;

                % mean squares
                BMS=0; % between subject
                WMS=0; % within subject 
                EMS=0; % error
                JMS=0; % session

                for sub=1:nr_subj,    
                    BMS = BMS + (subjmean(sub)-grandmean).^2;
                    for sess=1:runs
                        WMS = WMS + (data(sub,sess)-subjmean(sub)).^2;
                        EMS = EMS + (data(sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
                    end
                end;

                for sess=1:runs
                    JMS = JMS + (sessionmean(sess)-grandmean).^2;
                end;

                %define the true value of the mean square.
                BMS= runs.*BMS./(nr_subj-1);
                WMS= WMS./(runs-1)./nr_subj;
                JMS= nr_subj.*JMS./(runs-1);
                EMS= EMS./(runs-1)./(nr_subj-1); 

                %consistency agreement  
                if cons==1
                    voxICC_con=(BMS-EMS)./(BMS+(runs-1).*WMS); 
                    ICC_con(ind_x, ind_y, ind_z) = voxICC_con;
                    z_ICC_con(ind_x, ind_y, ind_z) = .5.*log((1+voxICC_con)./(1-voxICC_con));
                end;

                %absolute agreement 
                if abs==1
                    voxICC_abs=(BMS-EMS)./(BMS+(runs-1).*EMS + runs.* (JMS-EMS)./nr_subj);
                    ICC_abs(ind_x, ind_y, ind_z) = voxICC_abs;
                    z_ICC_abs(ind_x,ind_y,ind_z) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                end;
           end; 
        end;
     end;

disp('saving ICC images')
% save ICC maps
target_img = temp_img;
target_img.fileprefix = 'ICC_con.nii';
target_img.img = ICC_con;
save_nii(target_img,target_img.fileprefix); 

clear target_img;
target_img = temp_img;
target_img.fileprefix = 'ICC_abs.nii';
target_img.img = ICC_abs;
save_nii(target_img,target_img.fileprefix); 

target_img = temp_img;
target_img.fileprefix = 'z_ICC_con.nii';
target_img.img = z_ICC_con;
save_nii(target_img,target_img.fileprefix); 

clear target_img;
target_img = temp_img;
target_img.fileprefix = 'z_ICC_abs.nii';
target_img.img = z_ICC_abs;
save_nii(target_img,target_img.fileprefix); 

%compute z mean ICC
mean_z = mean(z_ICC_con(:),'omitnan');
summary(1,end+1)=mean_z;
cols{1,end+1} = 'mean_z_ICC_con';
%compute min and max ICC
[my_min, idx] = min(z_ICC_con(:));
[my_max, idx] = max(z_ICC_con(:));
summary(1,end+1)= my_min;
cols{1,end+1} = 'min_z_ICC_con';
summary(1,end+1) = my_max;
cols{1,end+1} = 'max_z_ICC_con'; 

%compute z mean ICC
mean_z = mean(z_ICC_abs(:),'omitnan');
summary(1,end+1)=mean_z;
cols{1,end+1} = 'mean_z_ICC_abs';
%compute min and max ICC
[my_min, idx] = min(z_ICC_abs(:));
[my_max, idx] = max(z_ICC_abs(:));
summary(1,end+1)= my_min;
cols{1,end+1} = 'min_z_ICC_abs';
summary(1,end+1) = my_max;
cols{1,end+1} = 'max_z_ICC_abs'; 

%based on split-half
elseif split == 1
     for ind_x = 1:x
      fprintf('scanning x = %d \n', ind_x);
        for ind_y = 1:y
           for ind_z = 1:z
                for ind_run = 1:runs
                    if runs == 1
                        eval(sprintf('data_%d = zeros(nr_subj,2)',single_run));
                        for ind_subj = 1:nr_subj
                           eval(sprintf('img%d_voxel_split1(ind_subj,1)= img_%d_split1 (ind_x, ind_y, ind_z, ind_subj);',single_run,single_run));
                           eval(sprintf('img%d_voxel_split2(ind_subj,1)= img_%d_split2 (ind_x, ind_y, ind_z, ind_subj);',single_run,single_run));
                        end;
                        eval(sprintf('data_%d(:,1)=img%d_voxel_split1;',single_run));
                        eval(sprintf('data_%d(:,2)=img%d_voxel_split2;',single_run));                        
                    else
                        eval(sprintf('data_%d = zeros(nr_subj,2)',ind_run));
                        for ind_subj = 1:nr_subj
                           eval(sprintf('img%d_voxel_split1(ind_subj,1)= img_%d_split1 (ind_x, ind_y, ind_z, ind_subj);',ind_run,ind_run));
                           eval(sprintf('img%d_voxel_split2(ind_subj,1)= img_%d_split2 (ind_x, ind_y, ind_z, ind_subj);',ind_run,ind_run));
                        end;
                        eval(sprintf('data_%d(:,1)=img%d_voxel_split1;',ind_run));
                        eval(sprintf('data_%d(:,2)=img%d_voxel_split2;',ind_run));
                    end;
                end;
                
                %ICC
                for ind_run = 1:runs
                    if runs == 1
                        ind_run = single_run;
                    end;
                    nsamples=nr_subj*2;
                        
                    grandmean=0;
                    for sub=1:nr_subj,     
                        for sess=1:2,
                           eval(sprintf('grandmean= grandmean + data_%d(sub,sess)',ind_run));
                        end
                    end;
                    grandmean=grandmean./nsamples;

                    sessionmean=zeros(2,1);
                    for sess=1:2
                        for sub=1:nr_subj,  
                            eval(sprintf('sessionmean(sess) = sessionmean(sess) + data_%d(sub,sess)',ind_run));
                        end
                        sessionmean(sess)=sessionmean(sess)./nr_subj;
                    end

                    subjmean=zeros(nr_subj,1);
                    for sub=1:nr_subj
                        for sess=1:2
                            eval(sprintf('subjmean(sub)=subjmean(sub) + data_%d(sub,sess)',ind_run));
                        end
                        subjmean(sub)=subjmean(sub)./2;
                    end

                    % mean squares
                    BMS=0; % between subject
                    WMS=0; % within subject 
                    EMS=0; % error
                    JMS=0; % session

                    for sub=1:nr_subj,    
                        BMS = BMS + (subjmean(sub)-grandmean).^2;
                        for sess=1:2
                            eval(sprintf('WMS = WMS + (data_%d(sub,sess)-subjmean(sub)).^2',ind_run));
                            eval(sprintf('EMS = EMS + (data_%d(sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2',ind_run));
                        end
                    end;

                    for sess=1:2
                        JMS=  JMS + (sessionmean(sess)-grandmean).^2;
                    end;

                    %define the true value of the mean square.
                    BMS= 2.*BMS./(nr_subj-1);
                    WMS= WMS./(2-1)./nr_subj;
                    JMS= nr_subj.*JMS./(2-1);
                    EMS= EMS./(2-1)./(nr_subj-1); 

                    %consistency agreement  
                    if cons==1
                        voxICC_con=(BMS-EMS)./(BMS+(2-1).*WMS); 
                        ICC_con(ind_x, ind_y, ind_z) = voxICC_con;
                        z_ICC_con(ind_x, ind_y, ind_z) = .5.*log((1+voxICC_con)./(1-voxICC_con));
                    end;

                    %absolute agreement 
                    if abs==1
                        voxICC_abs=(BMS-EMS)./(BMS+(2-1).*EMS + 2.* (JMS-EMS)./nr_subj);
                        ICC_abs(ind_x, ind_y, ind_z) = voxICC_abs;
                        z_ICC_abs(ind_x,ind_y,ind_z) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                    end;
                    
            disp('saving ICC images')
             % save ICC maps
            target_img = temp_img;
            file1 = sprintf('ICC_con_%d_split.nii',ind_run);
            target_img.fileprefix = file1;
            target_img.img = ICC_con;
            save_nii(target_img,target_img.fileprefix); 

            clear target_img;
            target_img = temp_img;
            file2 = sprintf('ICC_abs_%d_split.nii',ind_run);
            target_img.fileprefix = file2;
            target_img.img = ICC_abs;
            save_nii(target_img,target_img.fileprefix); 

            target_img = temp_img;
            file3 = sprintf('z_ICC_con_%d_split.nii',ind_run);
            target_img.fileprefix = file3;
            target_img.img = z_ICC_con;
            save_nii(target_img,target_img.fileprefix); 

            clear target_img;
            target_img = temp_img;
            file4 = sprintf('z_ICC_abs_%d_split.nii',ind_run);
            target_img.fileprefix = file4;
            target_img.img = z_ICC_abs;
            save_nii(target_img,target_img.fileprefix); 

            %compute z mean ICC
            mean_z = mean(z_ICC_con(:),'omitnan');
            summary(1,end+1)=mean_z;
            cols{1,end+1} = sprintf('mean_z_ICC_con_split%d',ind_run);
            %compute min and max ICC
            [my_min, idx] = min(z_ICC_con(:));
            [my_max, idx] = max(z_ICC_con(:));
            summary(1,end+1)= my_min;
            cols{1,end+1} = sprintf('min_z_ICC_con_split%d',ind_run);
            summary(1,end+1) = my_max;
            cols{1,end+1} = sprintf('max_z_ICC_con_split%d',ind_run); 

            %compute z mean ICC
            mean_z = mean(z_ICC_abs(:),'omitnan');
            summary(1,end+1)=mean_z;
            cols{1,end+1} = sprintf('mean_z_ICC_abs_split%d',ind_run);
            %compute min and max ICC
            [my_min, idx] = min(z_ICC_abs(:));
            [my_max, idx] = max(z_ICC_abs(:));
            summary(1,end+1)= my_min;
            cols{1,end+1} = sprintf('min_z_ICC_abs_split%d',ind_run);
            summary(1,end+1) = my_max;
            cols{1,end+1} = sprintf('max_z_ICC_abs_split%d',ind_run); 
                end;
           end; 
        end;
     end;
     
    if nr_para > 0
        fprintf('investigating parametric contrast')
        for ind_para = 1:nr_para
             for ind_x = 1:x
               fprintf('scanning x = %d \n', ind_x);
                for ind_y = 1:y
                   for ind_z = 1:z
                        for ind_run = 1:runs
                            if runs == 1
                                ind_run = single_run;
                            end;
                            eval(sprintf('data_%d = zeros(nr_subj,2)',ind_run));
                            for ind_subj = 1:nr_subj
                               eval(sprintf('img%d_par%d_voxel_split1(ind_subj,1)= img_%d_par%d_split1 (ind_x, ind_y, ind_z, ind_subj);',ind_run,ind_para,ind_run,ind_para));
                               eval(sprintf('img%d_par%d_voxel_split2(ind_subj,1)= img_%d_par%d_split2 (ind_x, ind_y, ind_z, ind_subj);',ind_run,ind_para,ind_run,ind_para));
                            end;
                            eval(sprintf('data_%d(:,1)=img%d_par%d_voxel_split1;',ind_run,ind_para));
                            eval(sprintf('data_%d(:,2)=img%d_par%d_voxel_split2;',ind_run,ind_para));
                        end;
                        %ICC
                        for ind_run = 1:runs
                            if runs == 1
                                ind_run = single_run;
                            end;
                            nsamples=nr_subj*runs;

                            grandmean=0;
                            for sub=1:nr_subj,     
                                for sess=1:2,
                                   eval(sprintf('grandmean= grandmean + data_%d(sub,sess)',ind_run));
                                end
                            end;
                            grandmean=grandmean./nsamples;

                            sessionmean=zeros(2,1);
                            for sess=1:2
                                for sub=1:nr_subj,  
                                    eval(sprintf('sessionmean(sess) = sessionmean(sess) + data_%d(sub,sess)',ind_run));
                                end
                                sessionmean(sess)=sessionmean(sess)./nr_subj;
                            end

                            subjmean=zeros(nr_subj,1);
                            for sub=1:nr_subj
                                for sess=1:2
                                    eval(sprintf('subjmean(sub)=subjmean(sub) + data_%d(sub,sess)',ind_run));
                                end
                                subjmean(sub)=subjmean(sub)./2;
                            end

                            % mean squares
                            BMS=0; % between subject
                            WMS=0; % within subject 
                            EMS=0; % error
                            JMS=0; % session

                            for sub=1:nr_subj,    
                                BMS = BMS + (subjmean(sub)-grandmean).^2;
                                for sess=1:2
                                    eval(sprintf('WMS = WMS + (data_%d(sub,sess)-subjmean(sub)).^2',ind_run));
                                    eval(sprintf('EMS = EMS + (data_%d(sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2',ind_run));
                                end
                            end;

                            for sess=1:2
                                JMS=  JMS + (sessionmean(sess)-grandmean).^2;
                            end;

                            %define the true value of the mean square.
                            BMS= 2.*BMS./(nr_subj-1);
                            WMS= WMS./(2-1)./nr_subj;
                            JMS= nr_subj.*JMS./(2-1);
                            EMS= EMS./(2-1)./(nr_subj-1); 

                            %consistency agreement  
                            if cons==1
                                voxICC_con=(BMS-EMS)./(BMS+(2-1).*WMS); 
                                ICC_con(ind_x, ind_y, ind_z) = voxICC_con;
                                z_ICC_con(ind_x, ind_y, ind_z) = .5.*log((1+voxICC_con)./(1-voxICC_con));
                            end;

                            %absolute agreement 
                            if abs==1
                                voxICC_abs=(BMS-EMS)./(BMS+(2-1).*EMS + 2.* (JMS-EMS)./nr_subj);
                                ICC_abs(ind_x, ind_y, ind_z) = voxICC_abs;
                                z_ICC_abs(ind_x,ind_y,ind_z) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                            end;
                            
                            disp('saving ICC images')
                             % save ICC maps
                            target_img = temp_img;
                            file1 = sprintf('ICC_con_%d_split_par%d.nii',ind_run,ind_para);
                            target_img.fileprefix = file1;
                            target_img.img = ICC_con;
                            save_nii(target_img,target_img.fileprefix); 

                            clear target_img;
                            target_img = temp_img;
                            file2 = sprintf('ICC_abs_%d_split_par%d.nii',ind_run,ind_para);
                            target_img.fileprefix = file2;
                            target_img.img = ICC_abs;
                            save_nii(target_img,target_img.fileprefix); 

                            target_img = temp_img;
                            file3 = sprintf('z_ICC_con_%d_split_par%d.nii',ind_run,ind_para);
                            target_img.fileprefix = file3;
                            target_img.img = z_ICC_con;
                            save_nii(target_img,target_img.fileprefix); 

                            clear target_img;
                            target_img = temp_img;
                            file4 = sprintf('z_ICC_abs_%d_split_par%d.nii',ind_run,ind_para);
                            target_img.fileprefix = file4;
                            target_img.img = z_ICC_abs;
                            save_nii(target_img,target_img.fileprefix); 
                        
                            %compute z mean ICC
                            mean_z = mean(z_ICC_con(:),'omitnan');
                            summary(1,end+1)=mean_z;
                            cols{1,end+1} = sprintf('mean_z_ICC_con_split%d_par%d',ind_run,ind_para);
                            %compute min and max ICC
                            [my_min, idx] = min(z_ICC_con(:));
                            [my_max, idx] = max(z_ICC_con(:));
                            summary(1,end+1)= my_min;
                            cols{1,end+1} = sprintf('min_z_ICC_con_split%d_par%d',ind_run,ind_para);
                            summary(1,end+1) = my_max;
                            cols{1,end+1} = sprintf('max_z_ICC_con_split%d_par%d',ind_run,ind_para); 

                            %compute z mean ICC
                            mean_z = mean(z_ICC_abs(:),'omitnan');
                            summary(1,end+1)=mean_z;
                            cols{1,end+1} = sprintf('mean_z_ICC_abs_split%d_par%d',ind_run,ind_para);
                            %compute min and max ICC
                            [my_min, idx] = min(z_ICC_abs(:));
                            [my_max, idx] = max(z_ICC_abs(:));
                            summary(1,end+1)= my_min;
                            cols{1,end+1} = sprintf('min_z_ICC_abs_split%d_par%d',ind_run,ind_para);
                            summary(1,end+1) = my_max;
                            cols{1,end+1} = sprintf('max_z_ICC_abs_split%d_par%d',ind_run,ind_para);                                     
                        end;
                   end; 
                end;
             end; 
        end;    
    end;

% two contrasts out of one statistic
elseif two_cons == 1
for i_run = 1:runs
    if runs == 1
        i_run = single_run;
    end;
    for ind_x = 1:x
        fprintf('scanning x = %d \n', ind_x);
        for ind_y = 1:y
           for ind_z = 1:z
                for ind_con = 1:2
                    for ind_subj = 1:nr_subj
                       eval(sprintf('i_con = con%d_count;',ind_con));
                       eval(sprintf('img%d_%d_voxel(ind_subj,1)= img_%d_%d (ind_x, ind_y, ind_z, ind_subj);',i_run,i_con,i_run,i_con));
                    end;
                    eval(sprintf('data(:,ind_con)=img%d_%d_voxel;',i_run,i_con));
                end;
                
                %ICC
                nsamples=nr_subj*2;

                grandmean=0;
                for sub=1:nr_subj,     
                    for sess=1:2,
                       grandmean= grandmean + data(sub,sess);
                    end

                end;
                grandmean=grandmean./nsamples;

                sessionmean=zeros(2,1);
                for sess=1:2
                    for sub=1:nr_subj,  
                        sessionmean(sess) = sessionmean(sess) + data(sub,sess);
                    end
                    sessionmean(sess)=sessionmean(sess)./nr_subj;
                end

                subjmean=zeros(nr_subj,1);
                for sub=1:nr_subj
                    for sess=1:2
                        subjmean(sub)=subjmean(sub) + data(sub,sess);
                    end
                      subjmean(sub)=subjmean(sub)./2;
                end

                % mean squares
                BMS=0; % between subject
                WMS=0; % within subject 
                EMS=0; % error
                JMS=0; % session

                for sub=1:nr_subj,    
                    BMS = BMS + (subjmean(sub)-grandmean).^2;
                    for sess=1:2
                        WMS = WMS + (data(sub,sess)-subjmean(sub)).^2;
                        EMS = EMS + (data(sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
                    end
                end;

                for sess=1:2
                    JMS=  JMS + (sessionmean(sess)-grandmean).^2;
                end;

                %define the true value of the mean square.
                BMS= 2.*BMS./(nr_subj-1);
                WMS= WMS./(2-1)./nr_subj;
                JMS= nr_subj.*JMS./(2-1);
                EMS= EMS./(2-1)./(nr_subj-1); 

                %consistency agreement  
                if cons==1
                voxICC_con=(BMS-EMS)./(BMS+(2-1).*WMS); 
                ICC_con(ind_x, ind_y, ind_z) = voxICC_con;
                z_ICC_con(ind_x, ind_y, ind_z) = .5.*log((1+voxICC_con)./(1-voxICC_con));
                end;

                %absolute agreement 
                if abs==1
                voxICC_abs=(BMS-EMS)./(BMS+(2-1).*EMS + ...
                                                   2.* (JMS-EMS)./nr_subj);
                ICC_abs(ind_x, ind_y, ind_z) = voxICC_abs;
                z_ICC_abs(ind_x,ind_y,ind_z) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                end;
            end; 
        end;
    end;
    
disp('saving ICC images')
 % save ICC maps
target_img = temp_img;
target_img.fileprefix = sprintf('ICC_con_%d.nii',i_run);
target_img.img = ICC_con;
save_nii(target_img,target_img.fileprefix); 

clear target_img;
target_img = temp_img;
target_img.fileprefix = sprintf('ICC_abs_%d.nii',i_run);
target_img.img = ICC_abs;
save_nii(target_img,target_img.fileprefix); 

target_img = temp_img;
target_img.fileprefix = sprintf('z_ICC_con_%d.nii',i_run);
target_img.img = z_ICC_con;
save_nii(target_img,target_img.fileprefix); 

clear target_img;
target_img = temp_img;
target_img.fileprefix = sprintf('z_ICC_abs_%d.nii',i_run);
target_img.img = z_ICC_abs;
save_nii(target_img,target_img.fileprefix); 

%compute z mean ICC
mean_z = mean(z_ICC_con(:),'omitnan');
summary(1,end+1)=mean_z;
cols{1,end+1} = sprintf('mean_z_ICC_con_%d',i_run);
%compute min and max ICC
[my_min, idx] = min(z_ICC_con(:));
[my_max, idx] = max(z_ICC_con(:));
summary(1,end+1)= my_min;
cols{1,end+1} = sprintf('min_z_ICC_con_%d',i_run);
summary(1,end+1) = my_max;
cols{1,end+1} = sprintf('max_z_ICC_con_%d',i_run); 

%compute z mean ICC
mean_z = mean(z_ICC_abs(:),'omitnan');
summary(1,end+1)=mean_z;
cols{1,end+1} = sprintf('mean_z_ICC_abs_%d',i_run);
%compute min and max ICC
[my_min, idx] = min(z_ICC_abs(:));
[my_max, idx] = max(z_ICC_abs(:));
summary(1,end+1)= my_min;
cols{1,end+1} = sprintf('min_z_ICC_abs_%d',i_run);
summary(1,end+1) = my_max;
cols{1,end+1} = sprintf('max_z_ICC_abs_%d',i_run);     

end;

end;
    
assignin('base','summary',summary);
assignin('base','cols',cols);

results_ICCvox=dataset({summary(1,:),cols{:}});
assignin('base','results_ICCvox',results_ICCvox);
    
cd(dir_results);             
save results_ICCvox.mat results_ICCvox      

disp('finished ICC voxelwise')
cd(boxpath);
