function varargout = ICC_roi(varargin)
% ICC_ROI MATLAB code for ICC_roi.fig
%      ICC_ROI, by itself, creates a new ICC_ROI or raises the existing
%      singleton*.
%
%      H = ICC_ROI returns the handle to a new ICC_ROI or the handle to
%      the existing singleton*.
%
%      ICC_ROI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ICC_ROI.M with the given input arguments.
%
%      ICC_ROI('Property','Value',...) creates a new ICC_ROI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ICC_roi_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ICC_roi_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ICC_roi

% Last Modified by GUIDE v2.5 03-Mar-2017 11:04:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ICC_roi_OpeningFcn, ...
                   'gui_OutputFcn',  @ICC_roi_OutputFcn, ...
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


% --- Executes just before ICC_roi is made visible.
function ICC_roi_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ICC_roi (see VARARGIN)

% Choose default command line output for ICC_roi
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ICC_roi wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ICC_roi_OutputFcn(hObject, eventdata, handles) 
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
% --- Executes on button press in run.



% --- Executes on button press in roi_dir.
function roi_dir_Callback(hObject, eventdata, handles)
% hObject    handle to roi_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
roi_dir = spm_select(1,'dir','choose roi directory');
assignin('base','roi_dir',roi_dir);


% --- Executes on button press in roi.
function roi_Callback(hObject, eventdata, handles)
% hObject    handle to roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function roi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to roi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cons.
function cons_Callback(hObject, eventdata, handles)
% hObject    handle to cons (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cons


% --- Executes on button press in abs.
function abs_Callback(hObject, eventdata, handles)
% hObject    handle to abs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of abs


% --- Executes on button press in con_def.
function con_def_Callback(hObject, eventdata, handles)
% hObject    handle to con_def (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
con_def = cellstr(spm_select(1,'mat','load contrast definition'));
load(con_def{1});
assignin('base','contrast_def',contrast_def);


function run_Callback(hObject, eventdata, handles)
% hObject    handle to run (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('starting calculation of ICC for ROI');
%% set parameters
%get information out of study_design
box=pwd;
study_design=evalin('base','study_design');
contrast_def=evalin('base','contrast_def');

runs=str2double(study_design.number_sessions); 
nr_subj=str2double(study_design.number_subjects);
load(study_design.subject_list);
stats_dir=study_design.stats_directory;
stats_path=study_design.stats_path;

if runs == 1
    single_run = str2double(study_design.identifier_session);
end;

%get GUI input
roi=get(handles.roi,'String');
roi_dir=evalin('base','roi_dir');
split = get(handles.split,'value');
cons = get(handles.cons,'value');
abs = get(handles.abs,'value');

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



dir_results = study_design.results_directory;
%% ROI settings
% load ROI
cd(roi_dir);
roi_ful = dir(sprintf('%s*',roi));
if isstruct(roi_ful)
    if length(roi_ful)==2
        roi_ful = roi_ful(2).name;
    else
        roi_ful = roi_ful(1).name;
    end;
end;
compl = sprintf('%s\\%s',roi_dir, roi_ful);

%% reslice ROI

disp('...reslicing ROI...');
stats_filled = sprintf(stats_dir,1);
temp = sprintf('%s\\%s\\%s\\%s,1',stats_path,vp{1},stats_filled,con);
matlabbatch{1}.spm.spatial.coreg.write.ref = {temp};
temp_1 = sprintf('%s,1',compl);
matlabbatch{1}.spm.spatial.coreg.write.source = {temp_1};
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 4;
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';

% save batch
save('reslice_roi','matlabbatch');

% run batch
spm_jobman('serial',matlabbatch);

% create index for ROI voxels
r_roi=dir(sprintf('r%s*',roi));
if length(r_roi)==2
    compl1 = sprintf('%s\\r%s.img',roi_dir,roi);
    movefile(compl1,dir_results,'f');
    compl2 = sprintf('%s\\r%s.hdr',roi_dir,roi);
    movefile(compl2,dir_results,'f');
    cd(dir_results);    
    r_roi = load_nii(sprintf('r%s.img',roi));
    r_roi_ind = r_roi.img==1;
else
    compl = sprintf('%s\\r%s.nii',roi_dir,roi);
    movefile(compl,dir_results,'f');
    cd(dir_results);    
    r_roi = load_nii(sprintf('r%s.nii',roi));
    r_roi_ind = r_roi.img==1;        
end;

cd(dir_results);


%% load images and set voxels outside ROI = 0
disp('...loading images and mask with ROI...');
if split == 0 && two_cons == 0
    for i = 1:runs
        if runs == 1
            i = single_run;
        end;
        img = load_nii(sprintf('4D_%d.nii',i));
        evalstr = sprintf('img_%d = img;',i);
        eval(evalstr);
        evalstr = sprintf('img_%d = img_%d.img;',i,i);
        eval(evalstr);
        for count_subj = 1:nr_subj
            str = sprintf('img_temp = img_%d(:,:,:,count_subj);',i);
            eval(str);
            img_temp(~r_roi_ind) = 0;
            evalstr = sprintf('img_%d(:,:,:,count_subj) = img_temp;',i);
            eval(evalstr);
        end;
    end;
    disp('...loads image dimensions..');
    stats_temp =sprintf(stats_dir,1);
    temp_img = sprintf('%s\\%s\\%s\\%s',stats_path,vp{1},stats_temp,con);
    temp_img=load_nii(temp_img);
    dim = size(temp_img.img);
    x = dim(1);
    y = dim(2);
    z = dim(3);
elseif split == 1
    nr_para = study_design.number_parametric;

    for i = 1:runs 
        if runs == 1
            i = single_run;
        end;
        img1 = load_nii(sprintf('4D_split1_%d.nii',i));
        evalstr = sprintf('img_%d_split1 = img1;',i);
        eval(evalstr);
        evalstr = sprintf('img_%d_split1 = img_%d_split1.img;',i,i);
        eval(evalstr);
        for count_subj = 1:nr_subj
            str = sprintf('img_temp = img_%d_split1(:,:,:,count_subj);',i);
            eval(str);
            img_temp(~r_roi_ind) = 0;
            evalstr = sprintf('img_%d_split1(:,:,:,count_subj) = img_temp;',i);
            eval(evalstr);
        end;
        
        img2 = load_nii(sprintf('4D_split2_%d.nii',i));
        evalstr = sprintf('img_%d_split2 = img2;',i);
        eval(evalstr);
        evalstr = sprintf('img_%d_split2 = img_%d_split2.img;',i,i);
        eval(evalstr);
        for count_subj = 1:nr_subj
            str = sprintf('img_temp = img_%d_split2(:,:,:,count_subj);',i);
            eval(str);
            img_temp(~r_roi_ind) = 0;
            evalstr = sprintf('img_%d_split2(:,:,:,count_subj) = img_temp;',i);
            eval(evalstr);
        end;
    end;
    
    if nr_para > 0
        for ind_para = 1:nr_para
               for i = 1:runs    
                   if runs == 1
                       i = single_run;
                   end;
                    nii = sprintf('4D_split1_par%d_%d.nii',ind_para,i);
                    img1 = load_nii(nii);
                    evalstr = sprintf('img_%d_par%d_split1 = img1;',i,ind_para);
                    eval(evalstr);
                    evalstr = sprintf('img_%d_par%d_split1 = img_%d_par%d_split1.img;',i,ind_para,i,ind_para);
                    eval(evalstr);
                    for count_subj = 1:nr_subj
                        str = sprintf('img_temp = img_%d_par%d_split1(:,:,:,count_subj);',i,ind_para);
                        eval(str);
                        img_temp(~r_roi_ind) = 0;
                        evalstr = sprintf('img_%d_par%d_split1(:,:,:,count_subj) = img_temp;',i,ind_para);
                        eval(evalstr);
                    end;                   
                    
                    nii2 = sprintf('4D_split2_par%d_%d.nii',ind_para,i);
                    img2 = load_nii(nii2);
                    evalstr = sprintf('img_%d_par%d_split2 = img2;',i,ind_para);
                    eval(evalstr);
                    evalstr = sprintf('img_%d_par%d_split2 = img_%d_par%d_split2.img;',i,ind_para,i,ind_para);
                    eval(evalstr);
                    for count_subj = 1:nr_subj
                        str = sprintf('img_temp = img_%d_par%d_split2(:,:,:,count_subj);',i,ind_para);
                        eval(str);
                        img_temp(~r_roi_ind) = 0;
                        evalstr = sprintf('img_%d_par%d_split2(:,:,:,count_subj) = img_temp;',i,ind_para);
                        eval(evalstr);
                    end;   
                end; 
        end;
    end;
    disp('...loads image dimensions..');
    stats_temp =sprintf(stats,1);
    temp_img = sprintf('%s\\%s\\%s\\%s.nii',path,vp{1},stats_temp,con);
    temp_img=load_nii(temp_img);
    dim = size(temp_img.img);
    x = dim(1);
    y = dim(2);
    z = dim(3);
elseif two_cons == 1
    for i = 1:runs    
        if runs == 1
            i = single_run;
        end;
        img1 = load_nii(sprintf('4D_%s_%d.nii',con1,i));
        evalstr = sprintf('img_%d_con1 = img1;',i);
        eval(evalstr);
        evalstr = sprintf('img_%d_con1 = img_%d_con1.img;',i,i);
        eval(evalstr);
        for count_subj = 1:nr_subj
            str = sprintf('img_temp = img_%d_con1(:,:,:,count_subj);',i);
            eval(str);
            img_temp(~r_roi_ind) = 0;
            evalstr = sprintf('img_%d_con1(:,:,:,count_subj) = img_temp;',i);
            eval(evalstr);
        end;
        
        img2 = load_nii(sprintf('4D_%s_%d.nii',con2,i));
        evalstr = sprintf('img_%d_con2 = img2;',i);
        eval(evalstr);
        evalstr = sprintf('img_%d_con2 = img_%d_con2.img;',i,i);
        eval(evalstr);
        for count_subj = 1:nr_subj
            str = sprintf('img_temp = img_%d_con2(:,:,:,count_subj);',i);
            eval(str);
            img_temp(~r_roi_ind) = 0;
            evalstr = sprintf('img_%d_con2(:,:,:,count_subj) = img_temp;',i);
            eval(evalstr);
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



data=zeros(nr_subj,runs);
ICC_con_ROI=zeros(x,y,z);
ICC_abs_ROI=zeros(x,y,z);
z_ICC_con_ROI=zeros(x,y,z);
z_ICC_abs_ROI=zeros(x,y,z);

%%initiating summary
summary = [];
cols = {};
%% calculating ICCs
disp('calculating ICCs')
if split == 0 && two_cons == 0
    for ind_x = 1:x
        fprintf('..scanning voxel x = %d ...\n',ind_x)
        for ind_y = 1:y
           for ind_z = 1:z
                for ind_run = 1:runs
                    if runs == 1
                        ind_run = single_run;
                    end;
                    for ind_subj = 1:nr_subj
                       estr = sprintf('img%d_voxel(ind_subj,1)= img_%d (ind_x, ind_y, ind_z, ind_subj);',ind_run,ind_run);
                       eval(estr);
                    end;
                    estr = sprintf('data(:,ind_run)=img%d_voxel;',ind_run);
                    eval(estr);
                end;
                    %ICC
                    nsamples=nr_subj*runs;

                    grandmean=0;
                    for sub=1:nr_subj,     
                        for sess=1:runs,
                           grandmean= grandmean + data(sub,sess);
                        end

                    end;
                    grandmean=grandmean./nsamples;

                    sessionmean=zeros(runs,1);
                    for sess=1:runs
                        for sub=1:nr_subj,  
                            sessionmean(sess) = sessionmean(sess) + data(sub,sess);
                        end
                        sessionmean(sess)=sessionmean(sess)./nr_subj;
                    end

                    subjmean=zeros(nr_subj,1);
                    for sub=1:nr_subj
                        for sess=1:runs
                            subjmean(sub)=subjmean(sub) + data(sub,sess);
                        end
                          subjmean(sub)=subjmean(sub)./runs;
                    end

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
                        JMS=  JMS + (sessionmean(sess)-grandmean).^2;
                    end;

                    %define the true value of the mean square.
                    BMS= runs.*BMS./(nr_subj-1);
                    WMS= WMS./(runs-1)./nr_subj;
                    JMS= nr_subj.*JMS./(runs-1);
                    EMS= EMS./(runs-1)./(nr_subj-1); 

                    %consistency agreement  
                    if cons==1
                    voxICC_con=(BMS-EMS)./(BMS+(runs-1).*WMS); 
                    ICC_con_ROI(ind_x, ind_y, ind_z) = voxICC_con;
                    z_ICC_con_ROI(ind_x, ind_y, ind_z) = .5.*log((1+voxICC_con)./(1-voxICC_con));
                    end;
                    
                    %absolute agreement 
                    if abs==1
                    voxICC_abs=(BMS-EMS)./(BMS+(runs-1).*EMS + ...
                                                       runs.* (JMS-EMS)./nr_subj);
                    ICC_abs_ROI(ind_x, ind_y, ind_z) = voxICC_abs;
                    z_ICC_abs_ROI(ind_x,ind_y,ind_z) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                    end;
           end; 
        end;
     end;
disp('saving ICC images')
 % save ICC maps
target_img = temp_img;
target_img.fileprefix = 'ICC_con_ROI.nii';
target_img.img = ICC_con_ROI;
save_nii(target_img,target_img.fileprefix); 

clear target_img;
target_img = temp_img;
target_img.fileprefix = 'ICC_abs_ROI.nii';
target_img.img = ICC_abs_ROI;
save_nii(target_img,target_img.fileprefix); 

target_img = temp_img;
target_img.fileprefix = 'z_ICC_con_ROI.nii';
target_img.img = z_ICC_con_ROI;
save_nii(target_img,target_img.fileprefix); 

clear target_img;
target_img = temp_img;
target_img.fileprefix = 'z_ICC_abs_ROI.nii';
target_img.img = z_ICC_abs_ROI;
save_nii(target_img,target_img.fileprefix); 

% computing ROI ICC
disp('...computing ROI total ICC...')
ROI_subj = zeros(nr_subj,runs);
for ind_subj = 1:nr_subj
    for ind_run = 1:runs
        if runs == 1
            ind_run = single_run;
        end;
        eval(sprintf('img%d_%d = img_%d (:,:,:,ind_subj);',ind_run,ind_subj,ind_run));
        eval(sprintf('ROI_subj(ind_subj,ind_run) = mean(img%d_%d(r_roi_ind));',ind_run,ind_subj));
    end;
end;                
                    %ICC
                    nsamples=nr_subj*runs;

                    grandmean=0;
                    for sub=1:nr_subj,     
                        for sess=1:runs,
                           grandmean= grandmean + ROI_subj(sub,sess);
                        end

                    end;
                    grandmean=grandmean./nsamples;

                    sessionmean=zeros(runs,1);
                    for sess=1:runs
                        for sub=1:nr_subj,  
                            sessionmean(sess) = sessionmean(sess) + ROI_subj(sub,sess);
                        end
                        sessionmean(sess)=sessionmean(sess)./nr_subj;
                    end

                    subjmean=zeros(nr_subj,1);
                    for sub=1:nr_subj
                        for sess=1:runs
                            subjmean(sub)=subjmean(sub) + ROI_subj(sub,sess);
                        end
                          subjmean(sub)=subjmean(sub)./runs;
                    end

                    % mean squares
                    BMS=0; % between subject
                    WMS=0; % within subject 
                    EMS=0; % error
                    JMS=0; % session
                    
                    for sub=1:nr_subj,    
                        BMS = BMS + (subjmean(sub)-grandmean).^2;
                        for sess=1:runs
                            WMS = WMS + (ROI_subj(sub,sess)-subjmean(sub)).^2;
                            EMS = EMS + (ROI_subj(sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
                        end
                    end;

                    for sess=1:runs
                        JMS=  JMS + (sessionmean(sess)-grandmean).^2;
                    end;

                    %define the true value of the mean square.
                    BMS= runs.*BMS./(nr_subj-1);
                    WMS= WMS./(runs-1)./nr_subj;
                    JMS= nr_subj.*JMS./(runs-1);
                    EMS= EMS./(runs-1)./(nr_subj-1); 

                    %consistency agreement  
                    if cons==1
                    voxICC_con=(BMS-EMS)./(BMS+(runs-1).*WMS); 
                    ICC_con_ROI_total = voxICC_con;
                    z_ICC_con_ROI_total = .5.*log((1+voxICC_con)./(1-voxICC_con));
                    cols{1,end+1}= 'ICC_con_ROI';
                    summary(1,end+1)=ICC_con_ROI_total;
                    cols{1,end+1}='z_ICC_con_ROI';
                    summary(1,end+1)=z_ICC_con_ROI_total;
                    
                    end;

                    %absolute agreement 
                    if abs==1
                    voxICC_abs=(BMS-EMS)./(BMS+(runs-1).*EMS + ...
                                                       runs.* (JMS-EMS)./nr_subj);
                    ICC_abs_ROI_total = voxICC_abs;
                    z_ICC_abs_ROI_total = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                    cols{1,end+1}= 'ICC_abs_ROI';
                    summary(1,end+1)=ICC_abs_ROI_total;
                    cols{1,end+1}='z_ICC_abs_ROI';
                    summary(1,end+1)=z_ICC_abs_ROI_total;                    
                    end;

%% based on split-half
elseif split == 1
for ind = 1:runs
     for ind_x = 1:x
          fprintf('...x = %d...',ind_x)
        for ind_y = 1:y
           for ind_z = 1:z
                for ind_run = 1:runs
                    if runs == 1
                        ind_run = single_run;
                    end;
                    eval(sprintf('data_%d = zeros(nr_subj,2);',ind_run));
                    for ind_subj = 1:nr_subj
                       estr = sprintf('img%d_voxel_split1(ind_subj,1)= img_%d_split1 (ind_x, ind_y, ind_z, ind_subj);',ind_run,ind_run);
                       eval(estr);
                       estr2 = sprintf('img%d_voxel_split2(ind_subj,1)= img_%d_split2 (ind_x, ind_y, ind_z, ind_subj);',ind_run,ind_run);
                       eval(estr2);                       
                    end;
                    estr = sprintf('data_%d(:,1)=img%d_voxel_split1;',ind_run,ind_run);
                    eval(estr);
                    estr2 = sprintf('data_%d(:,2)=img%d_voxel_split2;',ind_run,ind_run);
                    eval(estr2);                    
                end;
                    %ICC
    
                        nsamples=nr_subj*runs;
                        
                        grandmean=0;
                        for sub=1:nr_subj,     
                            for sess=1:2,
                               eva=sprintf('grandmean= grandmean + data_%d(sub,sess);',ind);
                               eval(eva);
                            end
                        end;
                        grandmean=grandmean./nsamples;

                        sessionmean=zeros(2,1);
                        for sess=1:2
                            for sub=1:nr_subj,  
                                e = sprintf('sessionmean(sess) = sessionmean(sess) + data_%d(sub,sess);',ind);
                                eval(e);
                            end
                            sessionmean(sess)=sessionmean(sess)./nr_subj;
                        end

                        subjmean=zeros(nr_subj,1);
                        for sub=1:nr_subj
                            for sess=1:2
                                ev = sprintf('subjmean(sub)=subjmean(sub) + data_%d(sub,sess);',ind);
                                eval(ev);
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
                                evs = sprintf('WMS = WMS + (data_%d(sub,sess)-subjmean(sub)).^2;',ind);
                                eval(evs);
                                evst = sprintf('EMS = EMS + (data_%d(sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;',ind);
                                eval(evst);
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
                        
                        disp('...calculate voxel ICC...')
                        %consistency agreement  
                        if cons==1
                        voxICC_con=(BMS-EMS)./(BMS+(2-1).*WMS); 
                        ICC_con_ROI(ind_x, ind_y, ind_z) = voxICC_con;
                        z_ICC_con_ROI(ind_x, ind_y, ind_z) = .5.*log((1+voxICC_con)./(1-voxICC_con));
                        end;

                        %absolute agreement 
                        if abs==1
                        voxICC_abs=(BMS-EMS)./(BMS+(2-1).*EMS + ...
                                                           2.* (JMS-EMS)./nr_subj);
                        ICC_abs_ROI(ind_x, ind_y, ind_z) = voxICC_abs;
                        z_ICC_abs_ROI(ind_x,ind_y,ind_z) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                        end;
                                                
            end;
         end; 
      end;
                            disp('saving ICC images')
                             % save ICC maps
                            target_img = temp_img;
                            file1 = sprintf('ICC_con_%d_split_ROI.nii',ind);
                            target_img.fileprefix = file1;
                            target_img.img = ICC_con_ROI;
                            save_nii(target_img,target_img.fileprefix); 

                            clear target_img;
                            target_img = temp_img;
                            file2 = sprintf('ICC_abs_%d_split_ROI.nii',ind);
                            target_img.fileprefix = file2;
                            target_img.img = ICC_abs_ROI;
                            save_nii(target_img,target_img.fileprefix); 

                            target_img = temp_img;
                            file3 = sprintf('z_ICC_con_%d_split_ROI.nii',ind);
                            target_img.fileprefix = file3;
                            target_img.img = z_ICC_con_ROI;
                            save_nii(target_img,target_img.fileprefix); 

                            clear target_img;
                            target_img = temp_img;
                            file4 = sprintf('z_ICC_abs_%d_split_ROI.nii',ind);
                            target_img.fileprefix = file4;
                            target_img.img = z_ICC_abs_ROI;
                            save_nii(target_img,target_img.fileprefix);              
end;
   
     

     % computing ROI ICC
disp('...computing ROI total ICC...')

for ind_runs = 1:runs
    if runs == 1
        ind_runs = single_run;
    end;
    ROI_subj = zeros(nr_subj,2);
    for ind_subj = 1:nr_subj
        for ind_split = 1:2
            eval(sprintf('img%d_%d = img_%d_split%d (:,:,:,ind_subj);',ind_split,ind_subj,ind_run,ind_split));
            eval(sprintf('ROI_subj(ind_subj,ind_split) = mean(img%d_%d(r_roi_ind));',ind_split,ind_subj));
        end;
    end;                
                    %ICC
                    nsamples=nr_subj*2;

                    grandmean=0;
                    for sub=1:nr_subj,     
                        for sess=1:2,
                           grandmean= grandmean + ROI_subj(sub,sess);
                        end

                    end;
                    grandmean=grandmean./nsamples;

                    sessionmean=zeros(2,1);
                    for sess=1:2
                        for sub=1:nr_subj,  
                            sessionmean(sess) = sessionmean(sess) + ROI_subj(sub,sess);
                        end
                        sessionmean(sess)=sessionmean(sess)./nr_subj;
                    end

                    subjmean=zeros(nr_subj,1);
                    for sub=1:nr_subj
                        for sess=1:2
                            subjmean(sub)=subjmean(sub) + ROI_subj(sub,sess);
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
                            WMS = WMS + (ROI_subj(sub,sess)-subjmean(sub)).^2;
                            EMS = EMS + (ROI_subj(sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
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
                    ICC_con_ROI_total = voxICC_con;
                    z_ICC_con_ROI_total = .5.*log((1+voxICC_con)./(1-voxICC_con));
                    cols{1,end+1}= sprintf('ICC_con_ROI_split%d',ind_runs);
                    summary(1,end+1)=ICC_con_ROI_total;
                    cols{1,end+1}=sprintf('z_ICC_con_ROI_split%d',ind_runs);
                    summary(1,end+1)=z_ICC_con_ROI_total;
                    end;
                    

                    
                    %absolute agreement 
                    if abs==1
                    voxICC_abs=(BMS-EMS)./(BMS+(2-1).*EMS + ...
                                                       2.* (JMS-EMS)./nr_subj);
                    ICC_abs_ROI_total = voxICC_abs;
                    z_ICC_abs_ROI_total = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                    cols{1,end+1}= sprintf('ICC_abs_ROI_split%d',ind_runs);
                    summary(1,end+1)=ICC_abs_ROI_total;
                    cols{1,end+1}=sprintf('z_ICC_abs_ROI_split%d',ind_runs);
                    summary(1,end+1)=z_ICC_abs_ROI_total;
                    end;


end;
if nr_para > 0
    fprintf('...ICC parametric modulator...\n')
        for ind_para = 1:nr_para
             for ind_x = 1:x
                 fprintf('...scanning x = %d ...',ind_x)
                for ind_y = 1:y
                   for ind_z = 1:z
                        for ind_run = 1:runs
                            if runs == 1
                                ind_run = single_run;
                            end;
                            eval(sprintf('data_%d = zeros(nr_subj,2);',ind_run));
                            for ind_subj = 1:nr_subj
                               eval(sprintf('img%d_par%d_voxel_split1(ind_subj,1)= img_%d_par%d_split1 (ind_x, ind_y, ind_z, ind_subj);',ind_run,ind_para,ind_run,ind_para));
                               eval(sprintf('img%d_par%d_voxel_split2(ind_subj,1)= img_%d_par%d_split2 (ind_x, ind_y, ind_z, ind_subj);',ind_run,ind_para,ind_run,ind_para));
                            end;
                            eval(sprintf('data_%d(:,1)=img%d_par%d_voxel_split1;',ind_run,ind_run,ind_para));
                            eval(sprintf('data_%d(:,2)=img%d_par%d_voxel_split2;',ind_run,ind_run,ind_para));
                        end;
                    end;
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
                                       eva=sprintf('grandmean= grandmean + data_%d(sub,sess);',ind_run);
                                       eval(eva);
                                    end
                                end;
                                grandmean=grandmean./nsamples;

                                sessionmean=zeros(2,1);
                                for sess=1:2
                                    for sub=1:nr_subj,  
                                        e = sprintf('sessionmean(sess) = sessionmean(sess) + data_%d(sub,sess);',ind_run);
                                        eval(e);
                                    end
                                    sessionmean(sess)=sessionmean(sess)./nr_subj;
                                end

                                subjmean=zeros(nr_subj,1);
                                for sub=1:nr_subj
                                    for sess=1:2
                                        ev = sprintf('subjmean(sub)=subjmean(sub) + data_%d(sub,sess);',ind_run);
                                        eval(ev);
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
                                        evs = sprintf('WMS = WMS + (data_%d(sub,sess)-subjmean(sub)).^2;',ind_run);
                                        eval(evs);
                                        evst = sprintf('EMS = EMS + (data_%d(sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;',ind_run);
                                        eval(evst);
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
                                ICC_con_ROI(ind_x, ind_y, ind_z) = voxICC_con;
                                z_ICC_con_ROI(ind_x, ind_y, ind_z) = .5.*log((1+voxICC_con)./(1-voxICC_con));
                                end;

                                %absolute agreement 
                                if abs==1
                                voxICC_abs=(BMS-EMS)./(BMS+(2-1).*EMS + ...
                                                                   2.* (JMS-EMS)./nr_subj);
                                ICC_abs_ROI(ind_x, ind_y, ind_z) = voxICC_abs;
                                z_ICC_abs_ROI(ind_x,ind_y,ind_z) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                                end;

                                             disp('saving ICC images')
                                     % save ICC maps
                                    target_img = temp_img;
                                    file1 = sprintf('ICC_con_%d_split_par%d_ROI.nii',ind_run,ind_para);
                                    target_img.fileprefix = file1;
                                    target_img.img = ICC_con_ROI;
                                    save_nii(target_img,target_img.fileprefix); 

                                    clear target_img;
                                    target_img = temp_img;
                                    file2 = sprintf('ICC_abs_%d_split_par%d_ROI.nii',ind_run,ind_para);
                                    target_img.fileprefix = file2;
                                    target_img.img = ICC_abs_ROI;
                                    save_nii(target_img,target_img.fileprefix); 

                                    target_img = temp_img;
                                    file3 = sprintf('z_ICC_con_%d_split_par%d_ROI.nii',ind_run,ind_para);
                                    target_img.fileprefix = file3;
                                    target_img.img = z_ICC_con_ROI;
                                    save_nii(target_img,target_img.fileprefix); 

                                    clear target_img;
                                    target_img = temp_img;
                                    file4 = sprintf('z_ICC_abs_%d_split_par%d_ROI.nii',ind_run,ind_para);
                                    target_img.fileprefix = file4;
                                    target_img.img = z_ICC_abs_ROI;
                                    save_nii(target_img,target_img.fileprefix); 
        end;
        end;   
        
for ind_para = 1:nr_para

for ind_runs = 1:runs
    if runs == 1
        ind_runs = single_run;
    end;
    ROI_subj = zeros(nr_subj,2);
    for ind_subj = 1:nr_subj
        for ind_split = 1:2
            eval(sprintf('img%d_%d = img_%d_par%d_split%d (:,:,:,ind_subj);',ind_split,ind_subj,ind_runs,ind_para,ind_split));
            eval(sprintf('ROI_subj(ind_subj,ind_split) = mean(img%d_%d(r_roi_ind));',ind_split,ind_subj));
        end;
    end;                
                    %ICC
                    nsamples=nr_subj*2;

                    grandmean=0;
                    for sub=1:nr_subj,     
                        for sess=1:2,
                           grandmean= grandmean + ROI_subj(sub,sess);
                        end

                    end;
                    grandmean=grandmean./nsamples;

                    sessionmean=zeros(2,1);
                    for sess=1:2
                        for sub=1:nr_subj,  
                            sessionmean(sess) = sessionmean(sess) + ROI_subj(sub,sess);
                        end
                        sessionmean(sess)=sessionmean(sess)./nr_subj;
                    end

                    subjmean=zeros(nr_subj,1);
                    for sub=1:nr_subj
                        for sess=1:2
                            subjmean(sub)=subjmean(sub) + ROI_subj(sub,sess);
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
                            WMS = WMS + (ROI_subj(sub,sess)-subjmean(sub)).^2;
                            EMS = EMS + (ROI_subj(sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
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
                    ICC_con_ROI_total = voxICC_con;
                    z_ICC_con_ROI_total = .5.*log((1+voxICC_con)./(1-voxICC_con));
                    
                    cols{1,end+1}= sprintf('ICC_con_ROI_split%d_par%d',ind_runs,ind_para);
                    summary(1,end+1)=ICC_con_ROI_total;
                    cols{1,end+1}=sprintf('z_ICC_con_ROI_split%d_par%d',ind_runs,ind_para);
                    summary(1,end+1)=z_ICC_con_ROI_total;                    
                    end;
                    
                    %absolute agreement 
                    if abs==1
                    voxICC_abs=(BMS-EMS)./(BMS+(2-1).*EMS + ...
                                                       2.* (JMS-EMS)./nr_subj);
                    ICC_abs_ROI_total = voxICC_abs;
                    z_ICC_abs_ROI_total = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                    
                    cols{1,end+1}= sprintf('ICC_abs_ROI_split%d_par%d',ind_runs,ind_para);
                    summary(1,end+1)=ICC_abs_ROI_total;
                    cols{1,end+1}=sprintf('z_ICC_abs_ROI_split%d_par%d',ind_runs,ind_para);
                    summary(1,end+1)=z_ICC_abs_ROI_total;                     
                    end;
                    
                  
                    
                    
end;
end;
end;

%% two contrasts out of one statistic
elseif two_cons == 1
for ind = 1:runs
    if runs == 1
        ind = single_run;
    end;
     eval(sprintf('data_%d = zeros(nr_subj,2);',ind));
    for ind_x = 1:x
          fprintf('...x = %d...\n',ind_x)
        for ind_y = 1:y
            for ind_z = 1:z
                    for ind_subj = 1:nr_subj
                       eval(sprintf('img%d_voxel_con1(ind_subj,1)= img_%d_con1 (ind_x, ind_y, ind_z, ind_subj);',ind,ind));
                       eval(sprintf('img%d_voxel_con2(ind_subj,1)= img_%d_con2 (ind_x, ind_y, ind_z, ind_subj);',ind,ind));
                    end;
                    eval(sprintf('data_%d(:,1)=img%d_voxel_con1;',ind,ind));
                    eval(sprintf('data_%d(:,2)=img%d_voxel_con2;',ind,ind));
                    %ICC
    
                        nsamples=nr_subj*runs;
                        
                        grandmean=0;
                        for sub=1:nr_subj,     
                            for sess=1:2,
                               eva=sprintf('grandmean= grandmean + data_%d(sub,sess);',ind);
                               eval(eva);
                            end
                        end;
                        grandmean=grandmean./nsamples;

                        sessionmean=zeros(2,1);
                        for sess=1:2
                            for sub=1:nr_subj,  
                                e = sprintf('sessionmean(sess) = sessionmean(sess) + data_%d(sub,sess);',ind);
                                eval(e);
                            end
                            sessionmean(sess)=sessionmean(sess)./nr_subj;
                        end

                        subjmean=zeros(nr_subj,1);
                        for sub=1:nr_subj
                            for sess=1:2
                                ev = sprintf('subjmean(sub)=subjmean(sub) + data_%d(sub,sess);',ind);
                                eval(ev);
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
                                evs = sprintf('WMS = WMS + (data_%d(sub,sess)-subjmean(sub)).^2;',ind);
                                eval(evs);
                                evst = sprintf('EMS = EMS + (data_%d(sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;',ind);
                                eval(evst);
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
                        ICC_con_ROI(ind_x, ind_y, ind_z) = voxICC_con;
                        z_ICC_con_ROI(ind_x, ind_y, ind_z) = .5.*log((1+voxICC_con)./(1-voxICC_con));
                        end;

                        %absolute agreement 
                        if abs==1
                        voxICC_abs=(BMS-EMS)./(BMS+(2-1).*EMS + ...
                                                           2.* (JMS-EMS)./nr_subj);
                        ICC_abs_ROI(ind_x, ind_y, ind_z) = voxICC_abs;
                        z_ICC_abs_ROI(ind_x,ind_y,ind_z) = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                        end;
                                                
            end;
        end; 
    end;
                            disp('saving ICC images')
                             % save ICC maps
                            target_img = temp_img;
                            file1 = sprintf('ICC_con_%d_con_ROI.nii',ind);
                            target_img.fileprefix = file1;
                            target_img.img = ICC_con_ROI;
                            save_nii(target_img,target_img.fileprefix); 

                            clear target_img;
                            target_img = temp_img;
                            file2 = sprintf('ICC_abs_%d_con_ROI.nii',ind);
                            target_img.fileprefix = file2;
                            target_img.img = ICC_abs_ROI;
                            save_nii(target_img,target_img.fileprefix); 

                            target_img = temp_img;
                            file3 = sprintf('z_ICC_con_%d_con_ROI.nii',ind);
                            target_img.fileprefix = file3;
                            target_img.img = z_ICC_con_ROI;
                            save_nii(target_img,target_img.fileprefix); 

                            clear target_img;
                            target_img = temp_img;
                            file4 = sprintf('z_ICC_abs_%d_con_ROI.nii',ind);
                            target_img.fileprefix = file4;
                            target_img.img = z_ICC_abs_ROI;
                            save_nii(target_img,target_img.fileprefix);              
end;
   
     

     % computing ROI ICC
disp('...computing ROI total ICC...')


for ind_runs = 1:runs
    if runs == 1
        ind_runs = single_run;
    end;
    ROI_subj = zeros(nr_subj,2);
  for ind_con = 1:2
    for ind_subj = 1:nr_subj
            eval(sprintf('img%d_%d = img_%d_con%d (:,:,:,ind_subj);',ind_con,ind_subj,ind_runs,ind_con));
            eval(sprintf('ROI_subj(ind_subj,ind_con) = mean(img%d_%d(r_roi_ind));',ind_con,ind_subj));
    end;
  end;
                    %ICC
                    nsamples=nr_subj*2;

                    grandmean=0;
                    for sub=1:nr_subj,     
                        for sess=1:2,
                           grandmean= grandmean + ROI_subj(sub,sess);
                        end

                    end;
                    grandmean=grandmean./nsamples;

                    sessionmean=zeros(2,1);
                    for sess=1:2
                        for sub=1:nr_subj,  
                            sessionmean(sess) = sessionmean(sess) + ROI_subj(sub,sess);
                        end
                        sessionmean(sess)=sessionmean(sess)./nr_subj;
                    end

                    subjmean=zeros(nr_subj,1);
                    for sub=1:nr_subj
                        for sess=1:2
                            subjmean(sub)=subjmean(sub) + ROI_subj(sub,sess);
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
                            WMS = WMS + (ROI_subj(sub,sess)-subjmean(sub)).^2;
                            EMS = EMS + (ROI_subj(sub,sess)-subjmean(sub)-sessionmean(sess)+grandmean).^2;
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
                    ICC_con_ROI_total = voxICC_con;
                    z_ICC_con_ROI_total = .5.*log((1+voxICC_con)./(1-voxICC_con));
                    cols{1,end+1}= sprintf('ICC_con_ROI_%d',ind_runs);
                    summary(1,end+1)=ICC_con_ROI_total;
                    cols{1,end+1}=sprintf('z_ICC_con_ROI_%d',ind_runs);
                    summary(1,end+1)=z_ICC_con_ROI_total;
                    end;
                    

                    
                    %absolute agreement 
                    if abs==1
                    voxICC_abs=(BMS-EMS)./(BMS+(2-1).*EMS + ...
                                                       2.* (JMS-EMS)./nr_subj);
                    ICC_abs_ROI_total = voxICC_abs;
                    z_ICC_abs_ROI_total = .5.*log((1+voxICC_abs)./(1-voxICC_abs));
                    cols{1,end+1}= sprintf('ICC_abs_ROI_%d',ind_runs);
                    summary(1,end+1)=ICC_abs_ROI_total;
                    cols{1,end+1}=sprintf('z_ICC_abs_ROI_%d',ind_runs);
                    summary(1,end+1)=z_ICC_abs_ROI_total;
                    end;


end;
end;

assignin('base','summary',summary);
assignin('base','cols',cols);

results_ICCroi=dataset({summary(1,:),cols{:}});
assignin('base','results_ICCroi',results_ICCroi);
    
cd(dir_results);             
save results_ICCroi.mat results_ICCroi    
       

disp('finished ICC ROI')

cd(box);
