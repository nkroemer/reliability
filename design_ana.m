function varargout = design_ana(varargin)
% DESIGN_ANA MATLAB code for design_ana.fig
%      DESIGN_ANA, by itself, creates a new DESIGN_ANA or raises the existing
%      singleton*.
%
%      H = DESIGN_ANA returns the handle to a new DESIGN_ANA or the handle to
%      the existing singleton*.
%
%      DESIGN_ANA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DESIGN_ANA.M with the given input arguments.
%
%      DESIGN_ANA('Property','Value',...) creates a new DESIGN_ANA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before design_ana_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to design_ana_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help design_ana

% Last Modified by GUIDE v2.5 13-Apr-2018 09:47:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @design_ana_OpeningFcn, ...
                   'gui_OutputFcn',  @design_ana_OutputFcn, ...
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


% --- Executes just before design_ana is made visible.
function design_ana_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to design_ana (see VARARGIN)

% Choose default command line output for design_ana
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes design_ana wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = design_ana_OutputFcn(hObject, eventdata, handles) 
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

function iter_Callback(hObject, eventdata, handles)
% hObject    handle to iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iter as text
%        str2double(get(hObject,'String')) returns contents of iter as a double

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
con_def = cellstr(spm_select(1,'mat','load contrast definition'));
load(con_def{1});
assignin('base','contrast_def',contrast_def);

% --- Executes during object creation, after setting all properties.
function iter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iter (see GCBO)
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
f = filesep;
box_path=pwd;
disp('Starting design analysis...');
study_design=evalin('base','study_design');
contrast_def = evalin('base','contrast_def');

%% get study design info
results_dir = study_design.results_directory;
runs = str2double(study_design.number_sessions);
nr_subj = str2double(study_design.number_subjects);
load(study_design.subject_list);
stats=study_design.stats_directory;
path=study_design.stats_path;
split_dir = study_design.split_directory;
nr_split = study_design.nr_split_reg;
nr_para = study_design.number_parametric;

%% get contrast info
two_cons = contrast_def.two_contrasts;
if two_cons == 1
    con1=contrast_def.contrast1;
    con2=contrast_def.contrast2;
end;


%% get GUI input

sim.it = str2double(get(handles.iter,'String')); %number of iterations for the simulation

%% adapted from: fmreli simulation, 2018
%contact: nils.kroemer@uni-tuebingen.de
%v0.2 Mar 13 2018
%for details, see preprint: https://doi.org/10.1101/215053 

%% define VARs

%sim.n = nr_subj;            %number of participants
sim.n = 126;

sim.contrast = nr_split+nr_para;       %number of contrasts to do simulations for
%load exemplary SPM.mat whole design
cd(path);
cd(id{1});             
cd(sprintf(split_dir,1));
load('SPM.mat');
sim.TR = SPM.nscan;           %number of TRs
        
%% calculate average amplitude and SD
cd(results_dir);
disp('...calculate average amplitude and SD...')
if two_cons==1
    %con1 
    Four_D = load_untouch_nii(sprintf('4D_%s_1.nii',con1));
    beta_subj = zeros(nr_subj,1);
    if ~isa(Four_D.img,'double')
        Four_D.img = double(Four_D.img);
    end
    for i_subj = 1:nr_subj
        subj=Four_D.img(:,:,:,i_subj);
        beta_subj(i_subj,1)=mean(subj(:),'omitnan');
    end;
    beta_avg1 = mean(beta_subj);
    beta_std1 = std(beta_subj);
    
    Four_D = load_untouch_nii(sprintf('4D_%s_1.nii',con2));
    beta_subj = zeros(nr_subj,1);
    %con2
    if ~isa(Four_D.img,'double')
        Four_D.img = double(Four_D.img);
    end    
    for i_subj = 1:nr_subj
        subj=Four_D.img(:,:,:,i_subj);
        beta_subj(i_subj,1)=mean(subj(:),'omitnan');
    end;
    beta_avg2 = mean(beta_subj);
    beta_std2 = std(beta_subj);
    
else
    
    Four_D = load_untouch_nii('4D_1.nii');
    beta_subj = zeros(nr_subj,1);
    if ~isa(Four_D.img,'double')
        Four_D.img = double(Four_D.img);
    end
    for i_subj = 1:nr_subj
        subj=Four_D.img(:,:,:,i_subj);
        beta_subj(i_subj,1)=mean(subj(:),'omitnan');
    end;
    beta_avg = mean(beta_subj);
    %beta_avg=1;
    beta_std = std(beta_subj);
    %beta_std=1;
end;

if nr_para > 0
    for i_par = 1:nr_para
    Four_D = load_untouch_nii(sprintf('4D_par%d_1.nii',i_par));
    beta_subj = zeros(nr_subj,1);
    if ~isa(Four_D.img,'double')
        Four_D.img = double(Four_D.img);
    end
    for i_subj = 1:nr_subj
        subj=Four_D.img(:,:,:,i_subj);
        beta_subj(i_subj,1)=mean(subj(:),'omitnan');
    end;
    eval(sprintf('beta_avg_par%d = mean(beta_subj);',i_par));
    eval(sprintf('beta_std_par%d = std(beta_subj);',i_par));
    end;
end;

%% do all simulations for each contrast and each parametric
for i_con = 1:nr_split
    if i_con == 1
        for i_reg = 1:length(SPM.xCon(1).c)
            if SPM.xCon(1).c(i_reg) == 1
                split(1)=i_reg;
            end;
            if SPM.xCon(2+nr_para).c(i_reg) == 1
                split(2)=i_reg;
            end;
        end;
    else
        for i_reg = 1:length(SPM.xCon(1).c)
            if SPM.xCon(((i_con-1)*2)+((i_con-1)*2*nr_para)).c(i_reg) == 1
                split(1)=i_reg;
            end;
            if SPM.xCon(((i_con-1)*2)+((i_con-1)*2*nr_para)+ 2 + nr_para).c(i_reg) == 1
                split(2)=i_reg;
            end;      
        end;
    end;

    fprintf('...do simulation for contrast %d...\n',i_con)
    if two_cons == 0
        sim.beta1 =   beta_avg;     %parameter for amplitude, half 1
        sim.beta2 =   beta_avg;     %parameter for amplitude, half 2
        sim.betaSD1 = beta_std;     %parameter for SD, half 1
        sim.betaSD2 = beta_std;     %parameter for SD, half 2
    else
        eval(sprintf('sim.beta1 =   beta_avg%d;',i_con));     %parameter for amplitude, half 1
        eval(sprintf('sim.beta2 =   beta_avg%d;',i_con));     %parameter for amplitude, half 2
        eval(sprintf('sim.betaSD1 = beta_std%d;',i_con));     %parameter for SD, half 1
        eval(sprintf('sim.betaSD2 = beta_std%d;',i_con));     %parameter for SD, half 2
    end;


%initialize VARs for performance
b1 = zeros(sim.n,1);
b2 = zeros(sim.n,1);
true_r = zeros(sim.it,1);
sim_r = zeros(sim.it,1);
del_r = zeros(sim.it,1);

%define range of beta coefficients to introduce dependence of splits
%r_beta_vec = -1:0.05:1;
r_beta_vec = 0;

%get filtered and prewhitened design matrix
DesMat = SPM.xX.xKXs.X;

%check correlation between split regressors
[res.r, res.p] = corrcoef(DesMat(:,split));

%loop through different levels of statistical dependence
for i_beta = 1:length(r_beta_vec)
    %fprintf('beta %d \n',r_beta_vec(i_beta))

    %iterations of the simulation
    for i_it = 1:sim.it
    
        %randomly sample from a gaussian for beta 1 and predict beta2 based
        %on the level of statistical dependence
        beta1 = normrnd(sim.beta1,sim.betaSD1,sim.n,1);
        beta2 = r_beta_vec(i_beta) * beta1 + normrnd(sim.beta2,sim.betaSD2,sim.n,1);

        %get correlation of the true betas
        true_r(i_it,1) = corr(beta1,beta2);

        %simulated IDs
        for i_ID = 1:sim.n

            %use of different design matrices for participants within a
            %sample, this drastically slows down the execution of the
            %simulation
%             i_dice = randi(4,1);
%             switch i_dice
%                 case 1
%                 load('SPM_split1.mat');
%                 case 2
%                 load('SPM_split2.mat');
%                 case 3
%                 load('SPM_split3.mat');
%                 case 4
%                 load('SPM_split4.mat');
%             end
%             DesMat = SPM.xX.xKXs.X;
            
            %build one "observed" timecourse based on the two beta and the
            %corresponding regressors in the design matrix and calculate
            %the betas
                timec1 = beta1(i_ID,1) .* DesMat(:,split(1)) + beta2(i_ID,1) .* DesMat(:,split(2)) + normrnd(0,1,sim.TR,1); 
                simbs1 = regress(timec1,DesMat(:,:));
                b1(i_ID,1) = simbs1(split(1));
                b2(i_ID,1) = simbs1(split(2));
                 
        end

        %calculate the correlation between the split betas from the
        %simulated timeseries and the difference between the "true" and the
        %obeserved correlation
        sim_r(i_it,1) = corr(b1,b2);
        del_r(i_it,1) = atanh(sim_r(i_it,1)) - atanh(true_r(i_it,1));

    end
    
    %plots for inspection within levels of statistical dependence
    % figure
    % histogram(del_r,'normalization','pdf')
    % vline(median(del_r),'r','median');
    % 
    % figure
    % histogram(true_r,'normalization','pdf')
    % hold on
    % histogram(sim_r,'normalization','pdf')
    % hold off
    % 
    % figure;
    % scatter(true_r,sim_r)
    % 
    % figure;
    % plot(DesMat(1:500,1))
    % hold on
    % plot(DesMat(1:500,3))

%store percentiles
CI_del_r(i_beta,1) = prctile(del_r,2.5);
Med_del_r(i_beta,1) = prctile(del_r,50);
CI_del_r(i_beta,2) = prctile(del_r,97.5);

Med_true_r(i_beta,1) = median(true_r);
CI_true_r(i_beta,1) = prctile(true_r,2.5);
CI_true_r(i_beta,2) = prctile(true_r,97.5);

Med_sim_r(i_beta,1) = median(sim_r);
end
sum_des_ana = [Med_del_r,CI_del_r];

%%plot and save results
cd(results_dir);

%plot bias curves
% figure;
% plot(Med_true_r,Med_del_r);
% hold on;
% plot(CI_true_r,CI_del_r);
% hline(0);
% vline(0);
% axis([-0.7 0.7 -0.7 0.7]);
% fig = gcf;
% fig.PaperPositionMode = 'auto';
 name=sprintf('median_CI_con%d',i_con);
% print(fig,name,'-dpng','-r0');
% hold off;


save([name '.mat'],'sum_des_ana','res','beta_avg');


% figure
% plot(Med_true_r,Med_sim_r);
% hold on
% hline(0);
% hline(abs(corr(DesMat(:,2),DesMat(:,4))));
% vline(0);
% axis([-0.7 0.7 -0.7 0.7])
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% name=sprintf('correlation_con%d',i_con);
% print(fig,name,'-dpng','-r0')
% hold off

%make plot for different noise levels
% load('sim_noise025.mat');
% figure
% plot(Med_true_r,Med_sim_r);
% hold on
% load('sim_noise1.mat');
% plot(Med_true_r,Med_sim_r);
% load('sim_noise4.mat');
% plot(Med_true_r,Med_sim_r);
% hline(0);
% hline(abs(corr(DesMat(:,1),DesMat(:,3))));
% vline(0);
% axis([-0.7 0.7 -0.7 0.7])
% hold off
if nr_para>0
    for i_par = 1:nr_para
                    fprintf('...do simulation for parametric contrast %d...\n',i_par)

        for i_reg = 1:length(SPM.xCon(1).c)
        if i_con == 1
            if SPM.xCon(1+i_par).c(i_reg) == 1
                split(1)=i_reg;
            end;
            if SPM.xCon(2+nr_para+i_par).c(i_reg) == 1
                split(2)=i_reg;
            end;
        else
            if SPM.xCon(2*(i_con+nr_para)+1+i_par).c(i_reg) == 1
                split(1)=i_reg;
            end;
            if SPM.xCon(2*(i_con+nr_para)+2+nr_para).c(i_reg) == 1
                split(2)=i_reg;
            end;            
        end;
        end;

    eval(sprintf('sim.beta1 =   beta_avg_par%d;',i_par));     %parameter for amplitude, half 1
    eval(sprintf('sim.beta2 =   beta_avg_par%d;',i_par));     %parameter for amplitude, half 2
    eval(sprintf('sim.betaSD1 = beta_std_par%d;',i_par));     %parameter for SD, half 1
    eval(sprintf('sim.betaSD2 = beta_std_par%d;',i_par));     %parameter for SD, half 2


    %initialize VARs for performance
    b1 = zeros(sim.n,1);
    b2 = zeros(sim.n,1);
    true_r = zeros(sim.it,1);
    sim_r = zeros(sim.it,1);
    del_r = zeros(sim.it,1);

    %define range of beta coefficients to introduce dependence of splits
    %r_beta_vec = -1:0.05:1;
    r_beta_vec = 0;
    %get filtered and prewhitened design matrix
    DesMat = SPM.xX.xKXs.X;

    %check correlation between split regressors
    [res.r, res.p] = corrcoef(DesMat(:,split));

    %loop through different levels of statistical dependence
    for i_beta = 1:length(r_beta_vec)

        %iterations of the simulation
        for i_it = 1:sim.it

            %randomly sample from a gaussian for beta 1 and predict beta2 based
            %on the level of statistical dependence
            beta1 = normrnd(sim.beta1,sim.betaSD1,sim.n,1);
            beta2 = r_beta_vec(i_beta) * beta1 + normrnd(sim.beta2,sim.betaSD2,sim.n,1);

            %get correlation of the true betas
            true_r(i_it,1) = corr(beta1,beta2);

            %simulated IDs
            for i_ID = 1:sim.n

                %use of different design matrices for participants within a
                %sample, this drastically slows down the execution of the
                %simulation
    %             i_dice = randi(4,1);
    %             switch i_dice
    %                 case 1
    %                 load('SPM_split1.mat');
    %                 case 2
    %                 load('SPM_split2.mat');
    %                 case 3
    %                 load('SPM_split3.mat');
    %                 case 4
    %                 load('SPM_split4.mat');
    %             end
    %             DesMat = SPM.xX.xKXs.X;

                %build one "observed" timecourse based on the two beta and the
                %corresponding regressors in the design matrix and calculate
                %the betas

                    timec1 = beta1(i_ID,1) .* DesMat(:,split(1)) + beta2(i_ID,1) .* DesMat(:,split(2)) + normrnd(0,1,sim.TR,1); 
                    simbs1 = regress(timec1,DesMat(:,:));
                    b1(i_ID,1) = simbs1(split(1));
                    b2(i_ID,1) = simbs1(split(2));


            end

            %calculate the correlation between the split betas from the
            %simulated timeseries and the difference between the "true" and the
            %obeserved correlation
            sim_r(i_it,1) = corr(b1,b2);
            del_r(i_it,1) = atanh(sim_r(i_it,1)) - atanh(true_r(i_it,1));

        end

        %plots for inspection within levels of statistical dependence
        % figure
        % histogram(del_r,'normalization','pdf')
        % vline(median(del_r),'r','median');
        % 
        % figure
        % histogram(true_r,'normalization','pdf')
        % hold on
        % histogram(sim_r,'normalization','pdf')
        % hold off
        % 
        % figure;
        % scatter(true_r,sim_r)
        % 
        % figure;
        % plot(DesMat(1:500,1))
        % hold on
        % plot(DesMat(1:500,3))

    %store percentiles
    CI_del_r(i_beta,1) = prctile(del_r,2.5);
    Med_del_r(i_beta,1) = prctile(del_r,50);
    CI_del_r(i_beta,2) = prctile(del_r,97.5);

    Med_true_r(i_beta,1) = median(true_r);
    CI_true_r(i_beta,1) = prctile(true_r,2.5);
    CI_true_r(i_beta,2) = prctile(true_r,97.5);
    Med_sim_r(i_beta,1) = median(sim_r);
    end

    %%plot and save results
    cd(results_dir);

    %plot bias curves
    % figure;
    % plot(Med_true_r,Med_del_r);
    % hold on;
    % plot(CI_true_r,CI_del_r);
    % hline(0);
    % vline(0);
    % axis([-0.7 0.7 -0.7 0.7]);
    % fig = gcf;
    % fig.PaperPositionMode = 'auto';
     name=sprintf('median_CI_par%d',i_par);
    % print(fig,name,'-dpng','-r0');
    % hold off;

    sum_des_ana = [Med_del_r,CI_del_r];
    name1=sprintf('beta_avg_par%d',i_par);
    save([name '.mat'],'sum_des_ana','res',name1);

    % figure
    % plot(Med_true_r,Med_sim_r);
    % hold on
    % hline(0);
    % hline(abs(corr(DesMat(:,2),DesMat(:,4))));
    % vline(0);
    % axis([-0.7 0.7 -0.7 0.7])
    % fig = gcf;
    % fig.PaperPositionMode = 'auto';
    %name=sprintf('correlation_par%d',i_con);
    % print(fig,name,'-dpng','-r0')
    % hold off

    %make plot for different noise levels
    % load('sim_noise025.mat');
    % figure
    % plot(Med_true_r,Med_sim_r);
    % hold on
    % load('sim_noise1.mat');
    % plot(Med_true_r,Med_sim_r);
    % load('sim_noise4.mat');
    % plot(Med_true_r,Med_sim_r);
    % hline(0);
    % hline(abs(corr(DesMat(:,1),DesMat(:,3))));
    % vline(0);
    % axis([-0.7 0.7 -0.7 0.7])
    % hold off
    end;
end;
end;

cd(box_path);
disp('DONE');





% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1
axes(hObject);
imshow('logo.png');