%%script for calculating split-half reliability of fMRI data
load ('vp_gesamt.mat');
%---modify me
run = 1;
for i = 1:length(vp)
%%---loads and modifies SPM for each participant
    dir_spm = sprintf('M:\SeSyN\019\BMBF_itech\Juliane\niftii\%06d\stats_%d_maskfix',vp(i),run);
    cd(dir_spm);
    load('SPM.mat');
    
    %split onsets
    onsets = SPM.Sess.U(1).ons;
    pos = randperm(length(onsets))';
    onsets1 = onsets(pos(1:end/2));
    onsets2 = onsets(pos(end/2:end));
    
    %split parametric value
    sv = SPM.Sess.U(1).P.P;
    sv1 = sv(pos(1:end/2));
    sv2 = sv(pos(end/2:end));
    
    %create new conditions
    SPM.Sess.U(5).name = 'half1';
    SPM.Sess.U(5).ons = onsets1;
    SPM.Sess.U(5).dur = zeros(length(onsets1),1);
    SPM.Sess.U(5).orth = 1;
    SPM.Sess.U(5).P.name = 'value';
    SPM.Sess.U(5).P.P = sv1;
    SPM.Sess.U(5).P.h = 1;
    SPM.Sess.U(5).P.i = [1,2];
    SPM.Sess.U(5).dt = SPM.Sess.U(1).dt;
    SPM.Sess.U(5).u = SPM.Sess.U(1).u;
    SPM.Sess.U(5).pst = SPM.Sess.U(1).pst;
    
    SPM.Sess.U(6).name = 'half2';
    SPM.Sess.U(6).ons = onsets2;
    SPM.Sess.U(6).dur = zeros(length(onsets2),1);
    SPM.Sess.U(6).orth = 1;
    SPM.Sess.U(6).P.name = 'value';
    SPM.Sess.U(6).P.P = sv2;
    SPM.Sess.U(6).P.h = 1;
    SPM.Sess.U(6).P.i = [1,2];
    SPM.Sess.U(6).dt = SPM.Sess.U(1).dt;
    SPM.Sess.U(6).u = SPM.Sess.U(1).u;
    SPM.Sess.U(6).pst = SPM.Sess.U(1).pst;
    
%% Convolve predictors with HRF
%-----half1
s = 1;
bf = SPM.xBF.bf;
V  = SPM.xBF.Volterra;
U = SPM.Sess.U(5);
[X,Xn,Fc] = spm_Volterra(U,bf,V);
k = SPM.nscan;
fMRI_T  = SPM.xBF.T;
fMRI_T0 = SPM.xBF.T0;
X  = X((0:(k - 1))*fMRI_T + fMRI_T0 + 32,:);
% apply high-pass filter + auto-decorrelate
xX = SPM.xX;
W  = SPM.xX.W;
xX.X(:,1:2) = X;
xX.xKXs   = spm_sp('Set',spm_filter(xX.K,W*xX.X));       % KWX
xX.xKXs.X = full(xX.xKXs.X);
HRFmK = xX.xKXs.X(:,1:2);
% regress BOLD signal on onsets 
%XX = [ones(length(Ybold),1) HRFmK];
%l = l - 1/2*log(2*pi*sigma^2) - 1/2*sum((Ybold-XX*beta').^2)/sigma^2;
%ll = -l;
    
%-----half2
s = 1;
bf = SPM.xBF.bf;
V  = SPM.xBF.Volterra;
U = SPM.Sess.U(6);
[X,Xn,Fc] = spm_Volterra(U,bf,V);
k = SPM.nscan;
fMRI_T  = SPM.xBF.T;
fMRI_T0 = SPM.xBF.T0;
X  = X((0:(k - 1))*fMRI_T + fMRI_T0 + 32,:);
% apply high-pass filter + auto-decorrelate
xX = SPM.xX;
W  = SPM.xX.W;
xX.X(:,1:2) = X;
xX.xKXs   = spm_sp('Set',spm_filter(xX.K,W*xX.X));       % KWX
xX.xKXs.X = full(xX.xKXs.X);
HRFmK = xX.xKXs.X(:,1:2);
% regress BOLD signal on onsets 
%XX = [ones(length(Ybold),1) HRFmK];
%l = l - 1/2*log(2*pi*sigma^2) - 1/2*sum((Ybold-XX*beta').^2)/sigma^2;
%ll = -l;

%% Calculate split-half reliability

end;
    
    