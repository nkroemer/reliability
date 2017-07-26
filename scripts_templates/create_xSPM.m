function xSPM = create_xSPM(SPM_path,box_dir,p,con_count,varargin)

%% adapted version of spm_getSPM to get xSPM without input
cd(SPM_path);
swd=pwd;
%-Load SPM.mat
%--------------------------------------------------------------------------
try
    load(fullfile(swd,'SPM.mat'));
catch
    error(['Cannot read ' fullfile(swd,'SPM.mat')]);
end
SPM.swd = swd;

cd(SPM.swd);

xX   = SPM.xX;                      %-Design definition structure
XYZ  = SPM.xVol.XYZ;                %-XYZ coordinates
S    = SPM.xVol.S;                  %-search Volume {voxels}
R    = SPM.xVol.R;                  %-search Volume {resels}
M    = SPM.xVol.M(1:3,1:3);         %-voxels to mm matrix
VOX  = sqrt(diag(M'*M))';           %-voxel dimensions

%contrasts
xCon = SPM.xCon; 
Ic        = con_count; % 3 = subjective value  
IcAdd = []; % no new contrasts added
nc        = length(Ic);  % Number of contrasts

%-Allow user to extend the null hypothesis for conjunctions
%
% n: conjunction number
% u: Null hyp is k<=u effects real; Alt hyp is k>u effects real
%    (NB Here u is from Friston et al 2004 paper, not statistic thresh).
%                  u         n
% Conjunction Null nc-1      1     |    u = nc-n
% Intermediate     1..nc-2   nc-u  |    #effects under null <= u
% Global Null      0         nc    |    #effects under alt  > u,  >= u+1
%----------------------------------+---------------------------------------
n = 1;

SPM.xCon = xCon;

    Im = [];
    pm = [];
    Ex = [];

%-Create/Get title string for comparison
%--------------------------------------------------------------------------
str  = xCon(Ic).name;

titlestr = str; 

%-Compute & store contrast parameters, contrast/ESS images, & SPM images
%==========================================================================
SPM.xCon = xCon;
if isnumeric(Im)
    SPM  = spm_contrasts(SPM, unique([Ic, Im, IcAdd]));
else
    SPM  = spm_contrasts(SPM, unique([Ic, IcAdd]));
end
xCon     = SPM.xCon;
STAT     = xCon(Ic(1)).STAT;
VspmSv   = cat(1,xCon(Ic).Vspm);

%-Check conjunctions - Must be same STAT w/ same df
%--------------------------------------------------------------------------
if (nc > 1) && (any(diff(double(cat(1,xCon(Ic).STAT)))) || ...
                any(abs(diff(cat(1,xCon(Ic).eidf))) > 1))
    error('illegal conjunction: can only conjoin SPMs of same STAT & df');
end


%-Degrees of Freedom and STAT string describing marginal distribution
%--------------------------------------------------------------------------
df     = [xCon(Ic(1)).eidf xX.erdf];
if nc > 1
    if n > 1
        str = sprintf('^{%d \\{Ha:k\\geq%d\\}}',nc,(nc-n)+1);
    else
        str = sprintf('^{%d \\{Ha:k=%d\\}}',nc,(nc-n)+1);
    end
else
    str = '';
end

switch STAT
    case 'T'
        STATstr = sprintf('%c%s_{%.0f}','T',str,df(2));
    case 'F'
        STATstr = sprintf('%c%s_{%.0f,%.0f}','F',str,df(1),df(2));
    case 'P'
        if strcmp(SPM.PPM.xCon(Ic).PSTAT,'T')
            STATstr = sprintf('%s^{%0.2f}','PPM',df(1));
        else
            STATstr='PPM';
        end
end


%-Compute (unfiltered) SPM pointlist for masked conjunction requested
%==========================================================================
fprintf('\t%-32s: %30s','SPM computation','...initialising')            %-#


%-Compute conjunction as minimum of SPMs
%--------------------------------------------------------------------------
Z     = Inf;
for i = Ic
    Z = min(Z,spm_data_read(xCon(i).Vspm,'xyz',XYZ));
end


%-Copy of Z and XYZ before masking, for later use with FDR
%--------------------------------------------------------------------------
XYZum = XYZ;
Zum   = Z;

%-Compute mask and eliminate masked voxels
%--------------------------------------------------------------------------
for i = 1:numel(Im)
    
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...masking')           %-#
    if isnumeric(Im)
        Mask = spm_data_read(xCon(Im(i)).Vspm,'xyz',XYZ);
        um   = spm_u(pm,[xCon(Im(i)).eidf,xX.erdf],xCon(Im(i)).STAT);
        if Ex
            Q = Mask <= um;
        else
            Q = Mask >  um;
        end
    else
        v = spm_data_hdr_read(Im{i});
        Mask = spm_data_read(v,'xyz',v.mat\SPM.xVol.M*[XYZ; ones(1,size(XYZ,2))]);
        Q = Mask ~= 0 & ~isnan(Mask);
        if Ex, Q = ~Q; end
    end
    XYZ   = XYZ(:,Q);
    Z     = Z(Q);
    if isempty(Q)
        fprintf('\n')                                                   %-#
        sw = warning('off','backtrace');
        warning('SPM:NoVoxels','No voxels survive masking at p=%4.2f',pm);
        warning(sw);
        break
    end
end


%==========================================================================
% - H E I G H T   &   E X T E N T   T H R E S H O L D S
%==========================================================================

u   = -Inf;        % height threshold
k   = 0;           % extent threshold {voxels}

%-Get FDR mode
%--------------------------------------------------------------------------
try
    topoFDR = spm_get_defaults('stats.topoFDR');
catch
    topoFDR = true;
end

if  spm_mesh_detect(xCon(Ic(1)).Vspm)
    G = export(gifti(SPM.xVol.G),'patch');
end

%-Height threshold - classical inference
%--------------------------------------------------------------------------
   %-Get height threshold
    %----------------------------------------------------------------------
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...height threshold')  %-#
 
    %% enable possibility for user to choose between thresholds
    thresDesc = 'none';
    % thresDesc = 'FWE';
    % thresDesc = 'FDR';
%     switch thresDesc
%         case 'FWE' % Family-wise false positive rate
%             %--------------------------------------------------------------
%             try
%                 u = xSPM.u;
%             catch
%                 u = spm_input('p value (FWE)','+0','r',0.05,1,[0,1]);
%             end
%             thresDesc = ['p<' num2str(u) ' (' thresDesc ')'];
%             u = spm_uc(u,df,STAT,R,n,S);
%             
%             
%         case 'FDR' % False discovery rate
%             %--------------------------------------------------------------
%             if topoFDR,
%                 fprintf('\n');                                          %-#
%                 error('Change defaults.stats.topoFDR to use voxel FDR');
%             end
%             try
%                 u = xSPM.u;
%             catch
%                 u = spm_input('p value (FDR)','+0','r',0.05,1,[0,1]);
%             end
%             thresDesc = ['p<' num2str(u) ' (' thresDesc ')'];
%             u = spm_uc_FDR(u,df,STAT,n,VspmSv,0);
%             
%         case 'none'  % No adjustment: p for conjunctions is p of the conjunction SPM
            %--------------------------------------------------------------
            u = p;
            if u <= 1
                thresDesc = ['p<' num2str(u) ' (unc.)'];
                u = spm_u(u^(1/n),df,STAT);
            else
                thresDesc = [STAT '=' num2str(u) ];
            end
                        
%         otherwise
%             %--------------------------------------------------------------
%             fprintf('\n');                                              %-#
%             error('Unknown control method "%s".',thresDesc);
%             
%     end % switch thresDesc
    
    %-Compute p-values for topological and voxel-wise FDR (all search voxels)
    %----------------------------------------------------------------------
    if ~topoFDR
        fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...for voxelFDR')  %-#
        Ps = spm_z2p(Zum,df,STAT,n);
    end
    
    %-Peak FDR
    %----------------------------------------------------------------------
    if ~spm_mesh_detect(xCon(Ic(1)).Vspm)
        [up,Pp] = spm_uc_peakFDR(0.05,df,STAT,R,n,Zum,XYZum,u);
    else
        [up,Pp] = spm_uc_peakFDR(0.05,df,STAT,R,n,Zum,XYZum,u,G);
    end
    
    %-Cluster FDR
    %----------------------------------------------------------------------
    if n == 1 %% && STAT == 'T'
        if ~spm_mesh_detect(xCon(Ic(1)).Vspm)
            V2R        = 1/prod(SPM.xVol.FWHM(SPM.xVol.DIM > 1));
            [uc,Pc,ue] = spm_uc_clusterFDR(0.05,df,STAT,R,n,Zum,XYZum,V2R,u);
        else
            V2R        = 1/prod(SPM.xVol.FWHM);
            [uc,Pc,ue] = spm_uc_clusterFDR(0.05,df,STAT,R,n,Zum,XYZum,V2R,u,G);
        end
    else
        uc  = NaN;
        ue  = NaN;
        Pc  = [];
    end
    
    %-Peak FWE
    %----------------------------------------------------------------------
    uu      = spm_uc(0.05,df,STAT,R,n,S);
    
    

%-Calculate height threshold filtering
%--------------------------------------------------------------------------
if spm_mesh_detect(xCon(Ic(1)).Vspm), str = 'vertices'; else str = 'voxels'; end
Q      = find(Z > u);

%-Apply height threshold
%--------------------------------------------------------------------------
Z      = Z(:,Q);
XYZ    = XYZ(:,Q);
if isempty(Q)
    fprintf('\n');                                                      %-#
    sw = warning('off','backtrace');
    warning('SPM:NoVoxels','No %s survive height threshold at u=%0.2g',str,u);
    warning(sw);
end


%-Extent threshold
%--------------------------------------------------------------------------
if ~isempty(XYZ)
    
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...extent threshold'); %-#
    
    %-Get extent threshold [default = 0]
    %----------------------------------------------------------------------
 k = 0; %% all voxels should be included!!!
    
    %-Calculate extent threshold filtering
    %----------------------------------------------------------------------
    if  ~spm_mesh_detect(xCon(Ic(1)).Vspm)
        A = spm_clusters(XYZ);
    else
        T = false(SPM.xVol.DIM');
        T(XYZ(1,:)) = true;
        A = spm_mesh_clusters(G,T)';
        A = A(XYZ(1,:));
    end
    Q     = [];
    for i = 1:max(A)
        j = find(A == i);
        if length(j) >= k, Q = [Q j]; end
    end
    
    % ...eliminate voxels
    %----------------------------------------------------------------------
    Z     = Z(:,Q);
    XYZ   = XYZ(:,Q);
    if isempty(Q)
        fprintf('\n');                                                  %-#
        sw = warning('off','backtrace');
        warning('SPM:NoVoxels','No %s survive extent threshold at k=%0.2g',str,k);
        warning(sw);
    end
    
else
    
    k = 0;
    
end % (if ~isempty(XYZ))

%==========================================================================
% - E N D
%==========================================================================
fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done')                %-#
spm('Pointer','Arrow')

%-Assemble output structures of unfiltered data
%==========================================================================
xSPM   = struct( ...
            'swd',      swd,...
            'title',    titlestr,...
            'Z',        Z,...
            'n',        n,...
            'STAT',     STAT,...
            'df',       df,...
            'STATstr',  STATstr,...
            'Ic',       Ic,...
            'Im',       {Im},...
            'pm',       pm,...
            'Ex',       Ex,...
            'u',        u,...
            'k',        k,...
            'XYZ',      XYZ,...
            'XYZmm',    SPM.xVol.M(1:3,:)*[XYZ; ones(1,size(XYZ,2))],...
            'S',        SPM.xVol.S,...
            'R',        SPM.xVol.R,...
            'FWHM',     SPM.xVol.FWHM,...
            'M',        SPM.xVol.M,...
            'iM',       SPM.xVol.iM,...
            'DIM',      SPM.xVol.DIM,...
            'VOX',      VOX,...
            'Vspm',     VspmSv,...
            'thresDesc',thresDesc);

%-RESELS per voxel (density) if it exists
%--------------------------------------------------------------------------
try, xSPM.VRpv = SPM.xVol.VRpv; end
try
    xSPM.units = SPM.xVol.units;
catch
    try, xSPM.units = varargin{1}.units; end
end

%-Topology for surface-based inference
%--------------------------------------------------------------------------
if spm_mesh_detect(xCon(Ic(1)).Vspm)
    xSPM.G     = G;
    xSPM.XYZmm = xSPM.G.vertices(xSPM.XYZ(1,:),:)';
end

%-p-values for topological and voxel-wise FDR
%--------------------------------------------------------------------------
try, xSPM.Ps   = Ps;             end  % voxel   FDR
try, xSPM.Pp   = Pp;             end  % peak    FDR
try, xSPM.Pc   = Pc;             end  % cluster FDR

%-0.05 critical thresholds for FWEp, FDRp, FWEc, FDRc
%--------------------------------------------------------------------------
try, xSPM.uc   = [uu up ue uc];  end

cd(box_dir);