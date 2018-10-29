function [BMS,WMS,JMS,EMS] = ICC_computation(nr_subj,runs,data)

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
    