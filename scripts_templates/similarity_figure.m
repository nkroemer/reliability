function [fig] = similarity_figure(r_mat,sim2mean,run1,run2,str)

fig = figure;
set(fig,'units','normalized','position',[0.25 0 0.50 1]);
%color matrix
subplot(2,2,1:2);
imagesc(r_mat,[0,1]);
colormap('inferno');
caxis([0,1]);
colorbar;
if sim2mean == 1
    xlabel(sprintf('subjects session %d (last column similarity to mean image)',run2));
else
    xlabel(sprintf('subjects session %d ',run2));
end;
ylabel(sprintf('subjects session %d',run1));
name1 = sprintf('similarity-%s-%d-%d',str,run1,run2);
title(name1);

% histograms
subplot(2,2,3);
triu_r_mat = triu(r_mat,1);
vec_off_diag_r_mat = triu_r_mat(triu_r_mat~=0);
h1 = histogram(diag(r_mat));
hold on;
h2 = histogram(vec_off_diag_r_mat);
h1.Normalization = 'probability';
h1.BinWidth = 0.1;
h2.Normalization = 'probability';
h2.BinWidth = 0.1;
xlabel('similarity');
ylabel('frequency in percentage');
title(sprintf('histogram-%s-%d-%d',str,run1,run2));

               
% ecdf - densitiy plots
subplot(2,2,4);
[f,x,flo,fup]=ecdf(diag(r_mat),'bounds','on');
plot(x,f,'LineWidth',2,'Color','blue')
hold on
plot(x,flo,'LineWidth',1,'Color','blue','LineStyle','--')
hold on
plot(x,fup,'LineWidth',1,'Color','blue','LineStyle','--')
hold on;
clear x f flo fup
[f,x,flo,fup]=ecdf(vec_off_diag_r_mat,'bounds','on');
plot(x,f,'LineWidth',2,'Color','red')
hold on
plot(x,flo,'LineWidth',1,'Color','red','LineStyle','--')
hold on
plot(x,fup,'LineWidth',1,'Color','red','LineStyle','--')
title(sprintf('cumulative-density-%s-%d-%d',str,run1,run2));





