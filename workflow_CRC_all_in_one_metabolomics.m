%% load additional scripts

addpath(genpath('additionalScripts'))

%% load metabolomics data (UK cohort)

C=readmatrix('data/CRC_metabolomics_data_anon.xlsx','Sheet','covariates','Range','B2:S161');
C_pid=readcell('data/CRC_metabolomics_data_anon.xlsx','Sheet','covariates','Range','A2:A161');
Clabs=readcell('data/CRC_metabolomics_data_anon.xlsx','Sheet','covariates','Range','B1:S1');
i_Tumour=readmatrix('data/CRC_metabolomics_data_anon.xlsx','Sheet','covariates','Range','T2:T161')==1;
M_data=readmatrix('data/CRC_metabolomics_data_anon.xlsx','Sheet','metabolomics_data','Range','B2:BD161');
M_pid=readcell('data/CRC_metabolomics_data_anon.xlsx','Sheet','metabolomics_data','Range','A2:A161');
M_labs=readcell('data/CRC_metabolomics_data_anon.xlsx','Sheet','metabolomics_labels','Range','A2:A56');

%% settings for analysis
alp=0.05; % alpha level

%% perform metabolome clustering - partial correlation - complete case analysis

% partial correlation
[r,p]=partialcorr(M_data(i_Tumour,:),C(i_Tumour,:),'Type','Spearman','Rows','pairwise');
nr=size(r,1);

% permutation testing (n>=100), scramble data within with variable randomly, then retest and find those that have less than 5% permuted correlations that are higher than the calculated correlation (also gets rid of negative correlations)
rng(1,'twister');
nrrb=r*0;
nperm=1000;
hw=waitbar(0,'Permutations to go');
for k=1:nperm
    waitbar(k/nperm,hw)
    Mr=M_data(i_Tumour,:);
    Cr=C(i_Tumour,:); % can this be kept constant as Mr is resampled already
    for j=1:nr
        Mr(:,j)=randsample(Mr(:,j),sum(i_Tumour));
    end
    rr=partialcorr(Mr,Cr,'Type','Spearman','Rows','pairwise');
    nrrb=nrrb+double(rr>r).*double(r>0)+double(rr<r).*double(r<0);
end
close(hw)
nrrb=nrrb/nperm;
r2=r;
r2(nrrb>alp)=0;
clear Mr Cr rr hw j

% cluster the data
pd=pdist(r2,'correlation');
tree=linkage(pd,'average');
figure;[~,~,op]=dendrogram(tree,0);
yl=[0 max(ylim)];
dt=flipud([0;tree(:,3);max(yl)]);
dt=(dt(1:(end-1))+dt(2:end))/2;
close;

% calculate modularity
[m,~,~,mb]=JMP_modularity(r2,0,tree,op);
m=m';
pval=ones(nr,1);
for k=1:nr
    [~,pval(k,1)]=ttest2(mb(:,k),m(k),'Tail','left');
end
pval(isnan(pval))=0;
pval=bhfdr(pval);
th=find(LocalMaxMin(m)==1 & pval==min(pval),1,'first'); % first local maximum that is most significantly different from random
if isempty(th)
    th=imax(m);th(pval(th)>alp)=[];
    if isempty(th)
        th=find(LocalMaxMin(m)==1 & pval<=alp,1,'first'); % first local maximum that is significantly different from random
        if isempty(th)
            th=imin(pval);
            if pval(th)>alp
                th=1;
            end
        end
    end
end
ci=cluster(tree,'maxclust',th);
cl=find(diff(ci(op))~=0)+0.5;

% create figure
np=5;
pmat=reshape(1:np^2,np,np)';
thCol=[1 1 1]*0.5;

f1=figure(3022);

pm=pmat(1:2,3:5)';h(2)=subplot(np,np,pm(:));
dendrogram(tree,0,'ColorThreshold',tree(nr-th,3)+eps);
dc=get(h(2),'Children');
for k=1:size(dc,1)
    dcc=get(dc(k),'Color');
    if dcc(1)~=dcc(2) || dcc(1)~=dcc(3) || dcc(2)~=dcc(3)
        set(dc(k),'LineWidth',2)
    end
end
hold on;plot(xlim,[1 1]*dt(th),'--','Color',thCol);hold off;
ylim(yl)
ylabel('Partial Spearman correlation distance')
set(gca,'XTickLabel',[],'YAxisLocation','Right')

pm=pmat(3:5,3:5)';h(4)=subplot(np,np,pm(:));
hold on;
imagesc(r2(op,op));cmap
for k=1:numel(cl)
    plot([1 1]*cl(k),[0.5 nr+0.5],'k--')
    plot([0.5 nr+0.5],[1 1]*cl(k),'k--')
end
xlim([0.5 nr+0.5])
ylim([0.5 nr+0.5])
hold off
set(gca,'XTick',1:numel(op),'XTickLabel',M_labs(op,end),'YTick',1:numel(op),'YTickLabel',M_labs(op,end),'FontSize',6)
XTickRotateJMP(M_labs(op,end),67.5)

pm=pmat(3:5,1)';h(4)=subplot(np,np,pm(:));
cmap
hc=colorbar;
axis off
set(get(hc,'YLabel'),'String','Partial Spearman correlation')

pm=pmat(1:2,1:2)';h(1)=subplot(np,np,pm(:));
hold on;
plot(m,dt,'g-','LineWidth',2);
plot(prctile(mb,0)',dt,'r-');
plot(xlim,[1 1]*dt(th),'--','Color',thCol);
plot(prctile(mb,0:2.5:100)',dt,'r-');
plot(m,dt,'g-','LineWidth',2);
set(gca,'XDir','reverse','XAxisLocation','Top')
xlim(lim([mb(:);m(:)]).*[1 1.01])
ylim(yl)
plot(xlim,[1 1]*dt(th),'--','Color',thCol)
hold off
xlabel('Modularity')
legend({'Actual clustering','Random clustering','Optimal number of clusters'},'Location','SouthWest')
set(gcf,'Position',[1372,1055,1300,1000])

%% Comparing metabolomics in T vs TPN

Mc_data=M_data;
Mcc_data=zeros(size(Mc_data,1),numel(unique(ci)));

LV_scal_param={};
for k=1:size(Mcc_data,2)
    [tr,~,mc,va]=JMP_scale(Mc_data(:,ci==k),[],1);
    LV_scal_param{k,1}=mc;
    LV_scal_param{k,2}=va;
    [U,D,V]=recursivePCA(tr,1);
    Mcc_data(:,k)=U*D;
    LV_scal_param{k,3}=V;
end

pvals=zeros(size(Mcc_data,2),1);
zval=zeros(size(Mcc_data,2),1);
hgrp={};
for k=1:size(Mcc_data,2)
    [pvals(k),~,stats]=signrank(Mcc_data(i_Tumour,k),Mcc_data(~i_Tumour,k));
    zval(k)=stats.zval;
    if median(Mcc_data(i_Tumour,k))>median(Mcc_data(~i_Tumour,k))
        hgrp{k,1}='T'; %#ok
    elseif median(Mcc_data(i_Tumour,k))<median(Mcc_data(~i_Tumour,k))
        hgrp{k,1}='TPN'; %#ok
    else
        hgrp{k,1}=' '; %#ok
    end
end
pord=unique(ci(op),'rows','stable');
pvals=pvals(pord);
pfdr=bhfdr(pvals);
zval=zval(pord);
hgrp=hgrp(pord);

gc=get(f1,'Children');
gc=gc(end);
H=get(gc,'Children');H(1)=[];
xcol=zeros(numel(ci),3);
for k=1:size(H,1)
    xd=get(H(k),'XData');
    yd=get(H(k),'YData');
    yd=xd(yd==yl(1));
    if ~isempty(yd)
        xcol(yd,:)=repmat(get(H(k),'Color'),numel(yd),1);
    end
end

xcol(sum(xcol,2)==0,1)=NaN;
xcol=unique(xcol,'rows','stable');
xcol(isnan(xcol(:,1)),1)=0;

f2=figure(3021);hold on;
for k=1:numel(pvals)
    patch([-1 1 1 -1 -1]*0.45+k,[0 0 -log10(pfdr(k))*sign(zval(k))*[1 1] 0],[1 1 1 1 1],'FaceColor',xcol(k,:),'EdgeColor',[0 0 0])
end
xlim([0.4 k+0.5])
plot(xlim,-log10(0.05)*[1 1],'--','Color',[1 1 1]*0.5)
plot(xlim,log10(0.05)*[1 1],'--','Color',[1 1 1]*0.5)
set(gca,'XTick',1:k)
xlabel('Cluster')
ylabel('-^1^0log p_F_D_R \times sign(z)')
title('Signed rank test for each cluster (T vs TPN)')

%% load metabolomics data (Czech cohort)

C=readmatrix('data/KRCA_metabolomics_data_anon.xlsx','Sheet','covariates','Range','B2:I53'); %
C_pid=readcell('data/KRCA_metabolomics_data_anon.xlsx','Sheet','covariates','Range','A2:A53');
Clabs=readcell('data/KRCA_metabolomics_data_anon.xlsx','Sheet','covariates','Range','B1:I1'); %
i_Tumour=readmatrix('data/KRCA_metabolomics_data_anon.xlsx','Sheet','covariates','Range','J2:J53')==1; %
M_data=readmatrix('data/KRCA_metabolomics_data_anon.xlsx','Sheet','metabolomics_data','Range','B2:BD53');
M_pid=readcell('data/KRCA_metabolomics_data_anon.xlsx','Sheet','metabolomics_data','Range','A2:A53');
M_labs=readcell('data/KRCA_metabolomics_data_anon.xlsx','Sheet','metabolomics_labels','Range','A2:A56');

%% perform metabolomics clustering - partial correlation - complete case analysis

% partial correlation
[r,p]=partialcorr(M_data(i_Tumour,:),C(i_Tumour,:),'Type','Spearman','Rows','pairwise');
nr=size(r,1);

% permutation testing (n>=100), scramble data within with variable randomly, then retest and find those that have less than 5% permuted correlations that are higher than the calculated correlation (also gets rid of negative correlations I think)
rng(1,'twister');
nrrb=r*0;
nperm=1000;
hw=waitbar(0,'Permutations to go');
for k=1:nperm
    waitbar(k/nperm,hw)
    Mr=M_data(i_Tumour,:);
    Cr=C(i_Tumour,:); % can this be kept constant as Mr is resampled already
    for j=1:nr
        Mr(:,j)=randsample(Mr(:,j),sum(i_Tumour));
    end
    rr=partialcorr(Mr,Cr,'Type','Spearman','Rows','pairwise');
    nrrb=nrrb+double(rr>r).*double(r>0)+double(rr<r).*double(r<0);
end
close(hw)
nrrb=nrrb/nperm;
r2=r;
r2(nrrb>alp)=0;
clear Mr Cr rr hw j

% calculate modularity - test UK clustering on Czech data
[m,~,~,mb]=JMP_modularity(r2,0,tree,op);
m=m';
pval=ones(nr,1);
for k=1:nr
    [~,pval(k,1)]=ttest2(mb(:,k),m(k),'Tail','left');
end
pval(isnan(pval))=0;
pval=bhfdr(pval);

% create figure
np=5;
pmat=reshape(1:np^2,np,np)';
thCol=[1 1 1]*0.5;

f3=figure(3020);

pm=pmat(1:2,3:5)';h(2)=subplot(np,np,pm(:));
dendrogram(tree,0,'ColorThreshold',tree(nr-th,3)+eps);
dc=get(h(2),'Children');
for k=1:size(dc,1)
    dcc=get(dc(k),'Color');
    if dcc(1)~=dcc(2) || dcc(1)~=dcc(3) || dcc(2)~=dcc(3)
        set(dc(k),'LineWidth',2)
    end
end
hold on;plot(xlim,[1 1]*dt(th),'--','Color',thCol);hold off;
ylim(yl)
ylabel('Partial Spearman correlation distance')
title('Clustering obtained from UK data applied on Czech data')
set(gca,'XTickLabel',[],'YAxisLocation','Right')

pm=pmat(3:5,3:5)';h(4)=subplot(np,np,pm(:));
hold on;
imagesc(r2(op,op));cmap
for k=1:numel(cl)
    plot([1 1]*cl(k),[0.5 nr+0.5],'k--')
    plot([0.5 nr+0.5],[1 1]*cl(k),'k--')
end
xlim([0.5 nr+0.5])
ylim([0.5 nr+0.5])
hold off
set(gca,'XTick',1:numel(op),'XTickLabel',M_labs(op,end),'YTick',1:numel(op),'YTickLabel',M_labs(op,end),'FontSize',6)
XTickRotateJMP(M_labs(op,end),67.5)

pm=pmat(3:5,1)';h(4)=subplot(np,np,pm(:));
cmap
hc=colorbar;
axis off
set(get(hc,'YLabel'),'String','Partial Spearman correlation')

pm=pmat(1:2,1:2)';h(1)=subplot(np,np,pm(:));
hold on;
plot(m,dt,'g-','LineWidth',2);
plot(prctile(mb,0)',dt,'r-');
plot(xlim,[1 1]*dt(th),'--','Color',thCol);
plot(prctile(mb,0:2.5:100)',dt,'r-');
plot(m,dt,'g-','LineWidth',2);
set(gca,'XDir','reverse','XAxisLocation','Top')
xlim(lim([mb(:);m(:)]).*[1 1.01])
ylim(yl)
plot(xlim,[1 1]*dt(th),'--','Color',thCol)
hold off
xlabel('Modularity')
legend({'Actual clustering (Czech data)','Random clustering (Czech data)','Optimal number of clusters (UK data)'},'Location','SouthWest')
set(gcf,'Position',[1372,1055,1300,1000])

%% load metabolomics data (Czech cohort)

K_pid=readcell('data/KRCA_metabolomics_data_paired_anon.xlsx','Sheet','covariates','Range','A2:A11');
i_Tumour=readmatrix('data/KRCA_metabolomics_data_paired_anon.xlsx','Sheet','covariates','Range','B2:B11')==1;
M_data=readmatrix('data/KRCA_metabolomics_data_paired_anon.xlsx','Sheet','metabolomics_data','Range','B2:BD11');
M_pid=readcell('data/KRCA_metabolomics_data_paired_anon.xlsx','Sheet','metabolomics_data','Range','A2:A11');
M_labs=readcell('data/KRCA_metabolomics_data_paired_anon.xlsx','Sheet','metabolomics_labels','Range','A2:A56');

%% Comparing metabolomics in T vs TPN - Czech cohort

Mc_data=M_data;
nte=size(Mc_data,1);
Mcc_data=zeros(nte,numel(unique(ci)));

for k=1:size(Mcc_data,2)
    te=Mc_data(:,ci==k)-LV_scal_param{k,1}(ones(nte,1),:);
    te=te./LV_scal_param{k,2}(ones(nte,1),:);
    Mcc_data(:,k)=te*LV_scal_param{k,3};
end

signifUK=find(pfdr<alp);

pvals=zeros(size(Mcc_data,2),1);
zval=zeros(size(Mcc_data,2),1);
hgrp={};
for k=1:size(Mcc_data,2)
    [pvals(k),~,stats]=signrank(Mcc_data(i_Tumour,k),Mcc_data(~i_Tumour,k),'method','approximate');
    zval(k)=stats.zval;
    if median(Mcc_data(i_Tumour,k))>median(Mcc_data(~i_Tumour,k))
        hgrp{k,1}='T'; %#ok
    elseif median(Mcc_data(i_Tumour,k))<median(Mcc_data(~i_Tumour,k))
        hgrp{k,1}='TPN'; %#ok
    else
        hgrp{k,1}=' '; %#ok
    end
end
pord=unique(ci(op),'rows','stable');
pvals=pvals(pord);
pfdr=bhfdr(pvals);
zval=zval(pord);
hgrp=hgrp(pord);

gc=get(f1,'Children');
gc=gc(end);
H=get(gc,'Children');H(1)=[];
xcol=zeros(numel(ci),3);
for k=1:size(H,1)
    xd=get(H(k),'XData');
    yd=get(H(k),'YData');
    yd=xd(yd==yl(1));
    if ~isempty(yd)
        xcol(yd,:)=repmat(get(H(k),'Color'),numel(yd),1);
    end
end

xcol(sum(xcol,2)==0,1)=NaN;
xcol=unique(xcol,'rows','stable');
xcol(isnan(xcol(:,1)),1)=0;

f4=figure(3019);hold on;
for k=1:numel(pvals)
    patch([-1 1 1 -1 -1]*0.45+k,[0 0 -log10(pfdr(k))*sign(zval(k))*[1 1] 0],[1 1 1 1 1],'FaceColor',xcol(k,:),'EdgeColor',[0 0 0])
end
xlim([0.4 k+0.5])
plot(xlim,-log10(0.05)*[1 1],'--','Color',[1 1 1]*0.5)
plot(xlim,log10(0.05)*[1 1],'--','Color',[1 1 1]*0.5)
set(gca,'XTick',1:k)
xlabel('Cluster')
ylabel('-^1^0log p_F_D_R \times sign(z)')
title('Signed rank test for each cluster (T vs TPN, Czech data)')

% only significant ones from UK
f5=figure(3018);hold on;
for k=1:numel(pvals)
    if sum(signifUK==k)==1
        patch([-1 1 1 -1 -1]*0.45+k,[0 0 -log10(pvals(k))*sign(zval(k))*[1 1] 0],[1 1 1 1 1],'FaceColor',xcol(k,:),'EdgeColor',[0 0 0])
    end
end
xlim([0.4 k+0.5])
plot(xlim,-log10(0.05)*[1 1],'--','Color',[1 1 1]*0.5)
plot(xlim,log10(0.05)*[1 1],'--','Color',[1 1 1]*0.5)
set(gca,'XTick',1:k)
xlabel('Cluster')
ylabel('-^1^0log p \times sign(z)')
title('Signed rank test for each significant UK cluster (T vs TPN, Czech data)')

