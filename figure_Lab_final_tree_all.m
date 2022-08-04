
load GM_tree_catching_error_all M G
% lblue = [1 1 1]*0.7;
% lred = [1 1 1]*0.7;

lblue = [0 53 200]./255;
lred = lblue;
% lblue = [137 207 240]./255;
% lred = [137 207 240]./255;

% lblue = [66 224 209]./255;
% lred = [66 224 209]./255;

load AGE age

Age = age;
NTb = [222 226 231 235 251 253 232 116 263 258 266 265 271 280 279 276 282];
NTg = [223 230 234 239 250 254 256 257 259 244 252 229 264];
% Valid = isfinite(Age) & isfinite(Bmat') & Age>7 & Bmat'>0 ...
%     & Age<13;
Valid = [NTb NTg];
figure
ASDb = [218 20 21 240 267 212 269 215 270 274];
ASDg = 233;

oldNTb = [222 226 253 258 266 279 280 271 232 116 276 282];
oldASDb = [215 21 218 240 267 269 270 274 20];

adultSetb = [2 7:9 15:16 17 500];
adultSetg = [1 3:4 6 18 19];
xOffset = median(age([adultSetb adultSetg]));

% for tree = 1:3,

for subj = [NTb NTg adultSetb adultSetg],
    Mselected = M{subj}(M{subj}~=0 & abs(M{subj})<mean(abs(M{subj}(M{subj}~=0)))+3*std(abs(M{subj}(M{subj}~=0))))*1000*1.1;
    %     Mselected = M{subj}(M{subj}~=0 & abs(M{subj})<0.30)*1000*1.1;
    
    if length(Mselected)>=3,
        Tree(subj,1) = median(abs(Mselected),'omitnan');
        %             Tree(subj,tree) = std(Mselected,'omitnan');
                            TreeSD(subj,1) = std(Mselected,'omitnan')./sqrt(length(Mselected));

    else
        Tree(subj,1) = NaN;
    end
    
    
end
% end

subjsel = [NTb NTg adultSetb adultSetg]; %MR20220721

Bmat = Tree;
incCr = mean(Bmat(Valid))+3*std(Bmat(Valid));
exClu = find(Bmat(Valid)>incCr)
AD = Bmat([adultSetb adultSetg]);
meanAD = mean(AD);
sdAD = std(AD);
Valid = setdiff(Valid, Valid(exClu));

for ag = 8:2:12,
    
    maSubj = Valid(Age(Valid)>=ag-1 & Age(Valid)<ag+1);
    maM(ag) = mean(Bmat(maSubj));
    tsma = tinv([0.025  0.975],length(maSubj)-1);
    SEMma = std(Bmat(maSubj))./sqrt(length(maSubj));
    maCI(ag) = tsma(2)*SEMma;
    %[~,tT(ag)] = ttest2(AD,Bmat(maSubj));
    [tT2(ag), ~] = ranksum(AD,Bmat(maSubj));
    [~,tT(ag),~,STAT{ag}] = ttest2(Bmat(maSubj),AD);
    %[tT(ag), ~] = ranksum(AD,Bmat(maSubj));
    
    L(ag) = length(Bmat(maSubj));
    cohenD(ag) = (maM(ag)-meanAD)./(sqrt(((L(ag)-1)*std(Bmat(maSubj)).^2 + (length(AD)-1)*sdAD.^2)/(L(ag)+length(AD)-2)));
    
    hold on
    line([ag ag],[maM(ag)-maCI(ag) maM(ag)+maCI(ag)],'color',[160 82 45]/255,'linewidth',5)
end
ag = 8:2:12;
tT(ag)
%patch([8:12 12:-1:8],[maM(8:12)+maCI(8:12) maM(12:-1:8)-maCI(12:-1:8)],[1 0.7 0.7]*0.9,'linestyle','none')
hold on
% plot(age([NTg adultSetg]),Bmat([NTg adultSetg]),'^','color',[1 1 1]*0.7,'linewidth',2,'markersize',12,'markerfacecolor',[1 1 1]*0.7,'markeredgecolor','none')
% plot(age([NTb adultSetb]),Bmat([NTb adultSetb]),'o','color',[1 1 1]*0.7,'linewidth',2,'markersize',12,'markerfacecolor',[1 1 1]*0.7,'markeredgecolor','none')

plot(age([NTb adultSetb]),Bmat([NTb adultSetb]),'o','color',lblue,'linewidth',2,'markersize',12,'markerfacecolor',lblue,'markeredgecolor','none')
plot(age([NTg adultSetg]),Bmat([NTg adultSetg]),'o','color',lred,'linewidth',2,'markersize',12,'markerfacecolor',lred,'markeredgecolor','none')

errorbar(age([NTb adultSetb]),Bmat([NTb adultSetb]),TreeSD([NTb adultSetb]),'linestyle','none','color',lblue,'linewidth',1)
errorbar(age([NTg adultSetg]),Bmat([NTg adultSetg]),TreeSD([NTg adultSetg]),'linestyle','none','color',lblue,'linewidth',1)

% plot(age([NTb adultSetb NTg adultSetg]),Bmat([NTb adultSetb NTg adultSetg]),'o','color',[1 1 1]*0.5,'linewidth',2,'markersize',10,'markerfacecolor','none','markeredgecolor',[1 1 1]*0.5)




[r p]=corr(age([NTb NTg])',Bmat([NTb NTg]))

%load Survey_318.mat
%Age = Age(1:317);

fcn = 'a*x +b';
fcnexp = 'a*exp(-b*x) +c';
[F GoF] = fit(Age(Valid)',Bmat(Valid),fcn,'start',[-1 1]);
[Fexp GoFexp] = fit(Age(Valid)',Bmat(Valid),fcnexp,'start',[1500 1/3 20]);

% figure,
% plot(Age(Valid),Bmat(Valid),'o')
Coef = coeffvalues(F);
a = Coef(1);
b = Coef(2);

Coef = coeffvalues(Fexp);
aexp = Coef(1);
bexp = Coef(2);
cexp = Coef(3);

CI = confint(F);
aU = CI(1,1);
bU = CI(1,2);
aL = CI(2,1);
bL = CI(2,2);


x = 7:0.1:13;
y = a*x +b;
yexp = aexp*exp(-bexp*x) +cexp;
N = length(Valid);
%Sxy = sqrt(sum(((Bmat(Valid)-mean(Bmat(Valid)))).^2)./sum((Age(Valid)-mean(Age(Valid))).^2));
Sigma = sqrt(sum(((Bmat(Valid)-mean(Bmat(Valid)))).^2)./(N-1));
SSx = sum((Age(Valid)-mean(Age(Valid))).^2);
%SSx = sum(sqrt(1/N+(Age(Valid)-mean(Age(Valid))).^2));
SEM = Sigma * sqrt(1/N+(x-mean(Age(Valid))).^2./SSx);
ts = tinv([0.025  0.975],N-1);
yU = y + ts(1).*SEM;
yL = y + ts(2).*SEM;


hold on
% plot(x,y,'color',[0 0.5 0],'linewidth',2)
plot(x,yexp,'color',[160 82 45]/255,'linewidth',5)

% plot(x,yexp,'b','linewidth',2)
% plot(x,yU,'b','linewidth',1,'linestyle','--')
% plot(x,yL,'b','linewidth',1,'linestyle','--')
plot(ag,maM(ag),'ko','linewidth',3,'markersize',15,'markerfacecolor',[160 82 45]/255,'markeredgecolor',[160 82 45]/255)






set(gca,'fontsize',20)
xlabel('Age (yr)')

% hold on,
% plotSpread(Bmat([adultSetb]),'distributionmarkers','o','distributionColors','k','xvalues',xOffset)
% plotSpread(Bmat([adultSetg]),'distributionmarkers','^','distributionColors','k','xvalues',xOffset)
% set(gca,'xtick',[7:2:13 xOffset],'xticklabel',{'7','9','11','13','Adults'})
Nad = length(AD);
% SEM = std(AD);
SEM = std(AD)/sqrt(length(AD));               % Standard Error
ts = tinv([0.025  0.975],length(AD)-1);      % T-Score
CI = mean(AD) + ts*SEM;
% CI = prctile(AD,[2.5 97.5]);
hold on
% line([xOffset-.5 xOffset+.5],[CI(1) CI(1)],'linewidth',2,'color','b')
% line([xOffset-.5 xOffset+.5],[CI(2) CI(2)],'linewidth',2,'color','b')
line([xOffset xOffset],[CI(1) CI(2)],'linewidth',5,'color',[160 82 45]/255)
plot(xOffset,mean(AD),'ko','linewidth',3,'markersize',15,'markerfacecolor',[160 82 45]/255,'markeredgecolor',[160 82 45]/255)
set(gca,'fontsize',28)
xlabel('Age (yr)','fontweight','bold')
%title(sprintf('N_c_h=%d, N_a_d=%d, r=%1.2f, p=%0.2g',N,Nad,r,p))
%ylim([0 0.35])
%ylabel('Median Absolute Error (mm)','fontweight','bold')
ylabel('Absolute Error (mm)','fontweight','bold')

ylim([0 160])
xlim([5 22])


% adultSetb = [2 7:9 12:13 15:16 17 500];
% adultSetg = [1 3:4 6 18 19];
%
% AD = Bmat([adultSetb adultSetg])
% SEM = std(AD)/sqrt(length(AD));               % Standard Error
% ts = tinv([0.025  0.975],length(AD)-1);      % T-Score
% CI = mean(AD) + ts*SEM;


%% prepare data for stats

% for i = 1:length(Bmat)
%     stat = [Bmat([Valid adultSetb adultSetg]) Age([Valid adultSetb adultSetg])' [Valid adultSetb adultSetg]'];
% end

%%
for ag = 8:2:12,
    
    group{ag} = Valid(Age(Valid)>=ag-1 & Age(Valid)<ag+1);
    %     maM(ag) = mean(Bmat(maSubj));
    %     tsma = tinv([0.025  0.975],length(maSubj)-1);
    %     SEMma = std(Bmat(maSubj))./sqrt(length(maSubj));
    %     maCI(ag) = tsma(2)*SEMma;
    %     %[~,tT(ag)] = ttest2(AD,Bmat(maSubj));
    %     [tT2(ag), ~] = ranksum(AD,Bmat(maSubj));
    %     [~,tT(ag),~,STAT{ag}] = ttest2(Bmat(maSubj),AD);
    %     %[tT(ag), ~] = ranksum(AD,Bmat(maSubj));
    %
    %     L(ag) = length(Bmat(maSubj));
    %     cohenD(ag) = (maM(ag)-meanAD)./(sqrt(((L(ag)-1)*std(Bmat(maSubj)).^2 + (length(AD)-1)*sdAD.^2)/(L(ag)+length(AD)-2)));
    %
    %     hold on
    %     line([ag ag],[maM(ag)-maCI(ag) maM(ag)+maCI(ag)],'color','k','linewidth',3)
end

%%
% for isubj = 1:length(Valid)
%    if ismember(Valid(isubj),group{8})
%       stat(isubj,4) = 8;
%    elseif ismember(Valid(isubj),group{10})
%            stat(isubj,4) = 10;
%        elseif ismember(Valid(isubj),group{12})
%            stat(isubj,4) = 12;
% %    else
% %        stat(isubj,4) = 18;
%    end
%
% end



%%
stat = [];
nsubj = length(subjsel);
for isubj = 1:nsubj,
    subj = subjsel(isubj);
    Mselected = M{subj}(M{subj}~=0 & abs(M{subj})<mean(abs(M{subj}(M{subj}~=0)))+3*std(abs(M{subj}(M{subj}~=0))))*1000*1.1;
    %     Mselected = M{subj}(M{subj}~=0 & abs(M{subj})<0.30)*1000*1.1;
        Gselected = G{subj}(M{subj}~=0 & abs(M{subj})<mean(abs(M{subj}(M{subj}~=0)))+3*std(abs(M{subj}(M{subj}~=0))));

    nsel = length(Mselected);
    
    if ismember(subj,group{8})
        groupsel = 8;
    elseif ismember(subj,group{10})
        groupsel = 10;
    elseif ismember(subj,group{12})
        groupsel = 12;
    else
        groupsel = 0;
    end
    
    if isubj == 1
        nrow = 0;
    else
    nrow = size(stat,1);
    
    end
    stat([1:nsel]+nrow,:) = [abs(Mselected) Gselected ones(nsel,1)*Age(subj) ones(nsel,1)*subj ones(nsel,1)*groupsel];
    
    
end
% end

%%

csvwrite(sprintf('mat_catching_lab_data_for_stat.csv'),stat)


%%
figure, hist(age(subjsel),30)