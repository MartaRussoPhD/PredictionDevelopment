% This script is to create figure for Catching results from children data
% collected from the Museum of Science Boston.
% Written by Se-Woong Park
% catching_MoS_new.mat and AgeGender449.mat should be in the same folder.
load catching_MoS_new G M
load AgeGender449.mat
%load Survey_318 Age Gender Hist*
% lblue = [102 178 255]/255;
% lred = [255 153 153]/255;
% lblue = [1 1 1]*0.7;
% lred = [1 1 1]*0.7;


colorLab = [0 53 200]./255;
colorMuseum = [0 171 240]./255;
lblue = colorMuseum;
lred = colorMuseum;
Male = find(Age>0 & Age<130 & Gender=='M' & Hist_Mov=='N' & Hist_Psych=='N')';
Female = find(Age>0 & Age<130 & Gender=='F' & Hist_Mov=='N' & Hist_Psych=='N')';
% subjsel = [Male Female]; %MR20220721

figure

for subj = [Male Female],
    
    Mselected = M{subj}(M{subj}~=0 & abs(M{subj})<mean(abs(M{subj}(M{subj}~=0)),'omitnan')+3*std(abs(M{subj}(M{subj}~=0)),'omitnan'));


    if length(Mselected)>=3,

        Tree(subj,1) = median(abs(Mselected),'omitnan');
    else
        Tree(subj,1) = NaN;
    end
    
end
% figure
% plot(Age(Male),Tree(Male,tree),'o','color',[0 0.447 0.741],'linewidth',2)
% 
% set(gca,'fontsize',20)
% xlabel('Age (yr)')
% 
% figure
% plot(Age(Female),Tree(Female,tree),'o','color',[0.9 0.447 0.041],'linewidth',2)
% 
% set(gca,'fontsize',20)
% xlabel('Age (yr)')

mERROR = Tree * 1000*1.1;
plot(Age(Male),mERROR(Male),'o','markerfacecolor',lblue,'markeredgecolor','none','linewidth',2,'markersize',12)
hold on
plot(Age(Female),mERROR(Female),'o','markerfacecolor',lred,'markeredgecolor','none','linewidth',2,'markersize',12)

AD = mERROR(Age>=18 & Age<22 & mERROR>0);
meanAD = mean(AD);
sdAD = std(AD);
hold on
CORR = find(Age>=5 & Age<13 & mERROR>0);
[r p] = corr(Age(CORR), mERROR(CORR))
fcn = 'a*x +b';
fcnexp = 'a*exp(-b*x) +c';
[F GoF] = fit(Age(CORR),mERROR(CORR),fcn,'start',[-1 1]);
[Fexp GoFexp] = fit(Age(CORR),mERROR(CORR),fcnexp,'start',[150 1/100 50]);

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


x = 5:0.1:13;
y = a*x +b;
yexp = aexp*exp(-bexp*x) +cexp;
%plot(x,y,'color',[0 0.5 0],'linewidth',2)
plot(x,yexp,'color',[160 85 42]/255,'linewidth',5)

% plot(Age(Male),mERROR(Male),'o','markerfacecolor','none','markeredgecolor',[1 1 1]*0.4,'linewidth',2,'markersize',10)
% hold on
% plot(Age(Female),mERROR(Female),'^','markerfacecolor','none','markeredgecolor',[1 1 1]*0.4,'linewidth',2,'markersize',10)


% plot(Age(Male),mERROR(Male),'o','markerfacecolor','none','markeredgecolor',[1 1 1]*0.5,'linewidth',2,'markersize',10)
% hold on
% plot(Age(Female),mERROR(Female),'o','markerfacecolor','none','markeredgecolor',[1 1 1]*0.5,'linewidth',2,'markersize',10)

% stepAge = 20;
% ageSubj = [Age,[1:449]'];
% sortAge = sortrows(ageSubj(Age>0 & isfinite(Age),:),1,'descend');
% seq = 0;
% for a = 1:stepAge:length(sortAge)-stepAge,
%     seq = seq + 1;
%     mAge(seq) = mean(sortAge(a:a+stepAge,1),'omitnan');
%     mRT(seq) = mean(mERROR(sortAge(a:a+stepAge,2)),'omitnan');
% end
% hold on
% plot(mAge,mRT,'r-','linewidth',3)


W = 1;
subjsel = [ ]; %MR20220721
for a = [6:2:16],
%     MALE = find(Age>a-W & Age<a+W & Gender=='M' & History_Mov=='N' & History_Psych=='N')';
%     FEMALE = find(Age>a-W & Age<a+W & Gender=='F' & History_Mov=='N' & History_Psych=='N')';
    SUBJ = find(Age>=a-W & Age<a+W & mERROR>0);
    %AGE(a) = mean(mean(ERROR([MALE FEMALE],:),2),'omitnan');
%     AGEl(a) = AGE(a) - std(ERROR([MALE FEMALE]),'omitnan')./sqrt(length([MALE FEMALE]));
%     AGEu(a) = AGE(a) + std(ERROR([MALE FEMALE]),'omitnan')./sqrt(length([MALE FEMALE]));
    AGE(a) = mean(mERROR(SUBJ),'omitnan');
    tsma = tinv([0.025  0.975],length([SUBJ])-1);  
    [~,tT(a),~,STAT{a}] = ttest2(mERROR(SUBJ),AD);
    tT2(a) = ranksum(mERROR(SUBJ),AD);
    L(a) = length(mERROR(SUBJ));
    cohenD(a) = (AGE(a)-meanAD)./(sqrt(((L(a)-1)*std(mERROR(SUBJ)).^2 + (length(AD)-1)*sdAD.^2)/(L(a)+length(AD)-2)));

    AGEl(a) = AGE(a) - std(mERROR([SUBJ]),'omitnan')./sqrt(L(a)).*tsma(2);
    AGEu(a) = AGE(a) + std(mERROR([SUBJ]),'omitnan')./sqrt(L(a)).*tsma(2);
        group{a} = SUBJ; %MR20220802
subjsel = [ subjsel SUBJ']; %MR20220721
    line([a a],[AGEu(a) AGEl(a)],'color','k','linewidth',3,'markersize',7,'markerfacecolor',[0 0 0])
    
end
% subjsel = [SUBJ]; %MR20220721

W = 2;
for a = [20],
%     MALE = find(Age>a-W & Age<a+W & Gender=='M' & History_Mov=='N' & History_Psych=='N')';
%     FEMALE = find(Age>a-W & Age<a+W & Gender=='F' & History_Mov=='N' & History_Psych=='N')';
    SUBJ = find(Age>=a-W & Age<a+W & mERROR>0);
    %AGE(a) = mean(mean(ERROR([MALE FEMALE],:),2),'omitnan');
%     AGEl(a) = AGE(a) - std(ERROR([MALE FEMALE]),'omitnan')./sqrt(length([MALE FEMALE]));
%     AGEu(a) = AGE(a) + std(ERROR([MALE FEMALE]),'omitnan')./sqrt(length([MALE FEMALE]));
    AGE(a) = mean(mERROR(SUBJ),'omitnan');
    tsma = tinv([0.025  0.975],length([SUBJ])-1);  
    tT2(a) = ranksum(mERROR(SUBJ),AD);
    L(a) = length(mERROR(SUBJ));
    AGEl(a) = AGE(a) - std(mERROR([SUBJ]),'omitnan')./sqrt(L(a)).*tsma(2);
    AGEu(a) = AGE(a) + std(mERROR([SUBJ]),'omitnan')./sqrt(L(a)).*tsma(2);
        group{a} = SUBJ; %MR20220802

    line([a a],[AGEu(a) AGEl(a)],'color','k','linewidth',3,'markersize',10)
end
subjsel = [ subjsel SUBJ']; %MR20220721

plot([6:2:13 20],AGE([6:2:13 20]),'ok','linewidth',3,'markersize',7,'markerfacecolor',[0 0 0])


set(gca,'fontsize',24,'box','off')
xlabel('Age (yr)','fontweight','bold')
xlim([5 22])
ylim([0 160])
% 
% figure
% plot(Age([Male Female]),mean(Tree([Male Female],:),2),'o','color',[0.9 0.447 0.041],'linewidth',2)
% 
% set(gca,'fontsize',20)
% xlabel('Age (yr)')
ylabel('Median Absolute Error (mm)','fontweight','bold')


%% MR20220802
% for ag = 6:2:12,
%     
%     group{ag} = subjsel(Age(subjsel)>=ag-1 & Age(subjsel)<ag+1);
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
    
    if ismember(subj,group{6})
        groupsel = 6;
    elseif ismember(subj,group{8})
        groupsel = 8;
    elseif ismember(subj,group{10})
        groupsel = 10;
    elseif ismember(subj,group{12})
        groupsel = 12;
         elseif ismember(subj,group{14})
        groupsel = 14;
         elseif ismember(subj,group{16})
        groupsel = 16;
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


%%
csvwrite(sprintf('mat_catching_museum_children_data_for_stat.csv'),stat)
