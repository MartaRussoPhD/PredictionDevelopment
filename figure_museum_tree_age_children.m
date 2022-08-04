%load GM_bouncing_latest.mat
%load B_bouncing B
% load GM_tree_catching_error_all M G
% load GM_mouse_catching M
%load GM_bouncing M
% load ball_bouncing E
% load GM_tree_catching_RT TC
% load GM_3_mouse G M
% load GM_window_button.mat
% load GM_mouse_catching_MT2.mat
% MT = M;
% load GM_mouse_catching.mat
% load MouseCatching_MoS.mat
load catching_MoS_new G M
% M = TC;
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
% Male = find(Age>0 & Age<130 & Gender=='M')';
Female = find(Age>0 & Age<130 & Gender=='F' & Hist_Mov=='N' & Hist_Psych=='N')';
% Female = find(Age>0 & Age<130 & Gender=='F')';
%M = E;
figure

for subj = [Male Female],
    %B{subj} = abs(median(M{subj}(M{subj}>0.15),'omitnan'));
    %B{subj} = median(M{subj}(M{subj}>0.15),'omitnan');
    %B{subj} = median(abs(M{subj}(G{subj}==3 & M{subj}>0)),'omitnan');
    %B{subj} = median(abs(M{subj}(G{subj}'==tree & M{subj}~=0)),'omitnan');
    %B{subj} = median(abs(M{subj}(G{subj}==tree & M{subj}~=0)),'omitnan');

    %B{subj} = std(M{subj}(M{subj}~=0),'omitnan');
    %B{subj} = size(M{subj},2)/60;
%     B{subj} = mean(M{subj});
%     B{subj} = median(MT{subj}(M{subj}>0.15 & M{subj}<0.5),'omitnan');
%     B{subj} = median(M{subj}(M{subj}>0.15 & M{subj}<0.5),'omitnan');
    %Mselected = P{subj}(G{subj}==tree & M{subj}~=0 & abs(M{subj})<0.15)*1000*1.1;
%     Mselected = M{subj}(G{subj}==tree & M{subj}~=0 & abs(M{subj})<mean(abs(M{subj}(M{subj}~=0)),'omitnan')+3*std(abs(M{subj}(M{subj}~=0)),'omitnan'));
    Mselected = M{subj}(M{subj}~=0 & abs(M{subj})<mean(abs(M{subj}(M{subj}~=0)),'omitnan')+3*std(abs(M{subj}(M{subj}~=0)),'omitnan'));
%     Mselected = M{subj}(M{subj}~=0 & abs(M{subj})<0.3);
%     Mselected = M{subj}(M{subj}~=0 & abs(M{subj})<0.4);

    if length(Mselected)>=3,
%     if length(Mselected(isfinite(Mselected)))>=1,

%         Tree(subj,tree) = median(abs(Mselected),'omitnan');
        Tree(subj,1) = median(abs(Mselected),'omitnan');
    else
%         Tree(subj,tree) = NaN;
        Tree(subj,1) = NaN;
    end
    
        
%     B{subj} = median(abs(M{subj}(G{subj}<tree & M{subj}~=0)),'omitnan');
    %B{subj} = median(abs(M{subj}(G{subj}<tree & M{subj}~=0)),'omitnan');

    %B{subj} = std(M{subj}(G{subj}<tree & M{subj}~=0),'omitnan');
    %plot(age(subj),abs(median(M{subj},'omitnan')),'o')
%     plot(Age(subj),B{subj},'o','color',[0.9 0.447 0.041],'linewidth',2)
%     hold on
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
for a = [6:2:12],
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
    
    line([a a],[AGEu(a) AGEl(a)],'color','k','linewidth',3,'markersize',7,'markerfacecolor',[0 0 0])
end

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
    
    line([a a],[AGEu(a) AGEl(a)],'color','k','linewidth',3,'markersize',10)
end

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
%ylabel('Mean Absolute Error (mm)','fontweight','bold')
