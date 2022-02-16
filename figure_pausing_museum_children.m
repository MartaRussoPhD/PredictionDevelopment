% This script is to create figure for Pausing results from children data
% collected from the Museum of Science Boston.
% Written by Se-Woong Park
% window_MoS.mat and AgeGender449.mat should be in the same folder.
load window_MoS

colorMuseum = [0 171 240]./255;
lblue = colorMuseum;
lred = colorMuseum;



load AgeGender449.mat
Male = find(Age>0 & Age<130 & Gender=='M' & Hist_Mov=='N' & Hist_Psych=='N')';
Female = find(Age>0 & Age<130 & Gender=='F' & Hist_Mov=='N' & Hist_Psych=='N')';

figure


for subj = [Male Female],

    Limit = mean(abs(M{subj}(M{subj}~=0)),'omitnan')+3*std(abs(M{subj}(M{subj}~=0)),'omitnan');
%     Limit = 0.25;
    Mselected = M{subj}(M{subj}~=0 & abs(M{subj})<Limit)*1000;

    if length(Mselected)>=3,
        ERROR(subj,1) = median(abs(Mselected),'omitnan');

    else
        ERROR(subj,1) = NaN;
    end

end




mERROR = ERROR;
mERROR(mERROR>140| mERROR==0)=NaN;

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

% plot(Age(Male),mERROR(Male),'o','markerfacecolor',lblue,'markeredgecolor','none','linewidth',2,'markersize',10)
% hold on
% plot(Age(Female),mERROR(Female),'o','markerfacecolor',lred,'markeredgecolor','none','linewidth',2,'markersize',10)
plot(Age(Female),mERROR(Female),'o','linewidth',2,'markersize',10,'markerfacecolor',lred,'markeredgecolor',lred)
hold on
plot(Age(Male),mERROR(Male),'o','linewidth',2,'markersize',10,'markerfacecolor',lblue,'markeredgecolor',lblue)

x = 5:0.1:13;
y = a*x +b;
yexp = aexp*exp(-bexp*x) +cexp;
% plot(x,yexp,'color',[0 0.5 0],'linewidth',3)
plot(x,yexp,'color',[160,82,45]/255,'linewidth',5)
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
    
    line([a a],[AGEu(a) AGEl(a)],'color','k','linewidth',3)
end
W = 2;
for a = 20,
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
    
    line([a a],[AGEu(a) AGEl(a)],'color','k','linewidth',3)
end

plot([6:2:12 20],AGE([6:2:12 20]),'ok','linewidth',3,'markersize',7,'markerfacecolor',[0 0 0])


set(gca,'fontsize',24)
xlabel('Age (yr)','fontweight','bold')
% ylabel('Absolute Error (ms)')
xlim([5 22])
ylim([0 140])
% 
% figure
% plot(Age([Male Female]),mean(Tree([Male Female],:),2),'o','color',[0.9 0.447 0.041],'linewidth',2)
% 
% set(gca,'fontsize',20)
% xlabel('Age (yr)')

