% This script is to create figure for Reacting results from children data
% collected from the Museum of Science Boston.
% Written by Se-Woong Park
% MouseCatching_MoS_MT.mat and AgeGender449.mat should be in the same folder.
load MouseCatching_MoS_MT.mat

load AgeGender449.mat
% lblue = [102 178 255]/255;
% lred = [255 153 153]/255;
colorLab = [0 53 200]./255;
colorMuseum = [0 171 240]./255;
lblue = colorMuseum;
lred = colorMuseum;
Male = find(Age>0 & Gender=='M' & (Hist_Mov=='N' & Hist_Psych=='N'))';

Female = find(Age>0 & Gender=='F'& (Hist_Mov=='N' & Hist_Psych=='N'))';

figure

tree = 4;

for subj = [Male Female],
    
    B{subj} = median(R{subj}(R{subj}>0.15 & R{subj}<0.5),'omitnan');

    hold on
    Bmat(subj) = B{subj};
end


W=1;
oldest = 60;
ERROR = Bmat' * 1000;
%& mERROR>0
mERROR = ERROR;

Male = intersect(Male,find(mERROR<400));
Female = intersect(Female,find(mERROR<400));
plot(Age(Male),mERROR(Male),'o','color',lblue,'markerfacecolor',lblue,'linewidth',2,'markersize',12)
hold on
plot(Age(Female),mERROR(Female),'o','color',lred,'markerfacecolor',lred,'linewidth',2,'markersize',12)



stepAge = 20;
ageSubj = [Age,[1:449]'];
sortAge = sortrows(ageSubj(mERROR>0 & Age>0 & isfinite(Age),:),1,'descend');
seq = 0;
binStart = 1:stepAge:length(sortAge)-stepAge;
for a = binStart,
    seq = seq + 1;
    if a == binStart(end),
        mAge(seq) = mean(sortAge(a:end,1),'omitnan');
        mRT(seq) = mean(mERROR(sortAge(a:end,2)),'omitnan');
    else
        mAge(seq) = mean(sortAge(a:a+stepAge-1,1),'omitnan');
        mRT(seq) = mean(mERROR(sortAge(a:a+stepAge-1,2)),'omitnan');
    end
end

% plot(mAge,mRT,'r-','linewidth',3)
CORR = find(Age>=5 & Age<13 & mERROR>0);
[r p] = corr(Age(CORR), mERROR(CORR))
fcn = 'a*x +b';
fcnexp = 'a*exp(-b*x) +c';
[F GoF] = fit(Age(CORR),mERROR(CORR),fcn,'start',[-1 1]);
[Fexp GoFexp] = fit(Age(CORR),mERROR(CORR),fcnexp,'start',[150 1/100 250]);

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
% plot(x,yexp,'color',[0 0.5 0],'linewidth',2)
plot(x,yexp,'color',[160 82 45]/255,'linewidth',5)

AD = mERROR(Age>=18 & Age<22 & mERROR>0);
meanAD = mean(AD);
sdAD = std(AD);
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

    AGEl(a) = AGE(a) - std(mERROR([SUBJ]),'omitnan')./sqrt(length([SUBJ])).*tsma(2);
    AGEu(a) = AGE(a) + std(mERROR([SUBJ]),'omitnan')./sqrt(length([SUBJ])).*tsma(2);
    line([a a],[AGEu(a) AGEl(a)],'color','k','linewidth',3,'markersize',10)

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
    AGEl(a) = AGE(a) - std(mERROR([SUBJ]),'omitnan')./sqrt(length([SUBJ])).*tsma(2);
    AGEu(a) = AGE(a) + std(mERROR([SUBJ]),'omitnan')./sqrt(length([SUBJ])).*tsma(2);
    line([a a],[AGEu(a) AGEl(a)],'color','k','linewidth',3)

end
plot([6:2:13 20],AGE([6:2:13 20]),'ok','linewidth',3,'markersize',7,'markerfacecolor',[0 0 0])

% hold on
% plot([5+W:oldest],AGE(5+W:end),'-r','linewidth',3)
% % plot([5+W:60],AGEl(5+W:end),'-r','linewidth',1)
% % plot([5+W:60],AGEu(5+W:end),'-r','linewidth',1)
set(gca,'fontsize',24)
xlabel('Age (yr)','fontweight','bold')
ylim([150 400])

xlim([5 22])
