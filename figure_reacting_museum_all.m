% This script is to create figure for Reacting results from individual data
% collected from the Museum of Science Boston.
% Written by Se-Woong Park
% MouseCatching_MoS_MT.mat and AgeGender449.mat should be in the same folder.

load MouseCatching_MoS_MT.mat

colorLab = [0 53 200]./255;
colorMuseum = [0 171 240]./255;
lblue = colorMuseum;
lred = colorMuseum;

load AgeGender449.mat
figure
Male = find(Age>=5 & Age<130 & Gender=='M' & Hist_Mov=='N' & Hist_Psych=='N')';
Female = find(Age>=5 & Age<130 & Gender=='F' & Hist_Mov=='N' & Hist_Psych=='N')';

for subj = [Male Female],
   
    B{subj} = median(R{subj}(R{subj}>0.15 & R{subj}<0.5),'omitnan');


    Bmat(subj) = B{subj};

end

hold on
%patch([5 10 10 5],[0 0 400 400],[1 1 1]*0.9,'edgecolor','none')
%patch([10 20 20 10],[0 0 160 160],[1 1 1]*0.8)
%patch([20 30 30 20],[0 0 400 400],[1 1 1]*0.9,'edgecolor','none')

W=1;
oldest = 60;
ERROR = Bmat'*1000;
%& mERROR>0
mERROR = ERROR;
p1 = plot(Age(Male),mERROR(Male),'o','color',lblue,'linewidth',2,'markersize',10,'markerfacecolor',lblue)
hold on
p2 = plot(Age(Female),mERROR(Female),'o','color',lred,'linewidth',2,'markersize',10,'markerfacecolor',lred)
% plot([5:10:55]-0.2,AGEM(5:10:55),'-','color',[0    0.4470    0.7410],'linewidth',3)
% plot([5:10:55]+0.2,AGEF(5:10:55),'-','color',[0.8500    0.3250    0.0980],'linewidth',3)

AD = mERROR(Age>=18 & Age<22 & mERROR>0);

xax = [7.5 15 25 60];
x = 0;
W = 5;
hold on


for a = 5:10:35,
    x = x + 1;
%     MALE = find(Age>a-W & Age<a+W & Gender=='M' & History_Mov=='N' & History_Psych=='N')';
%     FEMALE = find(Age>a-W & Age<a+W & Gender=='F' & History_Mov=='N' & History_Psych=='N')';
    if a<35,
        SUBJ = find(Age>=a-W & Age<a+W & mERROR>0 & Hist_Mov=='N' & Hist_Psych=='N');
        MALE = find(Gender == 'M' & Age>=a-W & Age<a+W & mERROR>0 & Hist_Mov=='N' & Hist_Psych=='N');
        FEMALE = find(Gender == 'F' & Age>=a-W & Age<a+W & mERROR>0 & Hist_Mov=='N' & Hist_Psych=='N');
        size(MALE)
        size(FEMALE)
    else
        SUBJ = find(Age>=a-W & mERROR>0 & Hist_Mov=='N' & Hist_Psych=='N');
        MALE = find(Gender == 'M' & Age>=a-W & mERROR>0 & Hist_Mov=='N' & Hist_Psych=='N');
        FEMALE = find(Gender == 'F' & Age>=a-W & mERROR>0 & Hist_Mov=='N' & Hist_Psych=='N');
        size(MALE)
        size(FEMALE)
    end
    %AGE(a) = mean(mean(ERROR([MALE FEMALE],:),2),'omitnan');
%     AGEl(a) = AGE(a) - std(ERROR([MALE FEMALE]),'omitnan')./sqrt(length([MALE FEMALE]));
%     AGEu(a) = AGE(a) + std(ERROR([MALE FEMALE]),'omitnan')./sqrt(length([MALE FEMALE]));
    AGE(a) = mean(mERROR(SUBJ),'omitnan');
    AGEM(a) = mean(mERROR(MALE),'omitnan');
    AGEF(a) = mean(mERROR(FEMALE),'omitnan');
tsmaS = tinv([0.025  0.975],length([SUBJ])-1);  
tsmaM = tinv([0.025  0.975],length([MALE])-1);  
tsmaF = tinv([0.025  0.975],length([FEMALE])-1);  

    tT2(a) = ranksum(mERROR(SUBJ),AD);
    AGEl(a) = AGE(a) - std(mERROR([SUBJ]),'omitnan')./sqrt(length([SUBJ])).*tsmaS(2);
    AGEu(a) = AGE(a) + std(mERROR([SUBJ]),'omitnan')./sqrt(length([SUBJ])).*tsmaS(2);
    AGElM(a) = AGEM(a) - std(mERROR([MALE]),'omitnan')./sqrt(length([MALE])).*tsmaM(2);
    AGEuM(a) = AGEM(a) + std(mERROR([MALE]),'omitnan')./sqrt(length([MALE])).*tsmaM(2);
    AGElF(a) = AGEF(a) - std(mERROR([FEMALE]),'omitnan')./sqrt(length([FEMALE])).*tsmaF(2);
    AGEuF(a) = AGEF(a) + std(mERROR([FEMALE]),'omitnan')./sqrt(length([FEMALE])).*tsmaF(2);
    line([xax(x) xax(x)],[AGEl(a) AGEu(a)],'color',[0 0 0],'linewidth',3)
%     line([xax(x)-0.4 xax(x)-0.4],[AGElM(a) AGEuM(a)],'color',[0 0 1],'linewidth',3)
%     line([xax(x)+0.4 xax(x)+0.4],[AGElF(a) AGEuF(a)],'color',[1 0 0],'linewidth',3)

end
p3 = plot([7.5 15 25 60],AGE(5:10:35),'-o','color',[0 0 0],'linewidth',3,'markersize',7,'markerfacecolor',[0 0 0])

% p3 = plot([7.5 15 25 60]-0.4,AGEM(5:10:35),'-','color',[0 0 1],'linewidth',3)
% p4 = plot([7.5 15 25 60]+0.4,AGEF(5:10:35),'-','color',[1 0 0],'linewidth',3)
% legend([p1,p2,p3,p4],'Male','Female','Male','Female')

% 
% stepAge = 20;
% ageSubj = [Age,[1:449]'];
% sortAge = sortrows(ageSubj(mERROR>0 & Age>0 & isfinite(Age),:),1,'descend');
% seq = 0;
% for a = 1:stepAge:length(sortAge)-stepAge,
%     seq = seq + 1;
%     mAge(seq) = mean(sortAge(a:a+stepAge-1,1),'omitnan');
%     mRT(seq) = mean(mERROR(sortAge(a:a+stepAge-1,2)),'omitnan');
% end
% hold on
% plot(mAge,mRT,'r-','linewidth',3)
set(gca,'fontsize',24)
xlabel('Age (yr)','fontweight','bold')
ylim([150 400])
% 
% figure
% plot(Age([Male Female]),mean(Tree([Male Female],:),2),'o','color',[0.9 0.447 0.041],'linewidth',2)
% 
% set(gca,'fontsize',20)
% xlabel('Age (yr)')
