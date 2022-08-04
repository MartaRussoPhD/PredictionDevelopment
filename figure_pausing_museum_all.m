% This script is to create figure for Pausing results from individual data
% collected from the Museum of Science Boston.
% Written by Se-Woong Park
% window_MoS.mat and AgeGender449.mat should be in the same folder.
colorMuseum = [0 171 240]./255;
lblue = colorMuseum;
lred = colorMuseum;

load window_MoS G M
load AgeGender449.mat

Male = find(Age>0 & Age<130 & Gender=='M' & Hist_Mov=='N' & Hist_Psych=='N')';
Female = find(Age>0 & Age<130 & Gender=='F' & Hist_Mov=='N' & Hist_Psych=='N')';
figure

for tree = 4,


for subj = [Male Female],
    
    Limit = mean(abs(M{subj}(M{subj}~=0)),'omitnan')+3*std(abs(M{subj}(M{subj}~=0)),'omitnan');
%     Limit = 0.25;
    Mselected = M{subj}(M{subj}~=0 & abs(M{subj})<Limit);
%     Mselected = M{subj}(M{subj}~=0 & abs(M{subj})<0.3);

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

end
ERROR = Tree * 1000;

mERROR = ERROR;
mERROR(mERROR>140| mERROR==0)=NaN;

hold on
%patch([5 10 10 5],[0 0 160 160],[1 1 1]*0.9,'edgecolor','none')
%patch([10 20 20 10],[0 0 160 160],[1 1 1]*0.8)
%patch([20 30 30 20],[0 0 160 160],[1 1 1]*0.9,'edgecolor','none')



p1 = plot(Age(Male),mERROR(Male),'o','color',lblue,'linewidth',2,'markersize',10,'markerfacecolor',lblue);
hold on
p2 = plot(Age(Female),mERROR(Female),'o','color',lred,'linewidth',2,'markersize',10,'markerfacecolor',lred);
% plot([5:10:55]-0.2,AGEM(5:10:55),'-','color',[0    0.4470    0.7410],'linewidth',3)
% plot([5:10:55]+0.2,AGEF(5:10:55),'-','color',[0.8500    0.3250    0.0980],'linewidth',3)

AD = mERROR(Age>=18 & Age<22 & mERROR>0);

W = 5;
x = 0;
xax = [7.5 15 25 60];
hold on
subjsel = [ ]; %MR20220721

for a = 5:10:35,
    x = x + 1;
%     MALE = find(Age>a-W & Age<a+W & Gender=='M' & History_Mov=='N' & History_Psych=='N')';
%     FEMALE = find(Age>a-W & Age<a+W & Gender=='F' & History_Mov=='N' & History_Psych=='N')';
    if a>5 & a<35,
        SUBJ = find(Age>=a-W & Age<a+W & mERROR>0 & Hist_Mov=='N' & Hist_Psych=='N');
        MALE = find(Gender == 'M' & Age>=a-W & Age<a+W & mERROR>0 & Hist_Mov=='N' & Hist_Psych=='N');
        FEMALE = find(Gender == 'F' & Age>=a-W & Age<a+W & mERROR>0 & Hist_Mov=='N' & Hist_Psych=='N');
        size(MALE)
        size(FEMALE)
    elseif a==5,
        SUBJ = find(Age>=a & Age<a+W & mERROR>0 & Hist_Mov=='N' & Hist_Psych=='N');
        MALE = find(Gender == 'M' & Age>=a & Age<a+W & mERROR>0 & Hist_Mov=='N' & Hist_Psych=='N');
        FEMALE = find(Gender == 'F' & Age>=a & Age<a+W & mERROR>0 & Hist_Mov=='N' & Hist_Psych=='N');
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
%     line([xax(x)-0.4 xax(x)-0.4],[AGElM(a) AGEuM(a)],'color',[0 0 1],'linewidth',3)
%     line([xax(x)+0.4 xax(x)+0.4],[AGElF(a) AGEuF(a)],'color',[1 0 0],'linewidth',3)
    line([xax(x) xax(x)],[AGEl(a) AGEu(a)],'color',[0 0 0],'linewidth',3)
    S{x} = SUBJ;
    
      group{a} = SUBJ; %MR20220802

subjsel = [ subjsel SUBJ']; %MR20220721
end
p3 = plot([7.5 15 25 60],AGE(5:10:35),'-o','color',[0 0 0],'linewidth',3,'markersize',7,'markerfacecolor',[0 0 0]);
%CAT = [S{1};S{2};S{3};S{4}];
%GROUP = [ones(length(S{1}),1); 2*ones(length(S{2}),1); 3*ones(length(S{3}),1); 4*ones(length(S{4}),1)];
%[P,ANOVATAB,STATS] = kruskalwallis(mERROR(CAT),GROUP);


% p3 = plot([7.5 15 25 60]-0.4,AGEM(5:10:35),'-','color',[0 0 1],'linewidth',3);
% p4 = plot([7.5 15 25 60]+0.4,AGEF(5:10:35),'-','color',[1 0 0],'linewidth',3);
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
ylim([0 140])
set(gca,'fontsize',24)
xlabel('Age (yr)','fontweight','bold')
% 
% figure
% plot(Age([Male Female]),mean(Tree([Male Female],:),2),'o','color',[0.9 0.447 0.041],'linewidth',2)
% 
% set(gca,'fontsize',20)
% xlabel('Age (yr)')

%% MR 20220802

% for ag = 5:10:35,
%     
%     group{ag} = subjsel(Age(subjsel)>=ag-W & Age(subjsel)<ag+W);
%     
% end


%%
stat = [];
nsubj = length(subjsel);
for isubj = 1:nsubj,
    subj = subjsel(isubj);
   Limit = mean(abs(M{subj}(M{subj}~=0)),'omitnan')+3*std(abs(M{subj}(M{subj}~=0)),'omitnan');
%     Limit = 0.25;
    Mselected = M{subj}(M{subj}~=0 & abs(M{subj})<Limit);
%     Mselected = M{subj}(M{subj}~=0 & abs(M{subj})<0.3);

    nsel = length(Mselected);
    
    if ismember(subj,group{5})
        groupsel = 5;
    elseif ismember(subj,group{15})
        groupsel = 15;
    elseif ismember(subj,group{25})
        groupsel = 25;
        elseif ismember(subj,group{35})
        groupsel = 35;
    else
        groupsel = 0;
    end
    
    if isubj == 1
        nrow = 0;
    else
    nrow = size(stat,1);
    
    end
    stat([1:nsel]+nrow,:) = [abs(Mselected) ones(nsel,1)*Age(subj) ones(nsel,1)*subj ones(nsel,1)*groupsel];
    
    
end
% end

%%

csvwrite(sprintf('mat_pausing_museum_all_data_for_stat.csv'),stat)
