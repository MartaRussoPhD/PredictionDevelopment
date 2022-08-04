% This script is to create figure for Choosing results from inidvidual data
% collected from the Museum of Science Boston.
% Written by Se-Woong Park
% 3mouse_button_MoS.mat and AgeGender449.mat should be in the same folder.
colorLab = [0 53 200]./255;
colorMuseum = [0 171 240]./255;
lblue = colorMuseum;
lred = colorMuseum;
load 3mouse_button_MoS


load AgeGender449.mat
W = 5;
x = 0;
xax = [7.5 15 25 60];

Male = find(Age>0 & Age<130 & Gender=='M' & Hist_Mov=='N' & Hist_Psych=='N')';
Female = find(Age>0 & Age<130 & Gender=='F' & Hist_Mov=='N' & Hist_Psych=='N')';
%M = E;
figure

% for speed = 1:3,

for subj = [Male Female],
    
    Mselected = M{subj};
    if length(Mselected)>=3,

        ERROR(subj,1) = mean(Mselected)*100;


    else
        ERROR(subj,1) = NaN;
    end

end
mERROR = ERROR;
hold on


p1 = plot(Age(Male),mERROR(Male),'o','color',lblue,'markerfacecolor',lblue,'linewidth',2,'markersize',10)
hold on
p2 = plot(Age(Female),mERROR(Female),'o','color',lred,'markerfacecolor',lred,'linewidth',2,'markersize',10)
% plot([5:10:55]-0.2,AGEM(5:10:55),'-','color',[0    0.4470    0.7410],'linewidth',3)
% plot([5:10:55]+0.2,AGEF(5:10:55),'-','color',[0.8500    0.3250    0.0980],'linewidth',3)

AD = mERROR(Age>=18 & Age<22 & mERROR>0);

hold on
for a = 5:10:35,
        x = x + 1;

%     MALE = find(Age>a-W & Age<a+W & Gender=='M' & History_Mov=='N' & History_Psych=='N')';
%     FEMALE = find(Age>a-W & Age<a+W & Gender=='F' & History_Mov=='N' & History_Psych=='N')';
    if a<55,
        SUBJ = find(Age>=a-W & Age<a+W & mERROR>0 & Hist_Mov=='N' & Hist_Psych=='N');
        MALE = find(Gender == 'M' & Age>=a-W & Age<a+W & mERROR>0 & Hist_Mov=='N' & Hist_Psych=='N');
        FEMALE = find(Gender == 'F' & Age>=a-W & Age<a+W & mERROR>0 & Hist_Mov=='N' & Hist_Psych=='N');
    else
        SUBJ = find(Age>=a-W & mERROR>0 & Hist_Mov=='N' & Hist_Psych=='N');
        MALE = find(Gender == 'M' & Age>=a-W & mERROR>0 & Hist_Mov=='N' & Hist_Psych=='N');
        FEMALE = find(Gender == 'F' & Age>=a-W & mERROR>0 & Hist_Mov=='N' & Hist_Psych=='N');
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
    [r p] = ttest2(mERROR(MALE(isfinite(MALE))),mERROR(FEMALE(isfinite(FEMALE))))
    tT2(a) = ranksum(mERROR(SUBJ),AD);
    AGEl(a) = AGE(a) - std(mERROR([SUBJ]),'omitnan')./sqrt(length([SUBJ])).*tsmaS(2);
    AGEu(a) = AGE(a) + std(mERROR([SUBJ]),'omitnan')./sqrt(length([SUBJ])).*tsmaS(2);
    AGElM(a) = AGEM(a) - std(mERROR([MALE]),'omitnan')./sqrt(length([MALE])).*tsmaM(2);
    AGEuM(a) = AGEM(a) + std(mERROR([MALE]),'omitnan')./sqrt(length([MALE])).*tsmaM(2);
    AGElF(a) = AGEF(a) - std(mERROR([FEMALE]),'omitnan')./sqrt(length([FEMALE])).*tsmaF(2);
    AGEuF(a) = AGEF(a) + std(mERROR([FEMALE]),'omitnan')./sqrt(length([FEMALE])).*tsmaF(2);
%     line([xax(x)-0.4 xax(x)-0.4],[AGElM(a) AGEuM(a)],'color',[0 0 1],'linewidth',3)
%     line([xax(x)+0.4 xax(x)+0.4],[AGElF(a) AGEuF(a)],'color',[1 0 0],'linewidth',3)

    line([xax(x) xax(x)],[AGEl(a) AGEu(a)],'color',[1 1 1]*0,'linewidth',3)
        group{a} = SUBJ; %MR20220802

    
end
p3 = plot([7.5 15 25 60],AGE(5:10:35),'-o','color',[1 1 1]*0,'linewidth',3,'markersize',7,'markerfacecolor',[0 0 0]);

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
set(gca,'fontsize',24,'box','off')
xlabel('Age (yr)','fontweight','bold')
xlim([0 100])

ylim([0 120])
set(gca,'ytick',0:20:100)

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
    if not(isempty(M{subj}))
    Mselected = M{subj};%(M{subj}~=0 & abs(M{subj})<mean(abs(M{subj}(M{subj}~=0)))+3*std(abs(M{subj}(M{subj}~=0))))*1000*1.1;
    %     Mselected = M{subj}(M{subj}~=0 & abs(M{subj})<0.30)*1000*1.1;
        Gselected = G{subj};%(M{subj}~=0 & abs(M{subj})<mean(abs(M{subj}(M{subj}~=0)))+3*std(abs(M{subj}(M{subj}~=0))));

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
    stat([1:nsel]+nrow,:) = [abs(Mselected) Gselected ones(nsel,1)*Age(subj) ones(nsel,1)*subj ones(nsel,1)*groupsel];
    
    end
end
% end

%%

csvwrite(sprintf('mat_choosing_museum_all_data_for_stat.csv'),stat)

