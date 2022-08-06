% This script is to create figure for Choosing results from children data
% collected from the Museum of Science Boston.
% Written by Se-Woong Park
% 3mouse_button_MoS.mat and AgeGender449.mat should be in the same folder.
load 3mouse_button_MoS


colorLab = [0 53 200]./255;
colorMuseum = [0 171 240]./255;
lblue = colorMuseum;
lred = colorMuseum;

load AgeGender449.mat
Male = find(Age>0 & Age<130 & Gender=='M' & Hist_Mov=='N' & Hist_Psych=='N')';
Female = find(Age>0 & Age<130 & Gender=='F' & Hist_Mov=='N' & Hist_Psych=='N')';
% subjsel = [Male Female]; %MR20220721
figure



for subj = [Male Female],
    
    Mselected = M{subj};
    if length(Mselected)>=3,
        ERROR(subj,1) = length(find(Mselected==1))/30*100;

    else
        ERROR(subj,1) = NaN;
    end

end


mERROR = ERROR;
plot(Age(Male),mERROR(Male),'o','markerfacecolor',lblue,'markeredgecolor','none','linewidth',2,'markersize',10)
hold on
plot(Age(Female),mERROR(Female),'o','markerfacecolor',lred,'markeredgecolor','none','linewidth',2,'markersize',10)


AD = mERROR(Age>=18 & Age<22 & mERROR>0);
meanAD = mean(AD);
sdAD = std(AD);
hold on
CORR = find(Age>=5 & Age<13 & mERROR>0);
[r p] = corr(Age(CORR), mERROR(CORR))
fcn = 'a*x +b';
fcnexp = 'a*exp(-b*x) +c';
[F GoF] = fit(Age(CORR),mERROR(CORR),fcn,'start',[-1 1]);
[Fexp GoFexp] = fit(Age(CORR),mERROR(CORR),fcnexp,'start',[-150 1/100 20]);

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
ag = 6:2:16;
subjsel = []; %MR20220721

for a = 6:2:16,
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
    subjsel = [subjsel SUBJ']; %MR20220721
    group{a} = SUBJ; %MR20220802

end

% subjsel = [SUBJ]; %MR20220721

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
    group{a} = SUBJ; %MR20220802
end
subjsel = [subjsel SUBJ']; %MR20220721

plot([6:2:13 20],AGE([6:2:13 20]),'ok','linewidth',3,'markersize',7,'markerfacecolor',[0 0 0])


set(gca,'fontsize',24,'box','off')
xlabel('Age (yr)','fontweight','bold')
xlim([5 22])
%ylabel('Success Rate','fontweight','bold')

ylim([0 120])
set(gca,'ytick',0:20:100)


%% MR20220802
% for a = ag,
    
%     group{a} = subjsel(Age(subjsel)>=a-1 & Age(subjsel)<a+1);
   
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
    stat([1:nsel]+nrow,:) = [Mselected Gselected ones(nsel,1)*Age(subj) ones(nsel,1)*subj ones(nsel,1)*groupsel];
    end
    
end


%%
csvwrite(sprintf('mat_choosing_museum_children_data_for_stat.csv'),stat)
