function file_name=Paris_create_posfile(event_file,new_pos_file) 

% event_file='E:\ALIZEE\EPILEPSY\DATA_SAFE\02141\02141\Neuralynx_files\02141_2014-03-25_10-30';
% new_pos_file='E:\ALIZEE\EPILEPSY\DATA_SAFE\02141\02141\Preproc_explicit\';
% mkdir(new_pos_file)
% This script extracts the behavior from Paris patients and compiles a .pos
% file with them
% It calls 2 scripts :
% - patient_XXXX to create patient's event file (needed to get samples)
% - Extract_behavior_Paris to extract behavior
% clc
% clear all
% origFormat = get(0, 'format');
% format('long');

% set(0,'DefaultFigureColor','w')
% scrsz = get(0,'ScreenSize');
% set(0, 'DefaultFigurePosition', [scrsz(3)/16 scrsz(4) scrsz(3)/8 scrsz(4)/2]);
%
% do_extract=1; % 1 to redo extraction, 0 to skip it
% curdir=pwd;
% Dir ='/mnt/data2/5-EPILEPSY';
% Data_dir =[ Dir filesep 'Data'];
%
% Sub_dir = dir([Data_dir filesep 'Sub*']);
%
% cd([Dir filesep 'Analysis' filesep 'Behavior'])
% load Sessions


load('E:\ALIZEE\EPILEPSY\DATA_PREPROCESS\Behavior_complete\all_sub.mat')
load('E:\ALIZEE\EPILEPSY\Analysis\Behavior\Extraction\Sessions.mat')
load('E:\ALIZEE\EPILEPSY\Analysis\Behavior\Extraction\Sub_names.mat')

index_name_sub=find(ismember(new_pos_file,'\'),5,'first');
name_sub=new_pos_file((index_name_sub(end-1)+1):index_name_sub(end)-1);
sub=find(strcmp([Sub_Names],name_sub));


path_neuralynx=event_file;

disp('--------------------------------------------------------------------------------------')
disp(' ==================                   Convert ncs to eeg               ================== ')
%  convert_cns2eeg(path_neuralynx)

disp('--------------------------------------------------------------------------------------')
disp(' ==================                   Getting Events               ================== ')

event=read_neuralynx_files(path_neuralynx);


evt_cell = struct2cell(event);
code     = cell2mat(evt_cell(2,:)');
% sample   = 512*double(cell2mat(evt_cell(3,:)'))/4000; % samples in 512 Hz for eeg2env2
sample   = 512*double(cell2mat(evt_cell(3,:)'))/4000; % samples in 512 Hz for eeg2env2

sample2  = double(cell2mat(evt_cell(3,:)'))/4; % samples in ms for comparison with .mat

% show sessions visually
code(code~=255)=0;
last=find(code==255,1,'last');
% figure; plot(sample2(1:last)/64, code(1:last));
sample(code~=255)=[];
sample2(code~=255)=[];
code(code~=255)=[];



%% TOTAL EVENT should be 10 (training) + 2*60 (age) + 10 (training) + 3*60 (pleasantness) + 3*30 (choices)
% Get the real number in behavioral data :
nb_trials=[length(Data_trials{1}{1}{sub}.Age_ID) length(Data_trials{1}{2}{sub}.Age_ID) ...
    length(Data_trials{2}{1}{sub}.Pleas_ID) length(Data_trials{2}{2}{sub}.Pleas_ID) length(Data_trials{2}{3}{sub}.Pleas_ID)...
    length(Data_trials{3}{1}{sub}.Left_choice) length(Data_trials{3}{2}{sub}.Left_choice) length(Data_trials{3}{3}{sub}.Left_choice)];


Onsets_rat = [Data_trials{1}{1}{sub}.Age_cross_onset,Data_trials{1}{2}{sub}.Age_cross_onset,...
    Data_trials{2}{1}{sub}.Pleas_cross_onset,Data_trials{2}{2}{sub}.Pleas_cross_onset,Data_trials{2}{3}{sub}.Pleas_cross_onset];

Onsets_ch = [Data_trials{3}{1}{sub}.Choice_cross_onset,Data_trials{3}{2}{sub}.Choice_cross_onset,Data_trials{3}{3}{sub}.Choice_cross_onset];

%% IDENTIFY SESSIONS WITH RT CORRELATIONS :
disp('--------------------------------------------------------------------------------')
disp(' ==================                Checking Events            ==================')
figure;
for task=1:5
    corr(2)=0;
    index=0;
    while corr(2)<=0.98
        index=index+1;
        corr=corrcoef(Onsets_rat(2:nb_trials(task),task)-Onsets_rat(1:nb_trials(task)-1,task),sample2(index+1:index+nb_trials(task)-1)-sample2(index:index+nb_trials(task)-2));
        
        %                  figure; scatter(Onsets_rat(2:60,task)-Onsets_rat(1:59,task),sample2(index+1:index+59)-sample2(index:index+58));
    end
    subplot(2,4,task)
    scatter(Onsets_rat(2:nb_trials(task),task)-Onsets_rat(1:nb_trials(task)-1,task),sample2(index+1:index+nb_trials(task)-1)-sample2(index:index+nb_trials(task)-2));
    indexes_task(task,:)=index:index+nb_trials(task)-1;
    onset_cross_task(task,:)=sample(index:index+nb_trials(task)-1);
    disp(['Task ' num2str(task) ' is ok'])
end

nb_trial_choice=nb_trials(6);
for task=1:3
    corr(2)=0;
    index=0;
    while corr(2)<=0.98
        index=index+1;
        corr=corrcoef(Onsets_ch(2:nb_trial_choice,task)-Onsets_ch(1:nb_trial_choice-1,task),sample2(index+1:index+nb_trial_choice-1)-sample2(index:index+nb_trial_choice-2));
        
        %                 figure; scatter(Onsets_rat(2:60,task)-Onsets_rat(1:59,task),sample2(index+1:index+59)-sample2(index:index+58));
    end
    subplot(2,4,task+5)

     scatter(Onsets_ch(2:nb_trial_choice,task)-Onsets_ch(1:nb_trial_choice-1,task),sample2(index+1:index+nb_trial_choice-1)-sample2(index:index+nb_trial_choice-2));
    indexes_task(task+5,:)=[index:index+nb_trial_choice-1 NaN(1,60-nb_trial_choice)];
    onset_cross_task(task+5,:)=[sample(index:index+nb_trial_choice-1)' NaN(1,60-nb_trial_choice)];
    disp(['Task ' num2str(task+5) ' is ok'])

end

close all
disp('--------------------------------------------------------------------------------')
disp(' ==================                Creating Pos File           =================')

data_pref1_sub=Data_trials{2}{1}{sub};
index1=data_pref1_sub.Pleas_ID;
data_age1_sub=Data_trials{1}{1}{sub};

data_pref2_sub=Data_trials{2}{2}{sub};
index2=data_pref2_sub.Pleas_ID;
data_age2_sub=Data_trials{1}{2}{sub};

for t=1:60
    pref1(t)=data_pref1_sub.Pleas(find(data_age1_sub.Age_ID(t)==index1));
end
data_pref2_sub=Data_trials{2}{2}{sub};
index2=data_pref2_sub.Pleas_ID;

for t=1:60
    pref2(t)=data_pref2_sub.Pleas(find(data_age2_sub.Age_ID(t)==index2));
end


pref_task=[pref1;pref2];
pref_task(3,:)=Data_trials{2}{1}{sub}.Pleas;
pref_task(4,:)=Data_trials{2}{2}{sub}.Pleas;
pref_task(5,:)=Data_trials{2}{3}{sub}.Pleas;

conf_task(1,:)=Data_trials{1}{1}{sub}.Age_conf;
conf_task(2,:)=Data_trials{1}{2}{sub}.Age_conf;
conf_task(3,:)=Data_trials{2}{1}{sub}.Pleas_conf;
conf_task(4,:)=Data_trials{2}{2}{sub}.Pleas_conf;
conf_task(5,:)=Data_trials{2}{3}{sub}.Pleas_conf;

RT_task(1,:)=Data_trials{1}{1}{sub}.Age_rt;
RT_task(2,:)=Data_trials{1}{2}{sub}.Age_rt;
RT_task(3,:)=Data_trials{2}{1}{sub}.Pleas_rt;
RT_task(4,:)=Data_trials{2}{2}{sub}.Pleas_rt;
RT_task(5,:)=Data_trials{2}{3}{sub}.Pleas_rt;

code_task(1,:)=11*ones(1,size(Data_trials{1}{1}{sub}.Age_conf,1));
code_task(2,:)=21*ones(1,size(Data_trials{1}{2}{sub}.Age_conf,1));
code_task(3,:)=31*ones(1,size(Data_trials{2}{1}{sub}.Pleas_conf,1));
code_task(4,:)=41*ones(1,size(Data_trials{2}{2}{sub}.Pleas_conf,1));
code_task(5,:)=51*ones(1,size(Data_trials{2}{3}{sub}.Pleas_conf,1));

start = 512*(Data_trials{1}{1}{sub}.Age_onset-Data_trials{1}{1}{sub}.Age_cross_onset)/1000;
onset_task(1,:)=onset_cross_task(1,:)+start';

start = 512*(Data_trials{1}{2}{sub}.Age_onset-Data_trials{1}{2}{sub}.Age_cross_onset)/1000;
onset_task(2,:)=onset_cross_task(2,:)+start';

start = 512*(Data_trials{2}{1}{sub}.Pleas_onset-Data_trials{2}{1}{sub}.Pleas_cross_onset)/1000;
onset_task(3,:)=onset_cross_task(3,:)+start';

start = 512*(Data_trials{2}{2}{sub}.Pleas_onset-Data_trials{2}{2}{sub}.Pleas_cross_onset)/1000;
onset_task(4,:)=onset_cross_task(4,:)+start';

start = 512*(Data_trials{2}{3}{sub}.Pleas_onset-Data_trials{2}{3}{sub}.Pleas_cross_onset)/1000;
onset_task(5,:)=onset_cross_task(5,:)+start';

for t=1:5
median_high(t,pref_task(t,:)>median(pref_task(t,:)))=101*ones(sum(pref_task(t,:)>median(pref_task(t,:))),1);
median_low(t,pref_task(t,:)<median(pref_task(t,:)))=100*ones(sum(pref_task(t,:)<median(pref_task(t,:))),1);
end
code_median=(median_high+median_low)';

onset_cross_task=onset_cross_task(1:5,:)';
onset_task=onset_task';
code_task=code_task';
pref_task=pref_task';
conf_task=conf_task';


pos1=[onset_task(:) code_task(:) pref_task(:) conf_task(:)];
pos2=[onset_task(:) code_median(:) pref_task(:) conf_task(:)];
pos=[pos1;pos2];


[a,b]=sort(pos(:,1));
pos_corrected=pos(b,:);%done: we now have all events of interest in the posfile :)


%%Raphaelle weird script : padding ?

%     final_pos(:,1)=final_pos(:,1)-final_pos(1,1)+padding*512+1;
% pos_corrected(:,1)=pos_corrected(:,1)-sample(find(sample~=0,1,'first'));
%  pos_corrected(:,1)=pos_corrected(:,1)-pos_corrected(1,1)+1;
% pos_corrected(:,1)=pos_corrected(:,1)-sample(1);

%     final_posZ(:,1)=final_posZ(:,1)-final_posZ(1,1)+padding*512+1;
%     final_posZ(:,1)=final_posZ(:,1)-sample(1);
%
%      pos_corrected=[pos_corrected(:,1)/8, pos_corrected(:,2:end)];



file_name=[new_pos_file filesep name_sub '.pos'];
dlmwrite([file_name],pos_corrected,'delimiter','\t','precision', '%2e');


