% This script extracts the behavior from Paris patients and compiles a .pos
% file with them
% It calls 2 scripts :
% - patient_XXXX to create patient's event file (needed to get samples)
% - Extract_behavior_Paris to extract behavior
clc
clear all
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
totaltrials=60;

subjects=4;
Subj_names={'','','02135','02141','02161','02171'};


% Flag System:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 11: âge tableau
% 12: confiance âge tableau

% 21: âge visage
% 22: confiance âge visage

% 31: valeur
% 32: confiance valeur tableau

% 41: valeur
% 42: confiance valeur visage

% 51: valeur
% 52: confiance valeur nourriture

% Final File
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Col   Type
% 1     Time
% 2     event code
% 3     Age
% 4     Age confidence
% 5     Pleas
% 6     Pleas confidence
% 7     AutoPleas (pleasantness of age items)
% 8     Autoconf (confidence in that pleasantness measure)
% 9     chosen or unchosen (in subsequent choice session)
load('E:\ALIZEE\EPILEPSY\DATA_PREPROCESS\Behavior_complete\all_sub.mat')
load('E:\ALIZEE\EPILEPSY\Analysis\Behavior\Extraction\Sessions.mat')
load('E:\ALIZEE\EPILEPSY\Analysis\Behavior\Extraction\Sub_names.mat')

pos=1;
n_sub=0;
for sub=subjects
    clear trial_length Trials_lenghts
%     toRun = fullfile([Dir filesep 'Analysis' filesep 'Neurons'], ['patient_' Subj_names{sub}]);
%     run(toRun)
path_neuralynx='E:\ALIZEE\EPILEPSY\DATA_SAFE\02141\02141\Neuralynx_files\02141_2014-03-25_10-30';

event=read_neuralynx_files(path_neuralynx);
    
    n_sub=n_sub+1;
%     cd([Dir filesep 'Analysis' filesep 'Behavior'])
%     Extract_behavior_Paris
    
%     data_sub=Data_trials{task}{stim}{sub};
    
%     cd([Data_dir filesep Sub_dir(sub).name filesep Subj_names{sub}])
    
%     load('Events.mat');
    evt_cell = struct2cell(event);
    code     = cell2mat(evt_cell(2,:)');
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
    nb_trials=[10 length(Data_trials{1}{1}{sub}.Age_ID) length(Data_trials{1}{2}{sub}.Age_ID) ...
        10 length(Data_trials{2}{1}{sub}.Pleas_ID) length(Data_trials{2}{2}{sub}.Pleas_ID) length(Data_trials{2}{3}{sub}.Pleas_ID)...
        length(Data_trials{3}{1}{sub}.Left_choice) length(Data_trials{3}{2}{sub}.Left_choice) length(Data_trials{3}{3}{sub}.Left_choice)];
    
    
    
    total_trials=sum(nb_trials);
    
    % find session start lines in sample vector
    sess_lines = find(sample2(2:length(sample2))-sample2(1:length(sample2)-1)>50000); %% NOT OK
    if length(sess_lines)==0
        error('Sessions are badly defined!!!')
    end
    
%     sess_lines(2:end)-sess_lines(1:end-1)
    sess_lines=sess_lines(2:end-1);
    
    Onsets_rat = [Data_trials{1}{1}{sub}.Age_cross_onset,Data_trials{1}{2}{sub}.Age_cross_onset,...
        Data_trials{2}{1}{sub}.Pleas_cross_onset,Data_trials{2}{2}{sub}.Pleas_cross_onset,Data_trials{2}{3}{sub}.Pleas_cross_onset];

    Onsets_ch = [Data_trials{3}{1}{sub}.Choice_cross_onset,Data_trials{3}{2}{sub}.Choice_cross_onset,Data_trials{3}{3}{sub}.Choice_cross_onset];
% 
%   Onsets_rat = [Data_trials{1}{1}{sub}.Age_onset,Data_trials{1}{2}{sub}.Age_onset,...
%         Data_trials{2}{1}{sub}.Pleas_onset,Data_trials{2}{2}{sub}.Pleas_onset,Data_trials{2}{3}{sub}.Pleas_onset];
% 
%     Onsets_ch = [Data_trials{3}{1}{sub}.Choice_onset,Data_trials{3}{2}{sub}.Choice_onset,Data_trials{3}{3}{sub}.Choice_onset];
%                
        %% IDENTIFY SESSIONS WITH RT CORRELATIONS : 
        
        for task=1:5
            corr(2)=0;
            index=0;
            while corr(2)<=0.98
                index=index+1;
                corr=corrcoef(Onsets_rat(2:60,task)-Onsets_rat(1:59,task),sample2(index+1:index+59)-sample2(index:index+58));
                
%                  figure; scatter(Onsets_rat(2:60,task)-Onsets_rat(1:59,task),sample2(index+1:index+59)-sample2(index:index+58));
            end
            figure; scatter(Onsets_rat(2:60,task)-Onsets_rat(1:59,task),sample2(index+1:index+59)-sample2(index:index+58));
            indexes_task(task,:)=index:index+59;
            onset_cross_task(task,:)=sample(index:index+59);

        end
                      
        nb_trial_choice=30;
        for task=1:3
            corr(2)=0;
            index=0;
            while corr(2)<=0.98
                index=index+1;
                corr=corrcoef(Onsets_ch(2:nb_trial_choice,task)-Onsets_ch(1:nb_trial_choice-1,task),sample2(index+1:index+nb_trial_choice-1)-sample2(index:index+nb_trial_choice-2));
                
%                 figure; scatter(Onsets_rat(2:60,task)-Onsets_rat(1:59,task),sample2(index+1:index+59)-sample2(index:index+58));
            end
            figure; scatter(Onsets_ch(2:nb_trial_choice,task)-Onsets_ch(1:nb_trial_choice-1,task),sample2(index+1:index+nb_trial_choice-1)-sample2(index:index+nb_trial_choice-2));
            indexes_task(task+5,:)=[index:index+nb_trial_choice-1 NaN(1,60-nb_trial_choice)];
            onset_cross_task(task+5,:)=[sample(index:index+nb_trial_choice-1)' NaN(1,60-nb_trial_choice)];
            
        end
    
        close all
    %% AGE TASK
    
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
    
onset_task=onset_task';
code_task=onset_task';
pref_task=pref_task';
conf_task=conf_task';


pos=[onset_task(:) code_task(:) pref_task(:) conf_task(:)];
  
    % find trials and check alignment
    for s=1:length(sess_lines)
        if s<length(sess_lines)
            trial_lines{s}  = sess_lines(s)+1:sess_lines(s+1)-1;
            trial_times{s}  = sample(trial_lines{s});  % in 512 Hz
            trial_times2{s} = sample2(trial_lines{s}); % in ms
        else
            trial_lines{s}  = sess_lines(s)+1:length(sample)-1;
            trial_times{s}  = sample(trial_lines{s});
            trial_times2{s} = sample2(trial_lines{s});
        end
        
        % check trial durations
        trial_length{s} = trial_times2{s}(2:length(trial_times2{s}))-trial_times2{s}(1:length(trial_times2{s})-1);
        
        if isnan(nsess{sub}(s))
            Etrial_lines{s}  = trial_lines{s};
            Etrial_times{s}  = trial_times{s};
            Etrial_length{s} = trial_length{s};
            
            trial_lines{s}  = [];
            trial_times{s}  = [];
            trial_length{s} = [];
            trial_times2{s} = [];
        else
            if nsess{sub}(s)==2
                Etrial_lines{s}  = trial_lines{s}(60:70);
                Etrial_times{s}  = trial_times{s}(60:70);
                Etrial_length{s} = trial_length{s}(60:70);
                
                trial_lines{s}(60:70)  = [];
                trial_times{s}(60:70)  = [];
                trial_times2{s}(60:70) = [];
                trial_length{s}(60:70) = [];
            end
            if and(sub>4, nsess{sub}(s)>5)
                trial_length{s}(30:length(trial_length{s}))=[];
            end
            sess = Order(sub,nsess{sub}(s))/10;
            if any(abs(trial_length{s}(1:29)-Trials_lenghts{sess}(1:29))>2)
                error(['Triggers and mat behavior are not aligned in session ', num2str(s)])
            end
            
            Pos.Code(trial_lines{s},1) = sess*10+1;
            if nsess{sub}(s)<3
                clear start
                start = 512*(Age_onset(:,n_sub,sess)-Age_cross_onset(:,n_sub,sess))/1000; % convert from ms to 512Hz
                Pos.Time(trial_lines{s})      = trial_times{s}+start;
               
                Pos.autoValZ(trial_lines{s})  = zscore(autoPleas(:,n_sub,sess));
                Pos.autoConfZ(trial_lines{s}) = zscore(autoConf(:,n_sub,sess));
                Pos.AgeZ(trial_lines{s})      = zscore(zAge(:,n_sub,sess));
                Pos.ConfageZ(trial_lines{s})  = zscore(Age_conf(:,n_sub,sess));
                
                Pos.autoVal(trial_lines{s})   = (autoPleas(:,n_sub,sess));
                Pos.autoConf(trial_lines{s})  = (autoConf(:,n_sub,sess));
                Pos.Age(trial_lines{s})       = (zAge(:,n_sub,sess));
                Pos.Confage(trial_lines{s})   = (Age_conf(:,n_sub,sess));
                
                Pos.Pref(trial_lines{s})      = Pref_age(:,n_sub,sess);
            elseif and(nsess{sub}(s)>2,nsess{sub}(s)<6)
                clear start
                start = 512*(Pleas_onset(:,n_sub,sess-2)-Pleas_cross_onset(:,n_sub,sess-2))/1000; % convert from ms to 512Hz
                Pos.Time(trial_lines{s})       = trial_times{s}+start;
                
                Pos.PleasZ(trial_lines{s})     = zscore(Pleas(:,n_sub,sess-2));
                Pos.ConfpleasZ(trial_lines{s}) = zscore(Pleas_conf(:,n_sub,sess-2));
                
                Pos.Pleas(trial_lines{s})      = (Pleas(:,n_sub,sess-2));
                Pos.Confpleas(trial_lines{s})  = (Pleas_conf(:,n_sub,sess-2));
                
                Pos.Pref(trial_lines{s})       = Pref_pleas(:,n_sub,sess-2);
            elseif nsess{sub}(s)>5
                clear start
                start = 512*(Choice_onset(1:length(trial_lines{s}),n_sub,sess-5)-Choice_cross_onset(1:length(trial_lines{s}),n_sub,sess-5))/1000; % convert from ms to 512Hz
                Pos.Time(trial_lines{s})       = trial_times{s}+start;
                Pos.autoConf(trial_lines{s})   = abs(Diff_value(1:length(trial_lines{s}), n_sub,sess-5));
                Pos.autoConfZ(trial_lines{s})  = abs(Diff_value(1:length(trial_lines{s}), n_sub,sess-5));
                Pos.Pref(trial_lines{s})       = Renv(1:length(trial_lines{s}), n_sub,sess-5);
            end
        end
    end
    
    final_pos  = [Pos.Time Pos.Code Pos.Age Pos.Confage Pos.Pleas Pos.Confpleas Pos.autoVal Pos.autoConf Pos.Pref];
    final_posZ = [Pos.Time Pos.Code Pos.AgeZ Pos.ConfageZ Pos.PleasZ Pos.ConfpleasZ Pos.autoValZ Pos.autoConfZ Pos.Pref];
    final_pos(isnan(Pos.Code),:)=[];
    final_posZ(isnan(Pos.Code),:)=[];
    
    final_pos(:,1)=final_pos(:,1)-final_pos(1,1)+padding*512+1;
    final_pos(:,1)=final_pos(:,1)-sample(1);
    
    final_posZ(:,1)=final_posZ(:,1)-final_posZ(1,1)+padding*512+1;
    final_posZ(:,1)=final_posZ(:,1)-sample(1);
    
    final_pos_ds8=[final_pos(:,1)/8, final_pos(:,2:9)];
    final_pos_ds8Z=[final_posZ(:,1)/8, final_posZ(:,2:9)];
    
    if pos==1 % BEWARE, saving to pos file in short format rounds your numbers!!!
        cd([Data_dir filesep Sub_dir(sub).name filesep Subj_names{sub} filesep 'eeg_rawdata'])
        dlmwrite([Subj_names{sub} '.pos'],final_pos,'delimiter','\t','precision', '%2e');
        
        cd([Data_dir filesep Sub_dir(sub).name filesep Subj_names{sub} filesep 'R_test'])
        dlmwrite([Subj_names{sub} '.pos'],final_pos,'delimiter','\t','precision', '%2e');
        
        
        cd([Data_dir filesep Sub_dir(sub).name filesep Subj_names{sub} filesep 'eeg_rawdata'])
        dlmwrite([Subj_names{sub} '_ds8_complete.pos'],final_pos_ds8,'delimiter','\t','precision', '%2e');
        
        cd([Data_dir filesep Sub_dir(sub).name filesep Subj_names{sub} filesep 'R_test'])
        dlmwrite([Subj_names{sub} '_ds8_complete.pos'],final_pos_ds8,'delimiter','\t','precision', '%2e');
        
        cd([Data_dir filesep Sub_dir(sub).name filesep Subj_names{sub} filesep 'eeg_rawdata'])
        dlmwrite([Subj_names{sub} '_ds8_completeZ.pos'],final_pos_ds8Z,'delimiter','\t','precision', '%2e');
        
        cd([Data_dir filesep Sub_dir(sub).name filesep Subj_names{sub} filesep 'R_test'])
        dlmwrite([Subj_names{sub} '_ds8_completeZ.pos'],final_pos_ds8Z,'delimiter','\t','precision', '%2e');
    end
    
    cd('/mnt/data2/5-EPILEPSY/Data')
    save(['Total_data_sub' num2str(sub)], ...
        'Age_ID','Age','Age_cross_onset','Age_onset','Age_cursor','Age_1st_move','Age_rt','zAge',...
        'Age_conf_onset','Age_conf','Age_conf_1st_move','Age_conf_cursor','Age_conf_rt','autoPleas', 'autoConf',...
        'Pleas_ID','Pleas','Pleas_cross_onset','Pleas_onset','Pleas_cursor','Pleas_1st_move','Pleas_rt',...
        'Pleas_conf_onset','Pleas_conf','Pleas_conf_cursor','Pleas_conf_1st_move','Pleas_conf_rt',...
        'Pref_pleas','Pref_age','Choice_cross_onset','Choice_onset','Choice_ID1','Choice_ID2','Choice_rt','Left_choice','Renv','Diff_value','Types');
end

% set(0,'format', origFormat);