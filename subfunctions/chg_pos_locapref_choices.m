function chg_pos_locapref_choices(a_oldposfile,a_logfile2)

%this function modifies .pos file to summarize LOCA_PREF data
%column 3 = values at stim onset for each event 51
%column 4 = squared value at stim onset = for each event 51
%column 5 = RT for each event 51
%FUNCTION NOT FINISHED: TO BE MODIFIED
%TO DO 1: recode events to get mediansplits during both blocks of trials
%also we need to add informations for session 2 (choice) for which I did
%not had time to finish the script

%EVENTS CODEs during block 1 (block 2)
%50: FIXATION
%51: STIMULUS (81 = pstim pair)
%5 choice confirmation


%extract INFO from LOGFILE %%%%%%%%%%%%%%%%%%%%
[a_code, v_ttime, choice, stim_id1, stim_id2, tmp]=extract_LOCA_pref_datafromlogfile_choices(a_logfile2);%ALIZEE: here I take events + timings from the .logfile (file from stimulus-presentation PC)
%from the logfile; we can quickly extract the timing of the 120 rating
%onset = time at which the stimulus is displayed by the presentation PC
cue_latency_pcstim=v_ttime(tmp);%timing in milliseconds for the PC stim clock
id_fix=find(strcmp(a_code,{'80'}));%index in the logfile corresponding to fixations
id_resplog=find(strcmp(a_code,{'confidencerate'}))-1;
fix_latency_pcstim=v_ttime(id_fix);%in milliseconds
resp_latency_pcstim=v_ttime(id_resplog);%in milliseconds

%we will need cue latency relative to fix latency in eeg samples later if
%there is a missing cue in the eeg file...
cue_rel_time_milliseconds=cue_latency_pcstim-fix_latency_pcstim;%should be around 1.5 second +- 100 ms

if sum((cue_rel_time_milliseconds>2000))
    cue_rel_time_milliseconds(cue_rel_time_milliseconds>2)=1500;
warning('two time corrections in log file');
end

if size(resp_latency_pcstim,1)>120
    %% HERE TO CREATE AN ERROR
    find(id_fix(2:end)-2~=id_resplog(1:end-2));
end

resp_rel_time_milliseconds=[fix_latency_pcstim(2:end)-resp_latency_pcstim(1:end-1); mean(fix_latency_pcstim(2:end)-resp_latency_pcstim(1:end-1))];



%WITH THESE INFO WE CAN CORRECT THE POSFILE (if necessary)%%%%%%%%%%%%

%Open + quick check of the POSFILE
a_newposfile =strrep(a_oldposfile,'.pos','_locapref2.pos');%
m_pos=load(a_oldposfile); %here are events from the brain PC (with missing events, sometimes...)

%before everything, We will need the sampling rate of the EEGfile
% EXTRACT SAMPLING RATE INFO FROM ENT FILE
a_entfile =strrep(a_oldposfile,'.pos','.eeg.ent');%
f_old=fopen([a_entfile ],'r');
for s_i=1:9
    a_line=fgetl(f_old);
end;
s_fs = 1/str2num(a_line);
s_fs = round(s_fs); % sampling frequency
fclose(f_old);
%DONE :%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%1. CHECK first block of trials: preference ratinngs
%"TRUE" latencies = from the brain files
%(1) Identify cue+ confirmation +Pref scores + keep the other events
id_fixpos=find(m_pos(:,2)==80);%OK: 120 events
id_cue=find(m_pos(:,2)==81);%almost OK (119 events): we could of course use .log data to back-up the missing event;
% id_resp=find(m_pos(:,2)==5);
cue_latency_eeg=m_pos(id_cue,1);%in samples... we need to get s_fs to convert samples into milliseconds: this info is in the .ent file
fix_latency_eeg=m_pos(id_fixpos,1);%
% resp_latency_eeg=m_pos(id_resp,1);


%% CHECK FOR TOO MANY TRIALS
pos_from_log=create_pos_from_log(a_code, v_ttime,fix_latency_pcstim,cue_latency_pcstim);

% find((m_pos(2:end,2)==51 & m_pos(1:end-1,2)~=5 & m_pos(1:end-1,2)~=50));


m_pos_corrected=m_pos;
% Nettoyage d'artefacts
m_pos_corrected(~ismember(m_pos_corrected(:,2),[50 51 2 4 3 1 5 80 81]),:)=[];
% m_pos_corrected(find(m_pos_corrected(1:end-1,2)==5 & m_pos_corrected(2:end,2)~=50)+1,:)=[];
%% RAPPORT DE CORRECTION DU POS FILE
%%
id_fixpos=find(m_pos_corrected(:,2)==80);%OK: 120 events
id_cue=find( m_pos_corrected(:,2)==81);%almost OK (119 events): we could of course use .log data to back-up the missing event;
% id_resp=find( m_pos_corrected(:,2)==5);

[rapport1, rap]=check_data(id_fixpos,id_cue,0);
disp(rapport1(1:2))


if rap(2)~=1
    %% MISSING CUES
    id_trial_no_cue=[];
    %let's correct that!
    if length(id_cue)<120
        %remove abnormal trials for a quick analysis: normally, the missing cue
        %should be just after a fixation (event 50)
        if length(id_fixpos)==120 %if we have all fixation: it's quite easy:
            id_trial_no_cue=find(~ismember(id_fixpos+1,id_cue));%identify trial number for which we did not have cue and pick up from logfile the relative latency of the stim event + resp event during this trial
            fix_latency_tmp=[fix_latency_eeg(id_trial_no_cue) fix_latency_eeg(id_trial_no_cue+1)];
            cue_rel_time_milliseconds_tmp=[cue_rel_time_milliseconds(id_trial_no_cue)  cue_rel_time_milliseconds(id_trial_no_cue+1) ];
            %we need that info in EEG samples: use s_fs now
            cue_rel_time_tmp=round(cue_rel_time_milliseconds_tmp/1000*s_fs);
            %             if isempty (find(cue_latency_eeg>fix_latency_tmp(1) & cue_latency_eeg<fix_latency_tmp(2))) %just check there is no cue event between these two fix events
            cue_onsets=fix_latency_tmp+cue_rel_time_tmp;
            %                 cue_latency_eeg(id_trial_no_cue);
            new_pos_line=[cue_onsets(:,1) 81*ones(length(cue_onsets(:,1)),1) zeros(length(cue_onsets(:,1)),1)];
            m_pos_corrected=[m_pos_corrected ;new_pos_line];
            %now re-order m_pos with respect to timings of column 1
            [a,b]=sort(m_pos_corrected(:,1));
            m_pos_corrected=m_pos_corrected(b,:);%done: we now have all events of interest in the posfile :)
            %             end
        end
        
        id_fixpos=find(m_pos_corrected(:,2)==80);%OK: 120 events
        id_cue=find( m_pos_corrected(:,2)==81);%almost OK (119 events): we could of course use .log data to back-up the missing event;
        id_resp=find( m_pos_corrected(:,2)==5);
        
    end
    % else
    %     m_pos_corrected=m_pos;
end
cue_latency_eeg=m_pos_corrected(id_cue,1);%in samples... we need to get s_fs to convert samples into milliseconds: this info is in the .ent file
fix_latency_eeg=m_pos_corrected(id_fixpos,1);%
if rap(1)~=1
    %% MISSING FIXATIONS
    id_trial_no_fix=[];
    %let's correct that!
    if length(id_fixpos)<120
        
        if length(id_cue)>=120 %if we have all cue: it's quite easy:
            id_trial_no_fix=find(~ismember(id_cue,id_fixpos+1));%identify trial number for which we did not have cue and pick up from logfile the relative latency of the stim event + resp event during this trial
            cue_latency_tmp=[cue_latency_eeg(id_trial_no_fix) cue_latency_eeg(id_trial_no_fix+1)];
            cue_rel_time_milliseconds_tmp=[cue_rel_time_milliseconds(id_trial_no_fix)  cue_rel_time_milliseconds(id_trial_no_fix+1) ];
            %we need that info in EEG samples: use s_fs now
            cue_rel_time_tmp=round(cue_rel_time_milliseconds_tmp/1000*s_fs);
            %             if isempty (find(cue_latency_eeg>cue_latency_tmp(1) & cue_latency_eeg<cue_latency_tmp(2))) %just check there is no cue event between these two fix events
            fix_onsets=cue_latency_tmp-cue_rel_time_tmp;
            %             fix_latency_eeg(id_trial_no_fix);
            new_pos_line=[fix_onsets(:,1) 80*ones(length(fix_onsets(:,1)),1) zeros(length(fix_onsets(:,1)),1)];
            m_pos_corrected=[m_pos_corrected ;new_pos_line];
            %now re-order m_pos with respect to timings of column 1
            [a,b]=sort(m_pos_corrected(:,1));
            m_pos_corrected=m_pos_corrected(b,:);%done: we now have all events of interest in the posfile :)
            %             end
        end
        
        id_fixpos=find(m_pos_corrected(:,2)==80);%OK: 120 events
        id_cue=find( m_pos_corrected(:,2)==81);%almost OK (119 events): we could of course use .log data to back-up the missing event;
    end
    
    % else
    %     m_pos_corrected=m_pos_corrected;
    
end
cue_latency_eeg=m_pos_corrected(id_cue,1);%in samples... we need to get s_fs to convert samples into milliseconds: this info is in the .ent file
fix_latency_eeg=m_pos_corrected(id_fixpos,1);%

%%


id_fixpos=find(m_pos_corrected(:,2)==80);%OK: 120 events
id_cue=find(m_pos_corrected(:,2)==81);%almost OK (119 events): we could of course use .log data to back-up the missing event;
[rapport2, rap]=check_data(id_fixpos,id_cue,0);
disp(rapport2(1:2))


%% CONTROL FOR CORRESPONDANCE BETWEEN POS AND LOG files (reaction times)

% give trial number to the pos file

index_base=80;

m_pos_corrected(:,5)=zeros(size(m_pos_corrected,1),1);
m_pos_corrected(find(m_pos_corrected(:,2)==index_base),5)=1:120;
trial=1;
for i=2:length(m_pos_corrected)
    if i>find(m_pos_corrected(:,2)==index_base,1,'first')
    if m_pos_corrected(i,5)~=trial+1
        m_pos_corrected(i-1,5)=trial;
    else
        m_pos_corrected(i-1,5)=trial;
        trial=trial+1;
    end
    end
end

% LAST CHECK: GET RT_pos and RT_log

for trial=1:119
    RT_pos(trial)=(m_pos_corrected(find(m_pos_corrected(:,5)==trial+1,1,'first'),1)-m_pos_corrected(find(m_pos_corrected(:,5)==trial,1,'first'),1))/s_fs;
    RT_log(trial)=(pos_from_log(find(pos_from_log(:,6)==trial+1,1,'first'),1)-pos_from_log(find(pos_from_log(:,6)==trial,1,'first'),1))/1000;
end
RT_pos(120)=(m_pos_corrected(end,1)-m_pos_corrected(find(m_pos_corrected(:,2)==80,1,'last'),1))/s_fs;
RT_log(120)=(pos_from_log(end-1,1)-pos_from_log(find(pos_from_log(:,2)==80,1,'last'),1))/1000;

% figure;scatter(RT_log(1:119),RT_pos(1:119))
[betas stat]=robustfit(RT_log(1:120),RT_pos(1:120));
if round(betas(2))==1 && round(betas(1))==0 && stat.p(1)>0.05
    outliers=find(abs(round(RT_pos*10)-round(RT_log*10))>2);
    display([num2str(length(outliers)) ' outliers, removed from pos file'])
else
    weird=1;
    figure;scatter(RT_log(1:120),RT_pos(1:120))
    % Identifity outliers
    %     trials_ok=find(abs(round(RT_pos*10)-round(RT_log*10))<2);
    %     figure;scatter(RT_log(trials_ok),RT_pos(trials_ok))
    %     outliers=find(abs(round(RT_pos*10)-round(RT_log*10))>2);
    bar([RT_log(outliers) ; RT_pos(outliers)]')
    error('Timing is unreliable between pos and log, careful with this subject, consider correct the pos file manually')
end

m_pos_corrected(ismember(m_pos_corrected(:,5),outliers),:)=[];
trials_included=~ismember(1:120,outliers);
choice=choice(trials_included);
RT_pos=RT_pos(trials_included);

id_cue=find(m_pos_corrected(:,2)==81);%almost OK (119 events): we could of course use .log data to back-up the missing event;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ALIZEE, on peut surement super améliorer le code très grandement ICI (en "réparant" le
%fichier .pos grâce au fichier .log" au lieu de simplement virer les
%essais... On en reparle ensemble !
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%NEXT: ADD INFO IN THE POSFILE to have a quick acess to different REGRESSORS
index_name_sub=find(ismember(a_oldposfile,'\'),2,'last');
name_sub=a_oldposfile((index_name_sub(1)+1):index_name_sub(2)-1);
load('E:\ALIZEE\EPILEPSY\DATA_PREPROCESS\Behavior_complete\all_sub.mat')
load('E:\ALIZEE\EPILEPSY\Analysis\Behavior\Extraction\Sub_names.mat')

index_behavior=find(strcmp([Sub_Names],name_sub));
data_choice_sub=Data_trials{3}{3}{index_behavior};

choice(choice==3)=ones(sum(choice==3),1);
choice(choice==4)=zeros(sum(choice==4),1);


if sum(data_choice_sub.Left_choice(trials_included)==choice')~=length(choice)
    error('behvaior downloaded does not match to computed data')
end

type=3;
sub=index_behavior;
ID1=data_choice_sub.Choice_ID1;
ID2=data_choice_sub.Choice_ID2;
if ~isempty(find(Data_trials{2}{type}{sub}.Pleas_ID==ID1(1)))
    
    clear V1 V2
    for t=1:length(ID1)
        if isempty(find(Data_trials{2}{type}{sub}.Pleas_ID==ID1(t))) || isempty(find(Data_trials{2}{type}{sub}.Pleas_ID==ID2(t)))
            V1(t)=NaN;
            V2(t)=NaN;
            
        else
            V1(t)=Data_trials{2}{type}{sub}.Pleas(find(Data_trials{2}{type}{sub}.Pleas_ID==ID1(t)));
            V2(t)=Data_trials{2}{type}{sub}.Pleas(find(Data_trials{2}{type}{sub}.Pleas_ID==ID2(t)));
        end
    end
    
end

V1=V1(trials_included);
V2=V2(trials_included);

Vchosen(choice==1)=V1(choice==1);
Vchosen(choice==0)=V2(choice==0);

Vunchosen(choice==1)=V2(choice==1);
Vunchosen(choice==0)=V1(choice==0);

clear m_pos
m_pos=m_pos_corrected(:,1:2);



m_pos_tmp=[];
% m_pos_tmp=m_pos(id_cue(11:end),:);
m_pos_tmp=m_pos(id_cue,:);
m_pos_tmp(:,3)=V1'; %LEFT 
m_pos_tmp(:,4)=V2'; % RIGHT
m_pos_tmp(:,5)=Vchosen'; %chosen 
m_pos_tmp(:,6)=Vunchosen'; % unchosen
m_pos_tmp(:,7)=V1-V2; %LEFT 
m_pos_tmp(:,8)=Vchosen-Vunchosen; % chosen-unchosen




% RT=RT./s_fs;
RT=RT_pos;
m_pos_tmp=[ m_pos_tmp RT']; %en colonne 5: la TR comme régresseur 3


%ce qui serait très bien c'est de modifier aussi le .pos pour avoir les
%median split (ajouter un event code pour les essais où la valeur est >median


%MEDIAN SPLIT EVENT CODES: create events 82 for low Dvalue trials and events
%83 for high Dvalue trials using the median of pref_score///
id_minpref=find(m_pos_tmp(:,8)<=median(m_pos_tmp(:,8)));
id_maxpref=find(m_pos_tmp(:,8)>median(m_pos_tmp(:,8)));
median_split_trials=m_pos_tmp(id_minpref,:);
median_split_trials(:,2)=82;
tmp=m_pos_tmp(id_maxpref,:);tmp(:,2)=83;
median_split_trials=[median_split_trials; tmp];
m_pos_tmp=[m_pos_tmp ;median_split_trials];
%now re-order m_pos with respect to timings of column 1
[a,b]=sort(m_pos_tmp(:,1));
m_pos_tmp=m_pos_tmp(b,:);%done: we now have all events of interest in the posfile :)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%RESTE A FAIRE: LA SESSION 2 avec les choix binaires: ALIZEE ??? J'ai déjà
%pas mal de codes de prêt mais pas eu le temps de nettoyer...


f=fopen(a_newposfile,'w');
for i=1:length(m_pos_tmp)
    fprintf(f,'%d\t%d\t%d\t%d\t%d\n',m_pos_tmp(i,:));
end


fclose(f);





