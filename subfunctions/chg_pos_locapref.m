function chg_pos_locapref(a_oldposfile,a_logfile)

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
[a_code v_ttime pref_score stim_id tmp]=extract_LOCA_pref_datafromlogfile(a_logfile);%ALIZEE: here I take events + timings from the .logfile (file from stimulus-presentation PC)
%from the logfile; we can quickly extract the timing of the 120 rating
%onset = time at which the stimulus is displayed by the presentation PC
cue_latency_pcstim=v_ttime(tmp);%timing in milliseconds for the PC stim clock
id_fix=find(strcmp(a_code,{'50'}));%index in the logfile corresponding to fixations
id_resplog=find(strcmp(a_code,{'5'}));
fix_latency_pcstim=v_ttime(id_fix);%in milliseconds
resp_latency_pcstim=v_ttime(id_resplog);%in milliseconds

%we will need cue latency relative to fix latency in eeg samples later if
%there is a missing cue in the eeg file...
cue_rel_time_milliseconds=cue_latency_pcstim-fix_latency_pcstim;%should be around 1.5 second +- 100 ms

if size(resp_latency_pcstim,1)>120
    find(id_fix(2:end)-3~=id_resplog(1:end-2),1,'first')
end

resp_rel_time_milliseconds=[fix_latency_pcstim(2:end)-resp_latency_pcstim(1:end-1); mean(fix_latency_pcstim(2:end)-resp_latency_pcstim(1:end-1))];



%WITH THESE INFO WE CAN CORRECT THE POSFILE (if necessary)%%%%%%%%%%%%

%Open + quick check of the POSFILE
a_newposfile =strrep(a_oldposfile,'.pos','_locapref.pos');%
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
id_fixpos=find(m_pos(:,2)==50);%OK: 120 events
id_cue=find(m_pos(:,2)==51);%almost OK (119 events): we could of course use .log data to back-up the missing event;
id_resp=find(m_pos(:,2)==5);
cue_latency_eeg=m_pos(id_cue,1);%in samples... we need to get s_fs to convert samples into milliseconds: this info is in the .ent file
fix_latency_eeg=m_pos(id_fixpos,1);%
resp_latency_eeg=m_pos(id_resp,1);


%% CHECK FOR TOO MANY TRIALS
pos_from_log=create_pos_from_log(a_code, v_ttime,fix_latency_pcstim,cue_latency_pcstim);

find((m_pos(2:end,2)==51 & m_pos(1:end-1,2)~=5 & m_pos(1:end-1,2)~=50)) 

%% RAPPORT DE CORRECTION DU POS FILE

[rapport, rap]=check_data(id_fixpos,id_cue,id_resp);

%%
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
            if isempty (find(cue_latency_eeg>fix_latency_tmp(1) & cue_latency_eeg<fix_latency_tmp(2))) %just check there is no cue event between these two fix events
                cue_onsets=fix_latency_tmp+cue_rel_time_tmp;
                cue_latency_eeg(id_trial_no_cue);
                new_pos_line=[cue_onsets(:,1) 51*ones(length(cue_onsets(:,1)),1) zeros(length(cue_onsets(:,1)),1)];
                m_pos_corrected=[m_pos ;new_pos_line];
                %now re-order m_pos with respect to timings of column 1
                [a,b]=sort(m_pos_corrected(:,1));
                m_pos_corrected=m_pos_corrected(b,:);%done: we now have all events of interest in the posfile :)
            end
        end
        
        id_fixpos=find(m_pos_corrected(:,2)==50);%OK: 120 events
        id_cue=find( m_pos_corrected(:,2)==51);%almost OK (119 events): we could of course use .log data to back-up the missing event;
        id_resp=find( m_pos_corrected(:,2)==5);
        
    end
else
    m_pos_corrected=m_pos;
end
cue_latency_eeg=m_pos(id_cue,1);%in samples... we need to get s_fs to convert samples into milliseconds: this info is in the .ent file
fix_latency_eeg=m_pos(id_fixpos,1);%
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
                new_pos_line=[fix_onsets(:,1) 50*ones(length(fix_onsets(:,1)),1) zeros(length(fix_onsets(:,1)),1)];
                m_pos_corrected=[m_pos_corrected ;new_pos_line];
                %now re-order m_pos with respect to timings of column 1
                [a,b]=sort(m_pos_corrected(:,1));
                m_pos_corrected=m_pos_corrected(b,:);%done: we now have all events of interest in the posfile :)
%             end
        end
        
        id_fixpos=find(m_pos_corrected(:,2)==50);%OK: 120 events
        id_cue=find( m_pos_corrected(:,2)==51);%almost OK (119 events): we could of course use .log data to back-up the missing event;
        id_resp=find( m_pos_corrected(:,2)==5);
    end
    
else
    m_pos_corrected=m_pos;
    
end
cue_latency_eeg=m_pos_corrected(id_cue,1);%in samples... we need to get s_fs to convert samples into milliseconds: this info is in the .ent file
fix_latency_eeg=m_pos_corrected(id_fixpos,1);%
if rap(3)~=1
    %% MISSING RESPONSE
    id_trial_no_resp=[];
    %let's correct that!
    if length(id_resp)<120
        
        if length(id_fixpos)>=120 %if we have all cue: it's quite easy:
            %             id_trial_no_resp=find(~ismember(id_fixpos(2:end)-1,id_resp));%identify trial number for which we did not have cue and pick up from logfile the relative latency of the stim event + resp event during this trial
            id_trial_no_resp=find(m_pos_corrected(id_fixpos(2:end)-1,2)~=5);
            resp_latency_tmp=[fix_latency_eeg(id_trial_no_resp(1:end)+1)];% fix_latency_eeg(id_trial_no_resp(1:end)+2)];
            resp_rel_time_milliseconds_tmp=[resp_rel_time_milliseconds(id_trial_no_resp(1:end))];%  resp_rel_time_milliseconds(id_trial_no_resp(1:end)+1) ];
            %we need that info in EEG samples: use s_fs now
            resp_rel_time_tmp=round(resp_rel_time_milliseconds_tmp*s_fs/1000);
            %                if isempty (find(fix_latency_eeg>resp_latency_tmp(1) & fix_latency_eeg<resp_latency_tmp(2))) %just check there is no cue event between these two fix events
            
            
            
            resp_onsets=resp_latency_tmp-200;
            %                 trucs=m_pos_corrected(id_fixpos(id_trial_no_resp(2:end))-1);
            %                 trucs2=m_pos_corrected(id_fixpos(id_trial_no_resp(2:end)));
            %                 val_diff=(trucs2-trucs)/2;
            %
            % %                 resp_onsets(find(resp_onsets(1:end-1,1)>=m_pos_corrected(id_fixpos(id_trial_no_resp(2:end))-2)))=trucs(find(resp_onsets(1:end-1,1)>=m_pos_corrected(id_fixpos(id_trial_no_resp(2:end))-2)))+val_diff(find(resp_onsets(1:end-1,1)>=m_pos_corrected(id_fixpos(id_trial_no_resp(2:end))-2)));
            %                 resp_onsets(find(resp_onsets(1:end-1,1)>=m_pos_corrected(id_fixpos(id_trial_no_resp(2:end))-1)))=trucs(find(resp_onsets(1:end-1,1)>=m_pos_corrected(id_fixpos(id_trial_no_resp(2:end))-1)))+162.5;
            
            %             fix_latency_eeg(id_trial_no_fix);
            new_pos_line=[resp_onsets(:,1) 5*ones(length(resp_onsets(:,1)),1) zeros(length(resp_onsets(:,1)),1) 1000*ones(length(resp_onsets(:,1)),1)];
            m_pos_corrected=[m_pos_corrected zeros(size(m_pos_corrected,1),1) ;new_pos_line];
            %now re-order m_pos with respect to timings of column 1
            [a,b]=sort(m_pos_corrected(:,1));
            m_pos_corrected=m_pos_corrected(b,:);%done: we now have all events of interest in the posfile :)
            %                end
        end
        
        id_fixpos=find(m_pos_corrected(:,2)==50);%OK: 120 events
        id_cue=find( m_pos_corrected(:,2)==51);%almost OK (119 events): we could of course use .log data to back-up the missing event;
        id_resp=find( m_pos_corrected(:,2)==5);
    else
            m_pos_corrected=[m_pos_corrected zeros(size(m_pos_corrected,1),1)];

    end
else
    m_pos_corrected=m_pos;   
    m_pos_corrected=[m_pos_corrected zeros(size(m_pos_corrected,1),1)];
end

%% 


%% CHECK IF LAST ANSWER IS MISSING OR NOT

if ~ismember(5,m_pos_corrected(find(m_pos_corrected(:,2)==50,1,'last'):find(m_pos_corrected(:,2)==80,1,'first'),2))  %% LAST answer is missing
    timing_stim=m_pos_corrected(find(m_pos_corrected(:,2)==51,1,'last'),1);
    approx_time_until_answer=(pos_from_log(find(pos_from_log(:,2)==5,1,'last'),1)-pos_from_log(find(pos_from_log(:,2)==51,1,'last'),1))*s_fs/1000;
    last_answer=[timing_stim+approx_time_until_answer 5 0 0];
    m_pos_corrected=[m_pos_corrected ; last_answer];
    [a,b]=sort(m_pos_corrected(:,1));
    m_pos_corrected=m_pos_corrected(b,:);
end

% MARK INFERED RESPONSE TIME
 m_pos_corrected(m_pos_corrected(:,4)==1000,1)=NaN(1,sum(m_pos_corrected(:,4)==1000));

 % Nettoyage d'artefacts
m_pos_corrected(~ismember(m_pos_corrected(:,2),[50 51 2 4 3 1 5 80 81]),:)=[];
m_pos_corrected(find(m_pos_corrected(1:end-1,2)==5 & m_pos_corrected(2:end,2)~=50)+1,:)=[];

%% CORRECTION FOR MISTAKE IN RESPONSE : look for doublon and delete the first one

m_pos_corrected(find(m_pos_corrected(1:end-1,2)==5 & m_pos_corrected(2:end,2)==5),:)=[];
%%
 
        id_fixpos=find(m_pos_corrected(:,2)==50);%OK: 120 events
        id_cue=find( m_pos_corrected(:,2)==51);%almost OK (119 events): we could of course use .log data to back-up the missing event;
        id_resp=find( m_pos_corrected(:,2)==5);
[rapport, rap]=check_data(id_fixpos,id_cue,id_resp);


%% CONTROL FOR CORRESPONDANCE BETWEEN POS AND LOG files (reaction times)

% give trial number to the pos file

index_base=50;

m_pos_corrected(:,5)=zeros(size(m_pos_corrected,1),1);
m_pos_corrected(find(m_pos_corrected(:,2)==index_base),5)=1:120;
trial=1;
for i=2:length(m_pos_corrected)
    if m_pos_corrected(i,5)~=trial+1
    m_pos_corrected(i-1,5)=trial;
    else
        m_pos_corrected(i-1,5)=trial;
        trial=trial+1;
    end
end

% LAST CHECK: GET RT_pos and RT_log

for trial=1:119
    RT_pos(trial)=(m_pos_corrected(find(m_pos_corrected(:,5)==trial+1,1,'first'),1)-m_pos_corrected(find(m_pos_corrected(:,5)==trial,1,'first'),1))/s_fs;
    RT_log(trial)=(pos_from_log(find(pos_from_log(:,6)==trial+1,1,'first'),1)-pos_from_log(find(pos_from_log(:,6)==trial,1,'first'),1))/1000;
    
end
% figure;scatter(RT_log(1:119),RT_pos(1:119))
[betas stat]=robustfit(RT_log(1:119),RT_pos(1:119));
if round(betas(2))==1 && round(betas(1))==0 && stat.p(1)>0.05
    display('All trials are ok, pos and log match')
else
    weird=1;
    figure;scatter(RT_log(1:119),RT_pos(1:119))
    % Identifity outliers
%     trials_ok=find(abs(round(RT_pos*10)-round(RT_log*10))<2);
%     figure;scatter(RT_log(trials_ok),RT_pos(trials_ok))
%     outliers=find(abs(round(RT_pos*10)-round(RT_log*10))>2);
error('Timing is unreliable between pos and log, careful with this subject, consider correct the pos file manually')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ALIZEE, on peut surement super améliorer le code très grandement ICI (en "réparant" le
%fichier .pos grâce au fichier .log" au lieu de simplement virer les
%essais... On en reparle ensemble !
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%NEXT: ADD INFO IN THE POSFILE to have a quick acess to different REGRESSORS

clear m_pos
m_pos=m_pos_corrected(:,1:2);



m_pos_tmp=[];
% m_pos_tmp=m_pos(id_cue(11:end),:);
m_pos_tmp=m_pos(id_cue,:);
m_pos_tmp(:,3)=pref_score'; %Ok this is one of the regressor we need... BUT NOT THE ONLY ONE !!
m_pos_tmp=[ m_pos_tmp pref_score'.^2]; %en colonne 4: la valeur au carré comem régresseur pour la confiance
%il ne manque plus que le TR:
for i=1:length(id_cue)
    resp=[];
    if i~=length(id_cue)
        resp=find(m_pos_corrected(id_cue(i):id_cue(i+1),2)==5);
    else
        if i==length(id_cue)
            resp=find(m_pos_corrected(id_cue(i):end,2)==5);
        end
    end
    
    if isempty(resp)
        %a possibility is that instead of event 5  you actually get a
        %prefscore measure (I had to add +10 to the actual prefscore so
        %that pref_socre =0 corresponds to event code 11
        if i~=length(id_cue)
            resp=find(m_pos_corrected(id_cue(i):id_cue(i+1),2)==50)-1;
        else
            resp=find(m_pos_corrected(id_cue(i):end,2)==50)-1;
        end
    end
    
    if m_pos_corrected(id_cue(i),2)==m_pos_corrected(id_cue(i)+1,2) %% One trial is missing
        resp_lat=NaN;
        RT(i)=NaN; %TR en sample = OK comme régresseur 3
    else
        resp_lat=m_pos_corrected(id_cue(i)-1+resp(1),1);
        RT(i)=resp_lat-m_pos_corrected(id_cue(i),1); %TR en sample = OK comme régresseur 3
    end
end

RT=RT./s_fs;
m_pos_tmp=[ m_pos_tmp RT']; %en colonne 5: la TR comme régresseur 3


%ce qui serait très bien c'est de modifier aussi le .pos pour avoir les
%median split (ajouter un event code pour les essais où la valeur est >median


%MEDIAN SPLIT EVENT CODES: create events 52 for low value trials and events
%53 for high value trials using the median of pref_score///
id_minpref=find(pref_score<=median(pref_score));
id_maxpref=find(pref_score>median(pref_score));
median_split_trials=m_pos_tmp(id_minpref,:);
median_split_trials(:,2)=52;
tmp=m_pos_tmp(id_maxpref,:);tmp(:,2)=53;
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





