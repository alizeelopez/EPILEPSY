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
id_fix=find(ismember(a_code,{'50'}));%index in the logfile corresponding to fixations
fix_latency_pcstim=v_ttime(id_fix);%in milliseconds
%we will need cue latency relative to fix latency in eeg samples later if
%there is a missing cue in the eeg file...
cue_rel_time_milliseconds=cue_latency_pcstim-fix_latency_pcstim;%should be around 1.5 second +- 100 ms

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

%% RAPPORT DE CORRECTION DU POS FILE

rapport=check_data(id_fixpos,id_cue,id_resp);
%%

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
            new_pos_line=[cue_onsets(1) 51 0];
            m_pos_corrected=[m_pos ;new_pos_line];
            %now re-order m_pos with respect to timings of column 1
            [a,b]=sort(m_pos_corrected(:,1));
            m_pos_corrected=m_pos_corrected(b,:);%done: we now have all events of interest in the posfile :)
        end
    end
    id_fixpos=find( m_pos_corrected(:,2)==50);%OK: 120 events
    id_cue=find( m_pos_corrected(:,2)==51);%almost OK (119 events): we could of course use .log data to back-up the missing event;
    id_resp=find( m_pos_corrected(:,2)==5);
else
    m_pos_corrected=m_pos;
    
end



% if length(id_fixpos)~=120
%     % Load raw posfile
%     a_rawposfile =[a_oldposfile(1:end-4),'_pref.pos'];%
%     raw_pos=load(a_rawposfile); %here are events from the brain PC (with missing events, sometimes...)
%     % pref_given=raw_pos(:
%
% end

if length(id_resp)>120 && length(id_fixpos)==120
    id_resp=[id_fixpos(2:end)-2; find(m_pos(:,2)==5,1,'last')];
end

rapport=check_data(id_fixpos,id_cue,id_resp);



%% MISE EN COMMUN POS(PC brain) ET LOG(PC stim)
i=0;
j=0;
for ii=1:length(a_code)
    m_pos_from_log(i+1,1)=0;
    if strcmp(a_code{ii},'50')
        i=i+1;
        m_pos_from_log(i,2)=50;
        m_pos_from_log(i,1)=cue_latency_pcstim(j+1);
        if ~ismember('N',a_code{ii+2})
            display('Image ID is missing in the log file')
            i=i+1;
            m_pos_from_log(i,2)=51;
            m_pos_from_log(i,1)=fix_latency_pcstim(j+1);
        end
    elseif ismember('N',a_code{ii})
        i=i+1;
        m_pos_from_log(i,2)=51;
        m_pos_from_log(i,1)=fix_latency_pcstim(j+1);
    elseif strcmp('5',a_code{ii})
        i=i+1;
        j=j+1;
        m_pos_from_log(i,2)=5;
        m_pos_from_log(i,3)=str2double(a_code{ii+1}(5:end));
        pref_from_log(j,1)=str2double(a_code{ii+1}(5:end));
        
    elseif strcmp('2',a_code{ii})
        i=i+1;
        m_pos_from_log(i,2)=2;
    elseif strcmp('4',a_code{ii})
        i=i+1;
        m_pos_from_log(i,2)=4;
    elseif strcmp('1',a_code{ii})
        i=i+1;
        m_pos_from_log(i,2)=1;
    elseif strcmp('3',a_code{ii})
        i=i+1;
        m_pos_from_log(i,2)=3;
    end
    if m_pos_from_log(i,1)==0;
     m_pos_from_log(i,1)=v_ttime(ii);
    end
    
end


new_m_pos=m_pos_corrected;
find(m_pos_corrected(:,2)==51)
m_pos_from_log(find(m_pos_from_log(:,2)==51),3)=pref_score;

if length(find(m_pos_corrected(:,2)==51))==length(find(m_pos_from_log(:,2)==51))
    m_pos_from_log(m_pos_from_log(:,2)==51,4)=m_pos_corrected(m_pos_corrected(:,2)==51,1);
%     figure;scatter(m_pos_from_log(m_pos_from_log(:,4)~=0,1),m_pos_from_log(m_pos_from_log(:,4)~=0,4))
end

time_brain_min_pc=m_pos_from_log(2,4)-m_pos_from_log(2,1);

m_pos_from_log(:,5)=m_pos_from_log(:,4)*1000/s_fs;


%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ALIZEE, on peut surement super améliorer le code très grandement ICI (en "réparant" le
%fichier .pos grâce au fichier .log" au lieu de simplement virer les
%essais... On en reparle ensemble !
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%NEXT: ADD INFO IN THE POSFILE to have a quick acess to different REGRESSORS





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





