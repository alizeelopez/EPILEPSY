function chg_pos_main(a_oldposfile,a_logfile)

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

code_event=[10 20 30 40 50];
nb_event=length(code_event);
%extract INFO from LOGFILE %%%%%%%%%%%%%%%%%%%%
[a_code, v_ttime, pref_score, ~, tmp]=extract_LOCA_pref_datafromlogfile(a_logfile);%ALIZEE: here I take events + timings from the .logfile (file from stimulus-presentation PC)
%from the logfile; we can quickly extract the timing of the 120 rating
%onset = time at which the stimulus is displayed by the presentation PC
cue_latency_pcstim=v_ttime(tmp);%timing in milliseconds for the PC stim clock



mat_to_order=[];
for c=1:nb_event
id_fix_ini{c}=find(strcmp(a_code,{num2str(code_event(c))}));%index in the logfile corresponding to fixations
mat_to_order=[mat_to_order id_fix_ini{c}(60)];
end
[~, order] = sort(mat_to_order);
code_event_real=code_event(order);
for c=1:nb_event
id_fix{c}=id_fix_ini{order(c)};
end
% for c=1:3
%     %% Too many fixations ?
%     if size(id_fix{c},1)>60 %% If age 30, 40 or 50 has been given, mistaken for a fix, so remove it
%         if ~isempty(id_fix{c}(find(ismember(a_code(id_resplog{c}(1:end)+1),num2str(code_event_real(c))))+1))
%         a_code{id_fix{c}(find(ismember(a_code(id_resplog{c}(1:end)+1),num2str(code_event_real(c))))+1)}=[num2str(code_event_real(c)) '.1'];
%         id_fix{c}=find(strcmp(a_code,{num2str(code_event_real(c))}));
%         end
%     end
% end

%% In case conf is 30,40 or 50 elsewhere :
if size(id_fix{3},1)>60
    if id_fix{3}(61:end)>id_fix{4}(1)
        id_fix{3}(id_fix{3}(61:end)>id_fix{4}(1) | id_fix{3}(61:end)<id_fix{4}(end))=[];
    end
end

for c=3:nb_event-1
if size(id_fix{c},1)>60
    if id_fix{c}(61:end)>id_fix{c+1}(1) | id_fix{c}(61:end)<id_fix{c-1}(end)
        id_fix{c}(id_fix{c}(61:end)>id_fix{c+1}(1) | id_fix{c}(61:end)<id_fix{c-1}(end))=[];
    end
end
end

if size(id_fix{5},1)>60
    if id_fix{5}(61:end)<id_fix{4}(end) + id_fix{5}(1:length(id_fix{5})-60)<id_fix{4}(end)
        id_fix{5}(id_fix{5}(61:end)<id_fix{4}(end)  | id_fix{5}(1:length(id_fix{5})-60)<id_fix{4}(end))=[];
    end
end


%% NOT OK !
% id_allresplog{1}=find(strcmp(a_code(id_fix{(1)}(1):id_fix{(2)}(1)),{'5'}));
% id_allresplog{2}=find(strcmp(a_code(id_fix{(2)}(1):id_fix{(3)}(1)),{'5'}));
% id_allresplog{3}=find(strcmp(a_code(id_fix{(3)}(1):id_fix{(4)}(1)),{'5'}));
% id_allresplog{4}=find(strcmp(a_code(id_fix{(4)}(1):id_fix{(5)}(1)),{'5'}));
% id_allresplog{5}=find(strcmp(a_code(id_fix{(5)}(1):end),{'5'}));

% for c=1:nb_event
%     %% Too many responses ?
%     if sum(id_allresplog{c}(2:end)-id_allresplog{c}(1:end-1)==1)
%         id_allresplog{c}(find(id_allresplog{c}(2:end)-id_allresplog{c}(1:end-1)==1)+1)=[];
%     end
% end

% if size(id_fix{2},1)>60 %% If age 20 has been given, mistaken for a fix, so remove it
%     index_false_code=find(ismember(a_code(id_allresplog{2}(1:end)+1),'20'))+1;
%     for i=1:length(index_false_code)
%     a_code{(index_false_code(i))}='20.1';
%     end
% end
% id_fix{2}=find(strcmp(a_code,{'20'}));%index in the logfile corresponding to fixations



% for c=1:nb_event
%     id_resplog{c}=id_allresplog{c}(1:2:length(id_fix{c})*2);
%     id_conflog{c}=id_allresplog{c}(2:2:length(id_fix{c})*2);
% end

% id_conflog{1}=find(cell2mat(regexp(a_code(id_fix{1}(1):id_fix{2}(1)),{'PS'}))==1);
% id_conflog{2}=find(cell2mat(regexp(a_code(id_fix{2}(1):id_fix{3}(1)),{'PS'}))==1);
% id_conflog{3}=find(cell2mat(regexp(a_code(id_fix{3}(1):id_fix{4}(1)),{'PS'}))==1);
% id_conflog{4}=find(cell2mat(regexp(a_code(id_fix{4}(1):id_fix{5}(1)),{'PS'}))==1);
% id_conflog{5}=find(cell2mat(regexp(a_code(id_fix{5}(1):end),{'PS'}))==1);
% 

for c=1:nb_event
%     resp_latency_pcstim{c}=v_ttime(id_resplog{c});%in milliseconds
    fix_latency_pcstim{c}=v_ttime(id_fix{c});%in milliseconds
%     cueconf_latency_pcstim{c}=v_ttime(id_conflog{c});%in milliseconds
    
end
%we will need cue latency relative to fix latency in eeg samples later if
%there is a missing cue in the eeg file...

cue_rel_time_milliseconds_total=cue_latency_pcstim-[fix_latency_pcstim{1};fix_latency_pcstim{2};fix_latency_pcstim{3} ;fix_latency_pcstim{4} ;fix_latency_pcstim{5}];%should be around 1.5 second +- 100 ms


if sum((cue_rel_time_milliseconds_total>2000)+(cue_rel_time_milliseconds_total<500))
    cue_rel_time_milliseconds_total(find(cue_rel_time_milliseconds_total>2000 | (cue_rel_time_milliseconds_total<500)))=1500;
    warning('two time corrections in log file');
end
ind=1:60:5*60;
for c=1:nb_event
cue_rel_time_milliseconds{c}=cue_rel_time_milliseconds_total(ind(c):ind(c)+59);
end


% for c=1:nb_event
%     resp_rel_time_milliseconds{c}=[fix_latency_pcstim{c}(2:end)-resp_latency_pcstim{c}(1:end-1); mean(fix_latency_pcstim{c}(2:end)-resp_latency_pcstim{c}(1:end-1))];
%     cueconf_rel_time_milliseconds{c}=[resp_latency_pcstim{c}(2:end)-cueconf_latency_pcstim{c}(1:end-1); mean(resp_latency_pcstim{c}(2:end)-cueconf_latency_pcstim{c}(1:end-1))];
%     
% end


%WITH THESE INFO WE CAN CORRECT THE POSFILE (if necessary)%%%%%%%%%%%%

%Open + quick check of the POSFILE
a_newposfile =strrep(a_oldposfile,'.pos','_locapref_all.pos');%
m_pos=load(a_oldposfile); %here are events from the brain PC (with missing events, sometimes...)

%before everything, We will need the sampling rate of the EEGfile
% EXTRACT SAMPLING RATE INFO FROM ENT FILE
[index_preproc, ind_end]=regexp((a_oldposfile),'\Preproc_explicit\');
if exist(strrep(a_oldposfile([1:index_preproc-1 (ind_end+1):end]),'.pos','_D1.eeg.ent'),'file')
    a_entfile =strrep(a_oldposfile([1:index_preproc-1 (ind_end+1):end]),'.pos','_D1.eeg.ent');%
else
    a_entfile =strrep(a_oldposfile([1:index_preproc-1 (ind_end+1):end]),'.pos','.eeg.ent');%
end
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

for c=1:nb_event
id_fixpos_ini{c}=find(m_pos(:,2)==code_event(c));
end


[~, order] = sort([id_fixpos_ini{1}(1) id_fixpos_ini{2}(1) id_fixpos_ini{3}(1)  id_fixpos_ini{4}(1)  id_fixpos_ini{5}(1)]);

code_event_pos=code_event(order);

for c=1:nb_event
id_fixpos{c}=id_fixpos_ini{order(c)};
id_cue{c}=find(m_pos(:,2)==code_event_pos(c)+1);%almost OK (119 events): we could of course use .log data to back-up the missing event;
% id_confcue{c}=find(m_pos(:,2)==code_event_pos(c)+3);
end

end_age_sess=find(m_pos(:,2)==60 | m_pos(:,2)==70 | m_pos(:,2)==80,1,'first');

% id_resp{1}=id_fixpos{1}(1)-1+find(m_pos(id_fixpos{1}(1):id_fixpos{2}(1),2)==5);
% id_resp{2}=id_fixpos{2}(1)-1+find(m_pos(id_fixpos{2}(1):id_fixpos{3}(1),2)==5);
% id_resp{3}=id_fixpos{3}(1)-1+find(m_pos(id_fixpos{3}(1):id_fixpos{4}(1),2)==5);
% id_resp{4}=id_fixpos{4}(1)-1+find(m_pos(id_fixpos{4}(1):id_fixpos{5}(1),2)==5);
% id_resp{5}=id_fixpos{5}(1)-1+find(m_pos(id_fixpos{5}(1):end_age_sess,2)==5);


for c=1:nb_event
%     if length(id_resp{c})>120
%         % Here to create an error
%         find(m_pos(1:end-1,2)==5 & (m_pos(2:end,2)==2 | m_pos(2:end,2)==4 | m_pos(2:end,2)==1 | m_pos(2:end,2)==3))
%         find(m_pos(1:end-1,2)==5 & ~(m_pos(2:end,2)~=23 | m_pos(2:end,2)~=24 | m_pos(2:end,2)~=22 | m_pos(2:end,2)~=20))
%     end
%     id_resp_conf{c}=id_resp{c}(2:2:end);
%     id_resp{c}=id_resp{c}(1:2:end);
    
    cue_latency_eeg{c}=m_pos(id_cue{c},1);%in samples... we need to get s_fs to convert samples into milliseconds: this info is in the .ent file
    fix_latency_eeg{c}=m_pos(id_fixpos{c},1);%
%     resp_latency_eeg{c}=m_pos(id_resp{c},1);
    
end


%%
nb_trial=60;
for c=1:nb_event
    [rapport{c}, rap{c}]=check_data_complete(id_fixpos{c},id_cue{c},0,nb_trial);
    disp(['Rating ' num2str(c) ': '])
    disp(rapport{c})
end
% for c=1:nb_event
%     disp(['CONF Rating' num2str(c) ': '])
%     
%     if length(id_resp_conf{c})~=nb_trial
%         disp(['Number of resp conf is not ' num2str(nb_trial) ' : ' num2str(length(id_resp_conf{c})) ', we try to correct that']);
%     else
%         disp(['Number of resp conf is ' num2str(nb_trial) ' : everything''s ok']);
%     end
%     if length(id_confcue{c})~=nb_trial
%         disp(['Number of cue conf is not ' num2str(nb_trial) ' : ' num2str(length(id_confcue{c})) ', we try to correct that']);
%     else
%         disp(['Number of cue conf is ' num2str(nb_trial) ' : everything''s ok']);
%     end
% end

pos_from_log=create_pos_from_log_main(a_code, v_ttime,[fix_latency_pcstim{1};fix_latency_pcstim{2};fix_latency_pcstim{3} ;fix_latency_pcstim{4} ;fix_latency_pcstim{5}],cue_latency_pcstim);

% find((m_pos(2:end,2)==51 & m_pos(1:end-1,2)~=5 & m_pos(1:end-1,2)~=50));

%% RAPPORT DE CORRECTION DU POS FILE

% [rapport1, rap]=check_data(id_fixpos,id_cue,id_resp);
% disp(rapport1)
old_pos=m_pos;
m_pos_corrected=m_pos;
m_pos=m_pos_corrected;
% Nettoyage d'artefacts
m_pos_corrected(~ismember(m_pos_corrected(:,2),[10 11 20 21 30 31 40 41 50 51 2 4 3 1 5 60 61 70 71 80 81]),:)=[];
% m_pos_corrected(find(m_pos_corrected(1:end-1,2)==5 & m_pos_corrected(2:end,2)~=50)+1,:)=[];

% %% CORRECTION FOR MISTAKE IN RESPONSE : look for doublon and delete the first one
%
% m_pos_corrected(find(m_pos_corrected(1:end-1,2)==5 & m_pos_corrected(2:end,2)==5),:)=[];

% Look for code between 5 and 50 and delete it
% m_pos_corrected(find(m_pos_corrected(1:end-2,2)==5 & m_pos_corrected(3:end,2)==50)+1,:)=[];

% find(id_resp_conf{1}(1:end-1)+2~=id_fixpos{1}(2:end))
% find(id_fixpos{3}(3:end-1)+1~=id_cue{3}(1:end))
for c=1:nb_event
    %% CUE
    
    if rap{c}(2)~=1
        %% MISSING CUES
        
        var_with_missing_trial=id_cue{c};
        time_ms_from_log=cue_rel_time_milliseconds{c};
        code=code_event_pos(c)+1;
        var_to_base_on=id_fixpos{c};
        latency_to_base_on=fix_latency_eeg{c};
        
        
        
        [m_pos_corrected, new_pos_line]=correct_from_pos_and_log(var_with_missing_trial,var_to_base_on,latency_to_base_on,...
            time_ms_from_log,s_fs,m_pos_corrected,code);
        
    end
end
m_pos=m_pos_corrected;

%% REDEFINE VAR



for c=1:nb_event
id_fixpos_ini{c}=find(m_pos(:,2)==code_event(c));
end


[~, order] = sort([id_fixpos_ini{1}(1) id_fixpos_ini{2}(1) id_fixpos_ini{3}(1)  id_fixpos_ini{4}(1)  id_fixpos_ini{5}(1)]);

code_event_pos=code_event(order);

for c=1:nb_event
id_fixpos{c}=id_fixpos_ini{order(c)};
id_cue{c}=find(m_pos(:,2)==code_event_pos(c)+1);%almost OK (119 events): we could of course use .log data to back-up the missing event;
% id_confcue{c}=find(m_pos(:,2)==code_event_pos(c)+3);
end

end_age_sess=find(m_pos(:,2)==60 | m_pos(:,2)==70 | m_pos(:,2)==80,1,'first');

% id_resp{1}=id_fixpos{1}(1)-1+find(m_pos(id_fixpos{1}(1):id_fixpos{2}(1),2)==5);
% id_resp{2}=id_fixpos{2}(1)-1+find(m_pos(id_fixpos{2}(1):id_fixpos{3}(1),2)==5);
% id_resp{3}=id_fixpos{3}(1)-1+find(m_pos(id_fixpos{3}(1):id_fixpos{4}(1),2)==5);
% id_resp{4}=id_fixpos{4}(1)-1+find(m_pos(id_fixpos{4}(1):id_fixpos{5}(1),2)==5);
% id_resp{5}=id_fixpos{5}(1)-1+find(m_pos(id_fixpos{5}(1):end_age_sess,2)==5);


for c=1:nb_event

%     id_resp_conf{c}=id_resp{c}(2:2:end);
%     id_resp{c}=id_resp{c}(1:2:end);
    
    cue_latency_eeg{c}=m_pos(id_cue{c},1);%in samples... we need to get s_fs to convert samples into milliseconds: this info is in the .ent file
    fix_latency_eeg{c}=m_pos(id_fixpos{c},1);%
%     resp_latency_eeg{c}=m_pos(id_resp{c},1);
    
end

%% FIX
for c=1:nb_event
    if rap{c}(1)~=1
        
        var_with_missing_trial=id_fixpos{c}+1;
        var_to_base_on=id_cue{c}-1;
        latency_to_base_on=cue_latency_eeg{c};
        time_ms_from_log=cue_rel_time_milliseconds{c};
        code=code_event_pos(c);
        
        [m_pos_corrected, new_pos_line]=correct_from_pos_and_log(var_with_missing_trial,var_to_base_on,latency_to_base_on,...
            time_ms_from_log,s_fs,m_pos_corrected,code);
        
    end
end

m_pos=m_pos_corrected;


%% REDEFINE VAR

for c=1:nb_event
id_fixpos_ini{c}=find(m_pos(:,2)==code_event(c));
end


[~, order] = sort([id_fixpos_ini{1}(1) id_fixpos_ini{2}(1) id_fixpos_ini{3}(1)  id_fixpos_ini{4}(1)  id_fixpos_ini{5}(1)]);

code_event_pos=code_event(order);

for c=1:nb_event
id_fixpos{c}=id_fixpos_ini{order(c)};
id_cue{c}=find(m_pos(:,2)==code_event_pos(c)+1);%almost OK (119 events): we could of course use .log data to back-up the missing event;
% id_confcue{c}=find(m_pos(:,2)==code_event_pos(c)+3);
end

end_age_sess=find(m_pos(:,2)==60 | m_pos(:,2)==70 | m_pos(:,2)==80,1,'first');

% id_resp{1}=id_fixpos{1}(1)-1+find(m_pos(id_fixpos{1}(1):id_fixpos{2}(1),2)==5);
% id_resp{2}=id_fixpos{2}(1)-1+find(m_pos(id_fixpos{2}(1):id_fixpos{3}(1),2)==5);
% id_resp{3}=id_fixpos{3}(1)-1+find(m_pos(id_fixpos{3}(1):id_fixpos{4}(1),2)==5);
% id_resp{4}=id_fixpos{4}(1)-1+find(m_pos(id_fixpos{4}(1):id_fixpos{5}(1),2)==5);
% id_resp{5}=id_fixpos{5}(1)-1+find(m_pos(id_fixpos{5}(1):end_age_sess,2)==5);
% 

for c=1:nb_event
%     id_resp_conf{c}=id_resp{c}(2:2:end);
%     id_resp{c}=id_resp{c}(1:2:end);
%     
    cue_latency_eeg{c}=m_pos(id_cue{c},1);%in samples... we need to get s_fs to convert samples into milliseconds: this info is in the .ent file
    fix_latency_eeg{c}=m_pos(id_fixpos{c},1);%
%     resp_latency_eeg{c}=m_pos(id_resp{c},1);    
end


%%
nb_trial=60;
for c=1:nb_event
    [rapport{c}, rap{c}]=check_data_complete(id_fixpos{c},id_cue{c},0,nb_trial);
    disp(['Rating ' num2str(c) ': '])
    disp(rapport{c})
end
% for c=1:nb_event
%     disp(['CONF Rating' num2str(c) ': '])
%     
%     if length(id_resp_conf{c})~=nb_trial
%         disp(['Number of resp conf is not ' num2str(nb_trial) ' : ' num2str(length(id_resp_conf{c})) ', we try to correct that']);
%     else
%         disp(['Number of resp conf is ' num2str(nb_trial) ' : everything''s ok']);
%     end
% 
% end

%% RESP
% for c=1:nb_event
% if rap{c}(3)~=1
%     error('not ready')
%     %% MISSING RESPONSE
%     id_trial_no_resp=[];
%     %let's correct that!
%     if length(id_resp)<120
%         
%         if length(id_fixpos)>=120 %if we have all cue: it's quite easy:
%             %             id_trial_no_resp=find(~ismember(id_fixpos(2:end)-1,id_resp));%identify trial number for which we did not have cue and pick up from logfile the relative latency of the stim event + resp event during this trial
%             id_trial_no_resp=find(m_pos_corrected(id_fixpos(2:end)-1,2)~=5);
%             resp_latency_tmp=[fix_latency_eeg(id_trial_no_resp(1:end)+1)];% fix_latency_eeg(id_trial_no_resp(1:end)+2)];
%             resp_rel_time_milliseconds_tmp=[resp_rel_time_milliseconds(id_trial_no_resp(1:end))];%  resp_rel_time_milliseconds(id_trial_no_resp(1:end)+1) ];
%             %we need that info in EEG samples: use s_fs now
%             resp_rel_time_tmp=round(resp_rel_time_milliseconds_tmp*s_fs/1000);
%             
%             resp_onsets=resp_latency_tmp-200;
%             
%             new_pos_line=[resp_onsets(:,1) 5*ones(length(resp_onsets(:,1)),1) zeros(length(resp_onsets(:,1)),1) 1000*ones(length(resp_onsets(:,1)),1)];
%             m_pos_corrected=[m_pos_corrected zeros(size(m_pos_corrected,1),1) ;new_pos_line];
%             %now re-order m_pos with respect to timings of column 1
%             [a,b]=sort(m_pos_corrected(:,1));
%             m_pos_corrected=m_pos_corrected(b,:);%done: we now have all events of interest in the posfile :)
%             %                end
%         end
%         
%         id_fixpos=find(m_pos_corrected(:,2)==50);%OK: 120 events
%         id_cue=find( m_pos_corrected(:,2)==51);%almost OK (119 events): we could of course use .log data to back-up the missing event;
%         id_resp=find( m_pos_corrected(:,2)==5);
%     else
%         m_pos_corrected=[m_pos_corrected zeros(size(m_pos_corrected,1),1)];
%         
%     end
%     % else
%     %     m_pos_corrected=m_pos;
%     m_pos_corrected=[m_pos_corrected zeros(size(m_pos_corrected,1),1)];
%     % MARK INFERED RESPONSE TIME
%     %  m_pos_corrected(m_pos_corrected(:,4)==1000,1)=NaN(1,sum(m_pos_corrected(:,4)==1000));
%     
% end
% end


%% CONTROL FOR CORRESPONDANCE BETWEEN POS AND LOG files (reaction times)
% give trial number to the pos file
code_event_real=[code_event_real m_pos_corrected(end_age_sess,2)];
end_age_sess=find(m_pos_corrected(:,2)==60 | m_pos_corrected(:,2)==70 | m_pos_corrected(:,2)==80,1,'first');
pos_from_log(end,2)=m_pos_corrected(end_age_sess,2);
code_event_pos=[code_event_pos m_pos_corrected(end_age_sess,2)];
    m_pos_corrected(:,5)=zeros(size(m_pos_corrected,1),1);
    pos_from_log(1:end-1,5)=zeros(size(pos_from_log,1)-1,1);
    trial=1;
for c=1:nb_event

index_base=code_event_pos(c);
    
m_pos_corrected(find(m_pos_corrected(:,2)==index_base),5)=1+(nb_trial*(c-1)):nb_trial*c;
 for i=find((m_pos_corrected(:,2)==code_event_pos(c)),1,'first')+1:find((m_pos_corrected(:,2)==code_event_pos(c+1)),1,'first')
        if m_pos_corrected(i,5)~=trial+1
            m_pos_corrected(i-1,5)=trial;
        else
            m_pos_corrected(i-1,5)=trial;
            trial=trial+1;
        end
 end
 trial=trial+1;
 
 pos_from_log(find(pos_from_log(:,2)==index_base),5)=1+(nb_trial*(c-1)):nb_trial*c;
 for i=find((pos_from_log(:,2)==code_event_pos(c)),1,'first')+1:find((pos_from_log(:,2)==code_event_pos(c+1)),1,'first')
        if pos_from_log(i,5)==0
                pos_from_log(i,5)=pos_from_log(i-1,5);
        end
 end
     
    
end


%% LAST CHECK: GET RT_pos and RT_log
% Remise en ordre des essais (entre pos et log)

for c=1:nb_event

for trial=(1+(c-1)*nb_trial):nb_trial*c-1
    RT_pos(trial)=(m_pos_corrected(find(m_pos_corrected(:,5)==trial+1,1,'first'),1)-m_pos_corrected(find(m_pos_corrected(:,5)==trial,1,'first'),1))/s_fs;
    RT_log(trial)=(pos_from_log(find(pos_from_log(:,5)==trial+1,1,'first'),1)-pos_from_log(find(pos_from_log(:,5)==trial,1,'first'),1))/1000;
end

first_code=find((m_pos_corrected(:,2)==code_event_pos(c)),1,'first');
first_next_code=find((m_pos_corrected(:,2)==code_event_pos(c+1)),1,'first');

code=code_event_pos(c);
last_answer=first_code-1+find(m_pos_corrected(first_code:first_next_code,2)==5,2,'last'); 
last_stim=first_code-1+find(m_pos_corrected(first_code:first_next_code,2)==code,1,'last');

RT_pos(nb_trial*c)=(m_pos_corrected(last_answer(1),1)-m_pos_corrected(last_stim,1))/s_fs;

first_code_log=find((pos_from_log(:,2)==code_event_pos(c)),1,'first');
first_next_code_log=find((pos_from_log(:,2)==code_event_real(find(code_event_real==code_event_pos(c))+1)),1,'first');

last_answer_log=first_code_log-1+find(pos_from_log(first_code_log:first_next_code_log,2)==5,2,'last'); 
last_stim_log=first_code_log-1+find(pos_from_log(first_code_log:first_next_code_log,2)==code,1,'last');

RT_log(nb_trial*c)=(pos_from_log(last_answer_log(1),1)-pos_from_log(last_stim_log,1))/1000;


end

outliers=find(abs(round(RT_pos*10)-round(RT_log*10))>2);

if length(outliers)>15
    figure;scatter(RT_log(1:end),RT_pos(1:end))
end
% figure;scatter(RT_log,RT_pos)
[betas stat]=robustfit(RT_log(1:end),RT_pos(1:end));
if round(betas(2))==1 && round(betas(1))==0 && stat.p(1)>0.05
    display([num2str(length(outliers)) ' outliers, removed from pos file'])
else
    weird=1;
    figure;scatter(RT_log(1:end),RT_pos(1:end))

    bar([RT_log(outliers) ; RT_pos(outliers)]')
    error('Timing is unreliable between pos and log, careful with this subject, consider correct the pos file manually')
end



m_pos_corrected(ismember(m_pos_corrected(:,5),outliers),:)=[];
trials_included=~ismember(1:nb_trial*nb_event,outliers);
RT_pos=RT_pos(trials_included);
id_cue=find(m_pos_corrected(:,2)==11 |m_pos_corrected(:,2)==21 | m_pos_corrected(:,2)==31 | m_pos_corrected(:,2)==41  | m_pos_corrected(:,2)==51);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 28/05/2015 19:57
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%NEXT: ADD INFO IN THE POSFILE to have a quick acess to different REGRESSORS

clear m_pos
m_pos=m_pos_corrected(:,1:2);

%NEXT: ADD INFO IN THE POSFILE to have a quick acess to different REGRESSORS
index_name_sub=find(ismember(a_oldposfile,'\'),5,'last');
name_sub=a_oldposfile((index_name_sub(1)+1):index_name_sub(2)-1);
load('E:\ALIZEE\EPILEPSY\DATA_PREPROCESS\Behavior_complete\all_sub.mat')
load('E:\ALIZEE\EPILEPSY\Analysis\Behavior\Extraction\Sub_names.mat')
index_behavior=find(strcmp([Sub_Names],name_sub));

%% AGE


if id_fixpos{2}(1)>id_fixpos{1}(1)
index_behavior=find(strcmp([Sub_Names],name_sub));
data_age1_sub=Data_trials{1}{1}{index_behavior};
age1=data_age1_sub.Age(:);
conf1=data_age1_sub.Age_conf(:);
data_age2_sub=Data_trials{1}{2}{index_behavior};
age2=data_age2_sub.Age(:);
conf2=data_age2_sub.Age_conf(:);
else
 index_behavior=find(strcmp([Sub_Names],name_sub));
data_age1_sub=Data_trials{1}{1}{index_behavior};
age1=data_age1_sub.Age(:);
conf1=data_age1_sub.Age_conf(:);
data_age2_sub=Data_trials{1}{2}{index_behavior};
age2=data_age2_sub.Age(:);
conf2=data_age2_sub.Age_conf(:);
end


data_pref1_sub=Data_trials{2}{1}{index_behavior};
index1=data_pref1_sub.Pleas_ID;

for t=1:60
    pref1(t)=data_pref1_sub.Pleas(find(data_age1_sub.Age_ID(t)==index1));
end
data_pref2_sub=Data_trials{2}{2}{index_behavior};
index2=data_pref2_sub.Pleas_ID;

for t=1:60
    pref2(t)=data_pref2_sub.Pleas(find(data_age2_sub.Age_ID(t)==index2));
end

if id_fixpos{2}(1)>id_fixpos{1}(1)% 20 after 10
ageA=[zscore(age1);zscore(age2)];
confA=[conf1;conf2];
prefA=[pref1';pref2'];
else
ageA=[zscore(age2);zscore(age1)];
confA=[conf2;conf1];
prefA=[pref2';pref1'];   
end


%% PREF

prefP=[];
confP=[];
[~,order_pos]=sort(code_event_pos(3:5));
[~,order2]=sort(order_pos);
for c=1:3
data_pref_sub=Data_trials{2}{order2(c)}{index_behavior};
prefP=[prefP;data_pref_sub.Pleas];
confP=[confP; data_pref_sub.Pleas_conf];
end


pref=[prefA; prefP];
conf=[confA; confP];

pref=pref(trials_included);
conf=conf(trials_included);

m_pos_tmp=[];
% m_pos_tmp=m_pos(id_cue(11:end),:);
m_pos_tmp=m_pos(id_cue,:);
m_pos_tmp(:,3)=pref'; %Ok this is one of the regressor we need... BUT NOT THE ONLY ONE !!
m_pos_tmp=[ m_pos_tmp conf];

RT=RT_pos;
m_pos_tmp=[ m_pos_tmp RT']; %en colonne 5: la TR comme régresseur 3


%ce qui serait très bien c'est de modifier aussi le .pos pour avoir les
%median split (ajouter un event code pour les essais où la valeur est >median


%MEDIAN SPLIT EVENT CODES: create events 52 for low value trials and events
%53 for high value trials using the median of pref_score///
id_minpref=find(pref<median(pref));
id_maxpref=find(pref>median(pref));
median_split_trials=m_pos_tmp(id_minpref,:);
median_split_trials(:,2)=100;
tmp=m_pos_tmp(id_maxpref,:);tmp(:,2)=101;
median_split_trials=[median_split_trials; tmp];
m_pos_tmp=[m_pos_tmp ;median_split_trials];
%now re-order m_pos with respect to timings of column 1
[a,b]=sort(m_pos_tmp(:,1));
m_pos_tmp=m_pos_tmp(b,:);%done: we now have all events of interest in the posfile :)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nb_col=[];
for col=1:size(m_pos_tmp,2)-1
    nb_col=[nb_col '%d\t'];
end
nb_col=[nb_col '%d\n'];


f=fopen(a_newposfile,'w');
for i=1:length(m_pos_tmp)
    fprintf(f,nb_col,m_pos_tmp(i,:));
end


fclose(f);





