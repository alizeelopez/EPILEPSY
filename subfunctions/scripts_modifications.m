%% CHANGEMENTS APORTES AUX SCRIPTS

%extract_LOCA_pref_datafromlogfile :

%Line 78-93 :
% Initial :
try
    if uch(1)=='N'      %stim IF begins by letter 'N'
        stim_id(s_count)=str2double(uch(2:end));
        s_count=s_count+1;
        tmp=[tmp s_i];
    end
catch
end
try
    if strcmp(uch(1:4),'CF1=')
        k=k+1;
        pref_score(k)=str2num(uch(5:end));
    end
catch
end
%             Replaced by :

if ~isempty(uch)
    if uch(1)=='N'      %stim IF begins by letter 'N'
        stim_id(s_count)=str2double(uch(2:end));
        s_count=s_count+1;
        tmp=[tmp s_i];
    end
    if length(uch)>=4
        
        if strcmp(uch(1:4),'CF1=')
            k=k+1;
            pref_score(k)=str2num(uch(5:end));
        end
    end
end

%% SCRIPT chg_pos_locapref, all id_fix for pos are renamed id_fixpos

%% SCRIPT chg_pos_locapref, lines 106-127
for i=1:length(id_cue)
    resp=[];
    try
    resp=find(m_pos_corrected(id_cue(i):id_cue(i+1),2)==5);
    catch
        if i==length(id_cue)
            resp=find(m_pos_corrected(id_cue(i):end,2)==5);
        end
    end
    if isempty(resp)
        %a possibility is that instead of event 5  you actually get a
        %prefscore measure (I had to add +10 to the actual prefscore so
        %that pref_socre =0 corresponds to event code 11
        resp=find(m_pos_corrected(id_cue(i):id_cue(i+1),2)==50)-1;
    end
    
    resp_lat=m_pos_corrected(id_cue(i)-1+resp(1),1);
    RT(i)=resp_lat-m_pos_corrected(id_cue(i),1); %TR en sample = OK comme régresseur 3
end
    % replaced by
    
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
    
    if m_pos_corrected(id_cue(i))==m_pos_corrected(id_cue(i+1)) %% One trial is missing
    resp_lat=NaN;
    RT(i)=NaN; %TR en sample = OK comme régresseur 3
    else        
    resp_lat=m_pos_corrected(id_cue(i)-1+resp(1),1);
    RT(i)=resp_lat-m_pos_corrected(id_cue(i),1); %TR en sample = OK comme régresseur 3
    end
end


%% eeg_loc_montage line 115 to 118

    clear v_rank;
    for s_i=1:length(v_id)
        v_rank(s_i)=find(v_y==v_id(s_i))+s_total;
    end
    
    replaced by 
    
        clear v_rank;
    for s_i=1:length(v_id)
        v_rank(s_i)=find(v_y==v_id(s_i),1,'first')+s_total;
    end
