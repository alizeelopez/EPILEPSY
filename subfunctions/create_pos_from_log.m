function m_pos_from_log=create_pos_from_log(a_code, v_ttime,fix_latency_pcstim,cue_latency_pcstim)

i=0;
j=0;
for ii=1:length(a_code)
    %% RATING
    m_pos_from_log(i+1,1)=0;
    if strcmp(a_code{ii},'50')
        i=i+1;
        j=j+1;
        m_pos_from_log(i,2)=50;
        m_pos_from_log(i,1)=fix_latency_pcstim(j);
        if ~ismember('N',a_code{ii+2})
            if ~ismember('N',a_code{ii+1})
            display('Image ID is missing in the log file')
            i=i+1;
            m_pos_from_log(i,2)=51;
            m_pos_from_log(i,1)=cue_latency_pcstim(j);
            end
        end
    elseif length(a_code{ii})>0 && strcmp(a_code{ii}(1),'N') && length(a_code{ii})<6
        i=i+1;
        m_pos_from_log(i,2)=51;
        m_pos_from_log(i,1)=cue_latency_pcstim(j);
    elseif strcmp('5',a_code{ii})
        i=i+1;
        
        m_pos_from_log(i,2)=5;
        m_pos_from_log(i,3)=str2double(a_code{ii+1}(5:end));
        pref_from_log(j,1)=str2double(a_code{ii+1}(5:end));
        
    end   
    %% CHOICES
    
    m_pos_from_log(i+1,1)=0;
    if strcmp(a_code{ii},'80')
        i=i+1;
        j=j+1;
        m_pos_from_log(i,2)=80;
        m_pos_from_log(i,1)=fix_latency_pcstim(j);
        if ~ismember('N',a_code{ii+2})
            if ~ismember('N',a_code{ii+1})
            display('Image ID is missing in the log file')
            i=i+1;
            m_pos_from_log(i,2)=81;
            m_pos_from_log(i,1)=cue_latency_pcstim(j);
            end
        end
    elseif length(a_code{ii})>0 && strcmp(a_code{ii}(1),'N') && length(a_code{ii})>5
        i=i+1;
        m_pos_from_log(i,2)=81;
        m_pos_from_log(i,1)=cue_latency_pcstim(j);

        
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
    m_pos_from_log(i,6)=j;
end

