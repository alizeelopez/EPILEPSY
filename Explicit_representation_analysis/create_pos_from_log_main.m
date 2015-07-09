function m_pos_from_log=create_pos_from_log_main(a_code, v_ttime,fix_latency_pcstim,cue_latency_pcstim)

i=0;
j=0;
for ii=1:length(a_code)
    %% RATING
    m_pos_from_log(i+1,1)=0;
    if strcmp(a_code{ii},'10') || strcmp(a_code{ii},'20') || strcmp(a_code{ii},'30') || strcmp(a_code{ii},'40') ||  strcmp(a_code{ii},'50')
        if ~strcmp(a_code{ii+1},'10') && ~strcmp(a_code{ii+1},'20') && ~strcmp(a_code{ii+1},'30') && ~strcmp(a_code{ii+1},'40') &&  ~strcmp(a_code{ii+1},'50')
            if length(a_code{ii+1})>0  && length(a_code{ii+1})<6 && ( strcmp(a_code{ii+1}(1),'T') || strcmp(a_code{ii+1}(1),'V')) || strcmp(a_code{ii+1}(1),'N')
            i=i+1;
            j=j+1;
            m_pos_from_log(i,2)=str2double(a_code{ii});
            m_pos_from_log(i,1)=fix_latency_pcstim(j);
            
            last_code=str2double(a_code{ii});
            if str2double(a_code{ii})==30
                code_30(j,:)=[i ii];
            elseif str2double(a_code{ii})==40
                code_40(j,:)=[i ii];
            elseif str2double(a_code{ii})==50
                code_50(j,:)=[i ii];
            end
            end
        end
    elseif length(a_code{ii})>0  && length(a_code{ii})<6 && ( strcmp(a_code{ii}(1),'T') || strcmp(a_code{ii}(1),'V')) || strcmp(a_code{ii}(1),'N')
        i=i+1;
        
        
        m_pos_from_log(i,2)=last_code+1;
        m_pos_from_log(i,1)=cue_latency_pcstim(j);
    elseif strcmp('5',a_code{ii})
        i=i+1;
        
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
    m_pos_from_log(i,6)=j;
end

% [a,b]=sort(m_pos_from_log(:,1));
% m_pos_from_log=m_pos_from_log(b,:);%done: we now have all events of interest in the posfile :)

