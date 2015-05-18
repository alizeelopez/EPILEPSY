function loc2_write_conf(a_conffile,v_filter,m_bip)

% this function will create a standard conf file, from a v_filter vector
% (length = nbchannel), that indicates for all the channel their visibility
% (1 or 0), and indicate the reference for each channel using m_bip, the
% matrix that indicates all the electrode pairs that can be formed
% note, v_filter is supposed to include the technical channels

s_nbchannel=length(v_filter);

f_conf=fopen(a_conffile,'w'); % this is the txt file that contains the hdr in the ent file

fprintf(f_conf,'%s\t%s\n','nb_channel',int2str(s_nbchannel));
fprintf(f_conf,'%s\t%s\n','sec_per_page',int2str(4));
fprintf(f_conf,'%s\t%s\n','amp_scale_type',int2str(1));
fprintf(f_conf,'%s\t%s\n','amp_scale_val',int2str(1));
fprintf(f_conf,'%s\n','channel_visibility'); %  JP ?????

for s_c=1:s_nbchannel
    s_line=v_filter(s_c);
    a_line=num2str(s_line);
    fprintf(f_conf,'%s\n',a_line); % 
end; % for s_c


fprintf(f_conf,'%s\n','channel_reference'); %  JP ?????

for s_c=1:s_nbchannel
    % at this point we have to check whether this elec is a big brother or not
    % since s_c here is a relative indice, we can look directly in the bipole matrix
    % if it is a bigbrother, then we can directly put the relative indice of the little brother (-1, that's the way it is)
    v_i=find(m_bip(:,1)==s_c);
    if (~isempty(v_i))
        s_little=max(m_bip(v_i(1),2)-1,-1); % yes, because values lower than -1 are not allowed
    else
        s_little=-1;
    end; % if 
    a_line=num2str(s_little);
    fprintf(f_conf,'%s\n',a_line); % 
end; % for s_c

fclose(f_conf);




