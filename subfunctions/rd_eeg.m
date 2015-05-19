function m_data_sc = rd_eeg(a_eegfilename,v_channumber,s_nbchannel)


f_in=fopen(a_eegfilename,'r','ieee-be');

s_datasizebyte=2;

% read the data
m_data_sc=[];
for s_c=1:length(v_channumber)
    s_chan_marker = v_channumber(s_c);
    fseek(f_in,(s_chan_marker-1)*s_datasizebyte,-1);
    v_data=fread(f_in,[1,inf],'int16',s_datasizebyte*(s_nbchannel-1));
    m_data_sc = [m_data_sc v_data(:)];    
end;
fclose(f_in);