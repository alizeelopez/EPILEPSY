function eeg2eeg_write_eeg_onechannel(m_data_sc,a_eegfile,s_channel,s_nbchannel)

% CREATE THE EEG FILE !

f_eeg=fopen(a_eegfile,'r+','ieee-be'); % this is the file with the data
m_data=m_data_sc';
s_datasizebyte=2;
fseek(f_eeg,(s_channel)*s_datasizebyte,-1); % This was the problem
s_nbsampwritten=fwrite(f_eeg,m_data(2:length(m_data)),'int16',s_datasizebyte*(s_nbchannel-1)); % we forget the first point, because of the skip field, which 'skips' 'write' 'skips' 'write' etc.
fclose(f_eeg);

