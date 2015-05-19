
function eeg2eeg_write_eeg_init(a_eegfile,s_nbsample,s_nbchannel)

% INITIALIZE THE EEG FILE !

f_eeg=fopen(a_eegfile,'w','ieee-be'); % this is the file with the data

for s_n = 1:s_nbsample
    s_nbsampwritten=fwrite(f_eeg,zeros(1,s_nbchannel),'int16');
end;
fclose(f_eeg);
