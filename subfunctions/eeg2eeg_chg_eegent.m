
function eeg2eeg_chg_eegent(a_oldent,a_newent,ss_bipole,s_newfs)



s_nbip = length(ss_bipole);

% then read the eeg.ent file
f_old=fopen(a_oldent,'r');
f_new=fopen(a_newent,'w');

for s_i=1:8
    a_line=fgetl(f_old);
    fprintf(f_new,'%s\n',a_line);
end;

% line 9 is the s_fs
a_line=fgetl(f_old);
fprintf(f_new,'%s\n',num2str(1/s_newfs));

% number of channels including the techinical ones (last two)
s_nbchannel = s_nbip+2;
a_nbchannel=int2str(s_nbchannel);
fprintf(f_new,'%s\n',a_nbchannel);

% write down the electrode names
for s_c=1:s_nbip
    % what is the absolute indice of the s_c-th visible channel?
    a_elecname=[deblank(ss_bipole(s_c).name) deblank(ss_bipole(s_c).afterpoint) ];   % JP?
    fprintf(f_new,'%s\n',a_elecname);
end; % for s_c
fprintf(f_new,'%s\n',['Num1']); % JP ?
fprintf(f_new,'%s\n',['Num2']); % JP ?

% what's the type of each channel? we set it to EEG by default
for s_c=1:s_nbip
    fprintf(f_new,'%s\n','Electrode EEG'); % JP ?
end; % for s_c
fprintf(f_new,'%s\n','dateur'); % JP ?
fprintf(f_new,'%s\n','event'); % JP ?

% what's the unit for each channel?
for s_c=1:s_nbip
    a_unit='microV'; % forced to microV for Grenoble's data
    fprintf(f_new,'%s\n',a_unit(1:min(6,length(a_unit)))); % JP ?
end; % for s_c
fprintf(f_new,'%s\n','sans'); % JP ?
fprintf(f_new,'%s\n','sans'); % JP ?

% what's the physical minimum for each channel?
for s_c=1:s_nbip
    s_line=-3200;
    a_line=num2str(s_line);
    fprintf(f_new,'%s\n',a_line); %
end; % for s_c
fprintf(f_new,'%s\n','-1'); %
fprintf(f_new,'%s\n','-1'); %

% what's the physical maximum for each channel?
for s_c=1:s_nbip
    s_line=3200;
    a_line=num2str(s_line);
    fprintf(f_new,'%s\n',a_line); %
end; % for s_c
fprintf(f_new,'%s\n','1'); %
fprintf(f_new,'%s\n','1'); %

% what's the digital minimum for each channel?
for s_c=1:s_nbip
    s_line=-32768;
    a_line=num2str(s_line);
    fprintf(f_new,'%s\n',a_line); %
end; % for s_c
fprintf(f_new,'%s\n','-1'); %
fprintf(f_new,'%s\n','-1'); %

% what's the digital maximum for each channel?
for s_c=1:s_nbip
    s_line=32768;
    a_line=num2str(s_line);
    fprintf(f_new,'%s\n',a_line); %
end; % for s_c
fprintf(f_new,'%s\n','1'); %  JP
fprintf(f_new,'%s\n','1'); %  JP


for s_c=1:s_nbip
    fprintf(f_new,'%s\n','filter'); %
end; % for s_c
fprintf(f_new,'%s\n','sans'); %  JP
fprintf(f_new,'%s\n','sans'); %  JP

for s_c=1:s_nbip+2
    fprintf(f_new,'%s\n','1'); %
end; % for s_c

for s_c=1:s_nbip+2
    fprintf(f_new,'%s\n','reserved'); %
end; % for s_c



fclose(f_new);
fclose(f_old);
