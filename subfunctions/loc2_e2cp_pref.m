function loc2_e2cp_pref(a_rootname, a_elecdat)

% this function takes as input a file a_rootname.eeg.ent, and creates a
% a_rootname.conf and a a_rootname_bipo.par and a a_rootname_bipo.par with
% the correct references

% the first thing is to get the names of the electrodes


% first adjust the elec.dat
chg_eegent(a_rootname,a_elecdat);


f_old=fopen([a_rootname '.eeg.ent'],'r');

for s_i=1:9
    a_line=fgetl(f_old);
end;

s_fs = 1/str2num(a_line);
s_fs = round(s_fs);
a_line=fgetl(f_old);

clear ss_channel;
% how many channels ?
s_nbchannel=str2num(a_line);
for s_c=1:s_nbchannel-2
    % read the name of the actual data channels
    s_indice=-1;
    a_line=fgetl(f_old);
    % the name is something like Cz.10;
    s_point=find(a_line=='.');
    if (~isempty(s_point))
        a_elecname=a_line(1:s_point(1)-1);
    end;
    a_elecname=upper(a_elecname);
    a_elecname=strrep(a_elecname,' ','');
    % now we have to find it in the list of elec.dat
    ss_channel.values(s_c).value=a_elecname;    
end; % for s_c

[v_neworder, v_neworder_rev, m_bipole]=loc2_montage(ss_channel);
fclose(f_old);

loc2_write_conf([a_rootname '.conf'],ones(s_nbchannel,1),[m_bipole; [s_nbchannel-1 -1;s_nbchannel -1]]);


loc2_write_par_visu_Wil(a_rootname,s_fs,s_nbchannel-2,m_bipole,[52 53 82 83]);%il faudrait retravailler un .pos pour avoir un TF pour les essais faible Valeur vs. forte valeur éventuellement ?

a_batch = [a_rootname '_batch.txt']
f_batch=fopen(a_batch,'w');
fprintf(f_batch,'%s\t%s\n','',['tfstat ' a_rootname '.eeg ' a_rootname '.pos ' a_rootname '_bipo.par']);
fclose(f_batch);









function loc2_write_par_visu_Wil(a_rootname,s_fs,s_nbchannel_notec,m_bipole,a_listevent)

% s_flag = 1 for kruskal
% s_flag = 0 for wilcoxon or tfavg

% at this point, we have the .conf, the .eeg.ent, we still need the .par for the mono and bipo and the eeg, and the pos
% FOR THE .PAR FILES
% what are the info that we need ?
% we need : the info about the events that are contained in the file
% that can be adapted after, may one should just put a 1, here, and adjust after
% we need : latencies in samples of the baseline for the evoked potential
% that should be in the dim(1), why not take that as a baseline ?
% for tomorrow, that's OK, I will adjust that. just write the bipoles
% mostly, I need to identify the number of bipoles that I want.
% for each of the channels, reordered, I need to say whether it is a 'big brother'
% or not. if it is not, then I set it to 0, and the ref also, if it is a big brother
% then I set it to 1, and I specify the little brother as the ref.

s_nbchannel = s_nbchannel_notec;

s_nombre_evt = length(a_listevent);
%a_listevent = (1:s_nombre_evt);

a_parfile_bipo=[a_rootname '_bipo.par'];
a_parfile_mono=[a_rootname '_mono.par'];

[a_path,a_thispre,a_ext,a_ver] = fileparts(a_rootname);
s_nbevtype = length(a_listevent);
v_evt=a_listevent;
v_evt=v_evt(:);

% start with the bipo
f_par=fopen(a_parfile_bipo,'w');

% TO BE ADJUSTED TO THE EXPERIMENT
fprintf(f_par,'%s\t%s\n','fileprefix',[a_thispre '_bipo']); %  JP ?????
fprintf(f_par,'%s\t%s\n','nb_eventcode',int2str(s_nbevtype)); %  JP ?????
fprintf(f_par,'%s\t%s\n','list_eventcode',int2str(v_evt')); %  JP ?????
fprintf(f_par,'%s\t%s\n','prestim_nbsample',int2str(s_fs*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','poststim_nbsample',int2str(s_fs*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','baseline_eventcode',int2str(v_evt')); %  JP ?????
fprintf(f_par,'%s\t','baseline_msec_start');
for s_i = 1:s_nbevtype
    fprintf(f_par,'%s\t','-300');
end;
fprintf(f_par,'%s\n','');

fprintf(f_par,'%s\t','baseline_msec_start');
for s_i = 1:s_nbevtype
    fprintf(f_par,'%s\t','-100');
end;
fprintf(f_par,'%s\n','');

% read the file with the bipoles.
% 1. read the excel file with all the conditions
v_bigbro=m_bipole(:,1);
v_litbro=m_bipole(:,2);

% v_bigpro contains the relative indices, after reordering, of all the big brothers
v_chan_flag=zeros(1,s_nbchannel+2);
v_chan_ref=zeros(1,s_nbchannel+2);
v_chan_flag(v_bigbro)=1; % all the big brothers are visible in this par file
v_chan_ref(v_bigbro)=v_litbro; % the corresponding little brother % only thing to change for the monopolar

fprintf(f_par,'%s\t','ep_channel_flag');
for s_ii=1:length(v_chan_flag)
    fprintf(f_par,'%s',int2str(v_chan_flag(s_ii)));
    if (rem(s_ii,10)==0)
        fprintf(f_par,'%s\n\t\t\t\t\t\t','');
    else
        fprintf(f_par,'%s\t','');
    end;
end;
fprintf(f_par,'%s\n\n','');

fprintf(f_par,'%s\t','ep_channel_rejtype');
for s_ii=1:length(v_chan_flag)
    fprintf(f_par,'%s',int2str(0));
    if (rem(s_ii,10)==0)
        fprintf(f_par,'%s\n\t\t\t\t\t\t','');
    else
        fprintf(f_par,'%s\t','');
    end;
end;
fprintf(f_par,'%s\n\n','');

fprintf(f_par,'%s\t','ep_channel_ref');
for s_ii=1:length(v_chan_flag)
    fprintf(f_par,'%s',int2str(max(0,v_chan_ref(s_ii))));
    if (rem(s_ii,10)==0)
        fprintf(f_par,'%s\n\t\t\t\t\t\t','');
    else
        fprintf(f_par,'%s\t','');
    end;
end;
fprintf(f_par,'%s\n\n','');

fprintf(f_par,'%s\t%s\n\n','erpa_positivity_up','1'); %  JP ?????

fprintf(f_par,'%s\t','tf_channel_flag');
for s_ii=1:length(v_chan_flag)
    fprintf(f_par,'%s',int2str(v_chan_flag(s_ii)));
    if (rem(s_ii,10)==0)
        fprintf(f_par,'%s\n\t\t\t\t\t\t','');
    else
        fprintf(f_par,'%s\t','');
    end;
end;
fprintf(f_par,'%s\n\n','');

fprintf(f_par,'%s\t','tf_channel_ref');
for s_ii=1:length(v_chan_flag)
    fprintf(f_par,'%s',int2str(max(v_chan_ref(s_ii),0)));
    if (rem(s_ii,10)==0)
        fprintf(f_par,'%s\n\t\t\t\t\t\t','');
    else
        fprintf(f_par,'%s\t','');
    end;
end;
fprintf(f_par,'%s\n\n','');

% to copy
fprintf(f_par,'%s\t%s\n','eegstat_type_stat',int2str(1*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','eegstat_time_hw',int2str(1*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','eegstat_time_step',int2str(1*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','eegstat_baseline_flag',int2str(0*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','eegstat_baseline_start',int2str(-300*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','eegstat_baseline_stop',num2str(-100*ones(1,s_nbevtype))); %  JP ?????

fprintf(f_par,'%s\t%s\n','tfstat_type_stat',int2str(1*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','tfstat_baseline_flag',int2str(0*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','tfstat_freq_hw',int2str(3*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','tfstat_freq_step',int2str(3*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','tfstat_time_hw',int2str(50*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','tfstat_time_step',int2str(50*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','tfstat_baseline_start',int2str(-300*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','tfstat_baseline_stop',num2str(-100*ones(1,s_nbevtype))); %  JP ?????

fprintf(f_par,'%s\t%s\n','tf_freq_start',int2str(1*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','tf_freq_stop',int2str(200*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','tf_freq_step',int2str(3*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','tf_nb_sample_blackman',int2str(30*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','tf_wavelet_type',int2str(1*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','tf_morlet_m',int2str(7*ones(1,s_nbevtype))); %  JP ?????

fclose(f_par);


f_par=fopen(a_parfile_mono,'w');

fprintf(f_par,'%s\t%s\n','fileprefix',[a_thispre '_mono']); %  JP ?????
fprintf(f_par,'%s\t%s\n','nb_eventcode',int2str(s_nbevtype)); %  JP ?????
fprintf(f_par,'%s\t%s\n','list_eventcode',int2str(v_evt')); %  JP ?????
fprintf(f_par,'%s\t%s\n','prestim_nbsample',int2str(s_fs*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','poststim_nbsample',int2str(s_fs*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','baseline_eventcode',int2str(v_evt')); %  JP ?????
fprintf(f_par,'%s\t','baseline_msec_start');
for s_i = 1:s_nbevtype
    fprintf(f_par,'%s\t','-300');
end;
fprintf(f_par,'%s\n','');

fprintf(f_par,'%s\t','baseline_msec_start');
for s_i = 1:s_nbevtype
    fprintf(f_par,'%s\t','-100');
end;
fprintf(f_par,'%s\n','');

% read the file with the bipoles.
% 1. read the excel file with all the conditions
v_bigbro=m_bipole(:,1);
v_litbro=m_bipole(:,2);

% v_bigpro contains the relative indices, after reordering, of all the big brothers
v_chan_flag=zeros(1,s_nbchannel+2);
v_chan_ref=zeros(1,s_nbchannel+2);
v_chan_flag(v_bigbro)=1; % all the big brothers are visible in this par file
%v_chan_ref(v_bigbro)=v_litbro; % the corresponding little brother % only thing to change for the monopolar

fprintf(f_par,'%s\t','ep_channel_flag');
for s_ii=1:length(v_chan_flag)
    fprintf(f_par,'%s',int2str(v_chan_flag(s_ii)));
    if (rem(s_ii,10)==0)
        fprintf(f_par,'%s\n\t\t\t\t\t\t','');
    else
        fprintf(f_par,'%s\t','');
    end;
end;
fprintf(f_par,'%s\n\n','');

fprintf(f_par,'%s\t','ep_channel_rejtype');
for s_ii=1:length(v_chan_flag)
    fprintf(f_par,'%s',int2str(0));
    if (rem(s_ii,10)==0)
        fprintf(f_par,'%s\n\t\t\t\t\t\t','');
    else
        fprintf(f_par,'%s\t','');
    end;
end;
fprintf(f_par,'%s\n\n','');

fprintf(f_par,'%s\t','ep_channel_ref');
for s_ii=1:length(v_chan_flag)
    fprintf(f_par,'%s',int2str(v_chan_ref(s_ii)));
    if (rem(s_ii,10)==0)
        fprintf(f_par,'%s\n\t\t\t\t\t\t','');
    else
        fprintf(f_par,'%s\t','');
    end;
end;
fprintf(f_par,'%s\n\n','');

fprintf(f_par,'%s\t%s\n\n','erpa_positivity_up','1'); %  JP ?????

fprintf(f_par,'%s\t','tf_channel_flag');
for s_ii=1:length(v_chan_flag)
    fprintf(f_par,'%s',int2str(v_chan_flag(s_ii)));
    if (rem(s_ii,10)==0)
        fprintf(f_par,'%s\n\t\t\t\t\t\t','');
    else
        fprintf(f_par,'%s\t','');
    end;
end;
fprintf(f_par,'%s\n\n','');

fprintf(f_par,'%s\t','tf_channel_ref');
for s_ii=1:length(v_chan_flag)
    fprintf(f_par,'%s',int2str(v_chan_ref(s_ii)));
    if (rem(s_ii,10)==0)
        fprintf(f_par,'%s\n\t\t\t\t\t\t','');
    else
        fprintf(f_par,'%s\t','');
    end;
end;
fprintf(f_par,'%s\n\n','');

% to copy
fprintf(f_par,'%s\t%s\n','eegstat_type_stat',int2str(1*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','eegstat_time_hw',int2str(1*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','eegstat_time_step',int2str(1*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','eegstat_baseline_flag',int2str(0*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','eegstat_baseline_start',int2str(-300*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','eegstat_baseline_stop',num2str(-100*ones(1,s_nbevtype))); %  JP ?????

fprintf(f_par,'%s\t%s\n','tfstat_type_stat',int2str(1*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','tfstat_baseline_flag',int2str(0*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','tfstat_freq_hw',int2str(3*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','tfstat_freq_step',int2str(3*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','tfstat_time_hw',int2str(50*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','tfstat_time_step',int2str(50*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','tfstat_baseline_start',int2str(-300*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','tfstat_baseline_stop',num2str(-100*ones(1,s_nbevtype))); %  JP ?????

fprintf(f_par,'%s\t%s\n','tf_freq_start',int2str(1*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','tf_freq_stop',int2str(200*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','tf_freq_step',int2str(3*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','tf_nb_sample_blackman',int2str(30*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','tf_wavelet_type',int2str(1*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','tf_morlet_m',int2str(7*ones(1,s_nbevtype))); %  JP ?????


fclose(f_par);




