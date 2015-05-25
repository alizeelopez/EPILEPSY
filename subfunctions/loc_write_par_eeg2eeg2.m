function loc_write_par_eeg2eeg2(a_npar,ss_bipole,s_fs,s_downsamp,a_listevent)

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

s_nbchannel = length(ss_bipole);

s_nombre_evt = length(a_listevent);
%a_listevent = (1:s_nombre_evt);

a_parfile_bipo=a_npar;
s_fs = round(s_fs/s_downsamp);

[a_path,a_thispre,a_ext] = fileparts(a_npar);
s_nbevtype = length(a_listevent);
v_evt=a_listevent;
v_evt=v_evt(:);

% start with the bipo
f_par=fopen(a_parfile_bipo,'w');

% TO BE ADJUSTED TO THE EXPERIMENT
fprintf(f_par,'%s\t%s\n','fileprefix',[a_thispre]); %  JP ?????
fprintf(f_par,'%s\t%s\n','nb_eventcode',int2str(s_nbevtype)); %  JP ?????
fprintf(f_par,'%s\t%s\n','list_eventcode',int2str(v_evt')); %  JP ?????
fprintf(f_par,'%s\t%s\n','prestim_nbsample',int2str(s_fs*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','poststim_nbsample',int2str(3*s_fs*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','baseline_eventcode',int2str(0*ones(1,s_nbevtype))); %  JP ?????

fprintf(f_par,'%s\t','baseline_msec_start');
for s_i = 1:s_nbevtype
    fprintf(f_par,'%s\t','-400');
end;
fprintf(f_par,'%s\n','');

fprintf(f_par,'%s\t','baseline_msec_start');
for s_i = 1:s_nbevtype
    fprintf(f_par,'%s\t','-100');
end;
fprintf(f_par,'%s\n','');



% v_bigpro contains the relative indices, after reordering, of all the big brothers
v_chan_flag=zeros(1,s_nbchannel+2);
v_chan_ref=zeros(1,s_nbchannel+2);
v_chan_flag(1:s_nbchannel)=1;

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

% to copy
fprintf(f_par,'%s\t%s\n','eegstat_type_stat',int2str(1*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','eegstat_time_hw',int2str(1*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','eegstat_time_step',int2str(1*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','eegstat_baseline_flag',int2str(0*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','eegstat_baseline_start',int2str(-400*ones(1,s_nbevtype))); %  JP ?????
fprintf(f_par,'%s\t%s\n','eegstat_baseline_stop',num2str(-100*ones(1,s_nbevtype))); %  JP ?????


fclose(f_par);
