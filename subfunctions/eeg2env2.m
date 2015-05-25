function eeg2env2(a_rootname, a_posfiles, v_freq, v_smoothing_ms)

% this function will take each channel of the eeg file, and transform it
% into a new eeg file where each signal is replaced with the enveloppe in
% the s_fmin:s_fmax frequency range.
% v_smoothing_ms: smoothing windows in ms
% by default, the signal will be filtered in bands of 5 Hz, then normalized
% relative to the mean and std over the entire interval. then averaged
% across frequencies.
%
% if s_smoothing_ms is non-zero, it will, as a last step, convolve the
% signal with a window of that length.
%
% note : all signals will have the same standard deviation at the end


%% 1. EXTRACT INFO FROM .ENT FILE

% read basic parameters, header file, etc.
a_eegfile = [a_rootname '.eeg'];
a_eegfile_orig = a_eegfile;
f_old=fopen([a_rootname '.eeg.ent'],'r');

a_new_rootname=[a_rootname(1:find(ismember(a_rootname,'\'),1,'last')) 'Preproc' a_rootname(find(ismember(a_rootname,'\'),1,'last'):end)];

for s_i=1:9
    a_line=fgetl(f_old);
end

s_fs = 1/str2num(a_line);
s_fs = round(s_fs); % sampling frequency
v_smoothing_sam = round(s_fs*v_smoothing_ms/1000);

% we will downsample the data to 64 samples per second
if (s_fs==512)
    s_downsamp = 8;
elseif (s_fs==1024)
    s_downsamp = 16;
elseif (s_fs==2048)
    s_downsamp = 32;
else
    s_downsamp = 4; % 256 Hz
end;


% here we can change the pos files, if not empty
for s_p = 1:size(a_posfiles,1)
    a_opos = ((a_posfiles(s_p,:)));
    a_npos = strrep(a_opos,'.pos',['_ds' int2str(s_downsamp) '.pos']);
    eeg2eeg_chg_pos(a_opos,a_npos,s_downsamp); % CHECK : make sure you get induced response
end


% extract channel names
a_line=fgetl(f_old);

clear ss_channel
% how many channels ?
s_nbchannel=str2num(a_line);
s_nbchannel_orig = s_nbchannel;
for s_c=1:s_nbchannel-2
    % read the name of the actual data channels
    s_indice=-1;
    a_line=fgetl(f_old);
    % the name is something like Cz.10/ v'1.1969;
    s_point=find(a_line=='.');
    if (~isempty(s_point))
        a_elecname=a_line(1:s_point(1)-1);
        a_afterpoint = a_line(s_point(1):length(a_line));
    else
        a_afterpoint = '.-1';
    end;
    a_elecname=upper(a_elecname);
    a_elecname=strrep(a_elecname,' ','');
    a_afterpoint=upper(a_afterpoint);
    a_afterpoint=strrep(a_afterpoint,' ','');
    
    % now we have to find it in the list of elec.dat
    ss_channel.values(s_c).value=a_elecname;
    ss_channel.values(s_c).afterpoint=a_afterpoint;
end

fclose(f_old);

% decide how to create bipoles from the electrode names
[v_neworder, v_neworder_rev, m_bipole]=eeg_loc_montage(ss_channel);

m_bipole_origid = v_neworder_rev(m_bipole); % for each bipole, specifies the indice of the big and litle brothers according to indexing in ss_channel and eeg file



%% 2. CREATE A STRUCTURE THAT CONTAINS DATA AND BIPOLE NAME FOR EACH BIPOLE

% now that we have that the idea is to process each bipole separately
clear ss_bipole;
s_nbip = size(m_bipole,1);

v_changroup = (1:5:s_nbip);
v_changroup = [v_changroup s_nbip+1];

for s_g = 1:(length(v_changroup)-1)
    %for s_g = 1:1
    
    disp(['processing bipoles ' int2str(v_changroup(s_g)) ' to ' int2str(v_changroup(s_g+1)-1)]);
    clear m_data_sc;
    disp('retrieving data');
    v_c = (v_changroup(s_g):(v_changroup(s_g+1)-1));
    s_nbc = length(v_c);
    
    
    
    % READ DATA AND INITIALIZE
    for s_i = 1:length(v_c)
        
        s_bip = v_c(s_i);
        
        disp(['reading data from bipole ' int2str(s_bip) ' out of ' int2str(s_nbip)]);
        
        ss_bipole(s_bip).name = [ss_channel.values(m_bipole_origid(s_bip,1)).value];
        ss_bipole(s_bip).afterpoint = [ss_channel.values(m_bipole_origid(s_bip,1)).afterpoint];
        m_data_tool = rd_eeg(a_eegfile_orig,m_bipole_origid(s_bip,:),s_nbchannel_orig);
        v_data_sc = m_data_tool(:,1)-m_data_tool(:,2);
        
        if (s_i==1)
            % initialize size
            m_data_sc = zeros(length(v_data_sc),length(v_c));
        end
        m_data_sc(:,s_i) =v_data_sc;
    end
    
    v_t = (1:size(m_data_sc,1))/s_fs;
    v_s = (1:s_downsamp:length(v_data_sc));
    s_nbsample = length(v_s);
    s_nbchannel = s_nbip+2;
    
    for sm=1:length(v_smoothing_ms)
        m_allsc_s{sm} = zeros(length(v_s),s_nbc);
    end
    
    
    % INITIALIZE THE EEG FILES
    if (s_g == 1)
        for sm=1:length(v_smoothing_ms)
            a_eegfile = [a_new_rootname '_f' int2str(min(v_freq)) 'f' int2str(max(v_freq)) '_ds' int2str(s_downsamp) '_sm' num2str(v_smoothing_ms(sm)) '.eeg'];
            eeg2eeg_write_eeg_init(a_eegfile,s_nbsample,s_nbchannel)
        end
        %         a_eegfile = [a_rootname '_f' int2str(min(v_freq)) 'f' int2str(max(v_freq)) '_ds' int2str(s_downsamp) '_sm250.eeg'];
        %         eeg2eeg_write_eeg_init(a_eegfile,s_nbsample,s_nbchannel)
    end
    
    
    % first shot : downsample
    v_ts = v_t(v_s);
    
    m_datae = zeros(length(v_s),size(m_data_sc,2));
    
    disp('filtering data');
    for s_f = 1:length(v_freq)-1
        s_fmi = v_freq(s_f);
        s_fma = v_freq(s_f+1);
        m_dataf=bpfilter(m_data_sc,s_fmi,s_fma,s_fs);
        m_dataf=bf_envhilb(m_dataf);
        m_dataf=m_dataf(v_s,:);
        m_datafm = mean(m_dataf(round(length(v_s)/4):3*round(length(v_s)/4),:)); % remove edge effects
        v_fff = find(m_datafm==0);
        if (~isempty(v_fff))
            m_datafm(v_fff)=1;
        end
        m_datafm = repmat(m_datafm,size(m_dataf,1),1);
        m_dataf = 100*m_dataf./m_datafm; % divide by mean for every frequency - m_dataf is now expressed in % (100 % = average power in the band)
        m_datae = m_datae + m_dataf;
    end
    
    m_datae = 10*m_datae/(length(v_freq)-1); % all values in per thousands (pour mille)  of mean value % DB
    %     m_datae;
    %     v_k = (round(size(m_datae,1)/3):2*round(size(m_datae,1)/3));
    %     m_tool = m_datae(v_k,:);
    %     v_tool_s = std(m_tool);
    %     v_tool_m = mean(m_tool);
    %     m_tool_s = repmat(v_tool_s,size(m_datae,1),1);
    %     m_tool_m = repmat(v_tool_m,size(m_datae,1),1);
    %     m_datae = (m_datae-m_tool_m)./m_tool_s;
    %     m_datae = m_datae*s_sigma; % this is simply to ensure that the data that will be written are centered and have a standard deviation compatible with elan
    
    % then we smooth the data to create 6 matrixes.
    for sm=1:length(v_smoothing_ms)
        if v_smoothing_ms(sm)==0
            m_allsc_s{sm}(:,1:length(v_c)) = m_datae;
        else
            m_datac = conv2(m_datae',ones(1,v_smoothing_sam(sm)/s_downsamp)/(v_smoothing_sam(sm)/s_downsamp),'same');
            m_datac = m_datac'; % now sample x channel, just like m_datae
            m_allsc_s{sm}(:,1:length(v_c)) = m_datac;
        end
    end
    %     % 250
    %     s_sa = 2;
    %     m_datac = conv2(m_datae',ones(1,v_smoothing_sam(s_sa)/s_downsamp)/(v_smoothing_sam(s_sa)/s_downsamp),'same');
    %     m_datac = m_datac'; % now sample x channel, just like m_datae
    %     m_allsc_s250(:,1:length(v_c)) = m_datac;
    
    
    
    % THEN WE WRITE TO EEG FILE
    disp('writing data');
    
    for sm=1:length(v_smoothing_ms)
        v_fff = find(isnan(m_allsc_s{sm}));
        if (~isempty(v_fff))
            m_allsc_s{sm}(v_fff)=0;
        end
        a_eegfile = [a_new_rootname '_f' int2str(min(v_freq)) 'f' int2str(max(v_freq)) '_ds' int2str(s_downsamp) '_sm' num2str(v_smoothing_ms(sm)) '.eeg'];
        for s_o = 1:length(v_c)
            eeg2eeg_write_eeg_onechannel(m_allsc_s{sm}(:,s_o),a_eegfile,v_c(s_o),s_nbchannel);
        end
    end
    %     v_fff = find(isnan(m_allsc_s250));
    %     if (~isempty(v_fff))
    %         m_allsc_s250(v_fff)=0;
    %     end
    %     a_eegfile = [a_rootname '_f' int2str(min(v_freq)) 'f' int2str(max(v_freq)) '_ds' int2str(s_downsamp) '_sm250.eeg'];
    %     for s_o = 1:length(v_c)
    %         eeg2eeg_write_eeg_onechannel(m_allsc_s250(:,s_o),a_eegfile,v_c(s_o),s_nbchannel);
    %     end
end

% write technical channels
v_c = [length(ss_bipole)+1 length(ss_bipole)+2];
v_0 = zeros(size(m_allsc_s{1},1),1);

for sm=1:length(v_smoothing_ms)
    a_eegfile = [a_new_rootname '_f' int2str(min(v_freq)) 'f' int2str(max(v_freq)) '_ds' int2str(s_downsamp) '_sm' num2str(v_smoothing_ms(sm)) '.eeg'];
    for s_o = 1:length(v_c)
        eeg2eeg_write_eeg_onechannel(v_0,a_eegfile,v_c(s_o),s_nbchannel);
    end
end

% % s250
% a_eegfile = [a_rootname '_f' int2str(min(v_freq)) 'f' int2str(max(v_freq)) '_ds' int2str(s_downsamp) '_sm250.eeg'];
% for s_o = 1:length(v_c)
%     eeg2eeg_write_eeg_onechannel(v_0,a_eegfile,v_c(s_o),s_nbchannel);
% end


%% WRITE TO ENT FILE
disp('writing data');

for sm=1:length(v_smoothing_ms)
    a_eegfile = [a_new_rootname '_f' int2str(min(v_freq)) 'f' int2str(max(v_freq)) '_ds' int2str(s_downsamp) '_sm' num2str(v_smoothing_ms(sm)) '.eeg'];
    eeg2eeg_chg_eegent([a_rootname '.eeg.ent'],[a_eegfile '.ent'],ss_bipole,(s_fs/s_downsamp));
end
% a_eegfile = [a_rootname '_f' int2str(min(v_freq)) 'f' int2str(max(v_freq)) '_ds' int2str(s_downsamp) '_sm250.eeg'];
% eeg2eeg_chg_eegent([a_rootname '.eeg.ent'],[a_eegfile '.ent'],ss_bipole,(s_fs/s_downsamp));


% Create a par file that is convenient to calculate evoked
% potentials for all files
if (~isempty(a_posfiles))
    a_npar = strrep(a_opos,'.pos',['_ds' int2str(s_downsamp) '.par']);
    
    % Get the list of all event types in the pos file
    a_tpos = ((a_posfiles(1,:)));
    m_pos = load(a_tpos);
    v_tev = m_pos(:,2);
    v_listevent = [];
    for s_i = 1:length(v_tev)
        s_id = v_tev(s_i);
        v_f = find(v_listevent==s_id);
        if (isempty(v_f))
            v_listevent = [v_listevent s_id];
        end
    end
    v_listevent = sort(v_listevent);
    
    s_limit_event = min(3,length(v_listevent));
    loc_write_par_eeg2eeg2(a_npar,ss_bipole,s_fs,s_downsamp,v_listevent(1:s_limit_event));
    % I use only the first three events so that it is easy to change by hand
end