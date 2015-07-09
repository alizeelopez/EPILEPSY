function convert_cns2eeg(folder_ncs)

% this script creates an eeg file from .cns file.
% clear all
% 
% cd /mnt/data2/5-EPILEPSY/Analysis/Neurons
% 
% subject = 4;
% pathstr = ['/mnt/data2/5-EPILEPSY/Data/Subj_' num2str(subject)];
% curdir = dir(pathstr);
% 
% %%% retrieve the subject's name
% subjectname = curdir(3).name;
% 
% %%% load subject's specific informations
% eval(subjectname);
% 
% toRun = fullfile(pwd, ['patient_' subjectname]);
% run(toRun)
% 
% cd(subjectdata.rawdata_dir)
cd(folder_ncs)
index_folder=find(ismember(folder_ncs,'\'),2,'last');
subject_folder=folder_ncs(1:index_folder(1));

index_name_subj=find(ismember(subject_folder,'\'),2,'last');
subjectname=subject_folder(index_name_subj(1)+1:end-1);
mkdir([subject_folder '/eeg_rawdata_AL/'])
filename_eeg = [subject_folder '/eeg_rawdata_AL/' subjectname '.eeg'];

elec_list = dir([folder_ncs '\' '*.ncs']);
elec_cell = struct2cell(elec_list);
names = elec_cell(1,:)';
[~,index] = sort(upper(names));
names = names(index);




f_id = fopen('elec.dat');
A = textscan(f_id,'%d%d%d%s');
fclose(f_id);
elec(:,1)=A{1}([540:558, 1001:1255]);
elec(:,2)=A{2}([540:558,  1001:1255]);
elec(:,3)=A{3}([540:558,  1001:1255]);

elen=A{4}([540:558,  1001:1255]);

    % Get all .ncs files in the current folder
%     nbr = 1;
%     for id = 1:length(elec_list)
%         % Get the file name (minus the extension)
%         old_name = elec_list(id).name;
%         [~, f] = fileparts(old_name);
%         underscore=strfind(f,'_');
%         if length(underscore)==3
%             elec_name{id} = f(underscore(3)+1:length(f));
%         elseif length(underscore)==4
%             elec_name{id} = f(underscore(3)+1:underscore(4)-1);
%         end
%         
%         if id==1
%             Names{1,1}= elec_name{id};
%         else
%             if strcmp(elec_name{id},elec_name{id-1})==0
%                 nbr = nbr+1;
%                 Names{1,nbr}=elec_name{id};
%             end
%         end
%         
%         elec_nbr(id) = nbr;
%         % Convert to number
%         num = str2double(elec_name{id});
%         if isnan(num)
%             % If is not numeric, rename
%             if length(underscore)==4
%                 new_name = [old_name(1:underscore(3)), num2str(nbr), old_name(underscore(4):length(old_name))];
%             else
%                 new_name = [old_name(1:23), num2str(nbr)];
%             end
%             movefile(old_name, new_name);
%             keep_trace{id}.old_name=old_name;
%             keep_trace{id}.new_name=new_name;
%             
%         end
%     end
%     save('keep_trace.mat',keep_trace)
%     for n=1:length(Names)
%     Names{3,n}=num2str(n);
%     Names{2,n}=input(['Name for ' Names{1,n} '?']);
%     end
%     save([subject_folder(1:end-length(subjectname)-1) 'Electrodes\electrode_numbers.mat'],'Names')
%%

load([subject_folder(1:end-length(subjectname)-1) 'Electrodes\electrode_numbers.mat'])
Names2=Names(2,:);

for ch = 1:length(names)
    file = cell2mat(names(ch));
    down = strfind(file,'_');
    if length(down)==4
        elec_nbr = str2double(file(down(3)+1:down(4)-1));
    else
        elec_nbr = str2double(file(down(3)+1:strfind(file,'.')-1));
    end
    electrode{ch} = Names2{elec_nbr};
    
    if ch>1 
        if ~strcmp(electrode{ch},electrode{ch-1}) % if this is a new electrode
            n_ch=1;
        else
            n_ch=n_ch+1;
        end
    else
        n_ch=1;
    end
    
    if strcmp(electrode{ch}, 'ECG1')
        channel ='ECG1.-1' ;
    elseif strcmp(electrode{ch}, 'Oz')
        channel = 'Oz.2237';
    else
        clear b;
        b=find(strncmp(electrode{ch},elen,1));
        channel = [elen{b(n_ch)} '.' num2str(elec(b(n_ch),1))];
    end
    
    % extract data
    cfg=[];
    cfg.datafile = names{ch};
%     cfg.trl = [subjectdata.begsample subjectdata.endsample 0];
    cfg.lpfilter      = 'no';
    cfg.hpfilter      = 'no';
    cfg.bpfilter      = 'no';
    cfg.bsfilter      = 'no';
    cfg.dftfilter     = 'no';
    cfg.medianfilter  = 'no';
    data = ft_preprocessing(cfg); % does not do anything appart from getting the data
    
    % resample the data to 512 Hz
    cfg.resamplefs = 512;
    cfg.detrend = 'no';
    cfg.demean  = 'no';    
    dat = ft_resampledata(cfg, data);
    
    % m_data should be 1 channel per line
    m_data(ch,:)  = cell2mat(data.trial);
    m_events      = [];
    s_fs          = dat.fsample;
    str_ori_file1 = 'Neuralynx';
    str_ori_file2 = 'Paris';
    v_label{ch}   = channel;
    v_channel_type{ch} = 'electrode EEG';
    v_channel_unit{ch} = 'microV';
end

% filename = ['/mnt/data2/5-EPILEPSY/Data/Subj_3/02135/eeg_rawdata/' subjectname '.eeg'];
 mat2eeg(m_data, filename_eeg, m_events, str_ori_file1, str_ori_file2, s_fs, v_label, v_channel_type,v_channel_unit)

% filename = ['/mnt/data2/5-EPILEPSY/Data/Subj_3/02135/R_test/test/' subjectname '.eeg'];
% mat2eeg(m_data, filename, m_events, str_ori_file1, str_ori_file2, s_fs, v_label, v_channel_type,v_channel_unit)

% load('/mnt/data2/5-EPILEPSY/Data/Subj_3/02135/R_test/02135.pos');
% m_pos=s_fs*X02135(:,1)/512;


% % find session start lines in sample vector
% evt_cell  = struct2cell(event);
% code      = cell2mat(evt_cell(2,:)');
% sample    = cell2mat(evt_cell(3,:)'); % samples in 1/4000 s
% sample2   = cell2mat(evt_cell(3,:)')/4; % samples in ms for comparison with .mat
% sample(code~=255)=[];
% sample2(code~=255)=[];
% code(code~=255)=[];
% sess_lines = find(sample2(2:length(sample2))-sample2(1:length(sample2)-1)<500);
% cd('/mnt/data2/5-EPILEPSY/Data')
% load('Total_data_sub3')
% for t=1:30
%     m_pos(t)=512*(sample(sess_lines(3)+t)+4*(Pleas_onset(t,1,1)-Pleas_cross_onset(t,1,1)))/4000;
% %     m_pos(t)=512*(sample(sess_lines(1)+t))/4000;
% end
% 
% figure
% hold
% for n=1:1
%     subplot(1,1,n)
%     hold on
%     for trig=1:30
%         plot([m_pos(trig) m_pos(trig)],[-300, 200],'r')
%     end
%     
%     % channel 42 for sub 2
%     plot(m_data(n,1:m_pos(trig)))
%     axis([m_pos(1) m_pos(trig) -400 200])
% end