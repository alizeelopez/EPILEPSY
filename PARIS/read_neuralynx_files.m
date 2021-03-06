function event=read_neuralynx_files(path_neuralynx)

% function event=read_neuralynx_files(path_neuralynx,electrode_names)

% ensure that we don't mix up subjects

% prior to reading the data, make sure you rename the files so elec numbers
% are 01,02...,10...

curdir = pwd;
do_conversion = 1;

% Directory of the files
% d = '/mnt/data2/5-EPILEPSY/Data/Subj_4/02141/Neuralynx_files/02141_2014-03-25_10-15_old';
% d = '/mnt/data2/5-EPILEPSY/Data/Subj_4/02141/Neuralynx_files/02141_2014-03-25_10-30';
d=path_neuralynx;

% from name to number
if do_conversion==1
    cd(d)
    % Get all .ncs files in the current folder
    files = dir('*.ncs');
    nbr = 1;
    for id = 1:length(files)
        % Get the file name (minus the extension)
        old_name = files(id).name;
        [~, f] = fileparts(old_name);
        underscore=strfind(f,'_');
        if length(underscore)==3
            elec_name{id} = f(underscore(3)+1:length(f));
        elseif length(underscore)==4
            elec_name{id} = f(underscore(3)+1:underscore(4)-1);
        end
        
        if id==1
            Names{1}= elec_name{id};
        else
            if strcmp(elec_name{id},elec_name{id-1})==0
                nbr = nbr+1;
                Names{nbr}=elec_name{id};
            end
        end
        
        elec_nbr(id) = nbr;
        % Convert to number
        num = str2double(elec_name{id});
        if isnan(num)
            % If is not numeric, rename
            if length(underscore)==4
                new_name = [old_name(1:underscore(3)), num2str(nbr), old_name(underscore(4):length(old_name))];
            else
                new_name = [old_name(1:23), num2str(nbr)];
            end
            movefile(old_name, new_name);
        end
    end
    
    %            1      2      3      4     5     6     7       8     9      10     11     12    13     14     15   16    17
    % Names = {'AmT2','ECG1','FT10','FT9','Fp1','Fp2','HaT1','HmT2','InT1','Insa','Insm','Insp','OrF','OrFa','Oz','TBm','Tpol'};
    
%     Names2=electrode_names;
     Names2  = {'A',   'ECG1','FT10','FT9','Fp1','Fp2', 'B',   'C',   'D',   'E',    'G',  'H',   'J',  'O',  'Oz', 'T',   'V'};
    save('File_conversion','elec_name','elec_nbr','Names', 'Names2');
    
    % convert event to matfile
    event_file=dir([path_neuralynx filesep '*Events.nev']);
%     event=ft_read_event('02141_2014-03-25_09-34_Events.nev');
    event=ft_read_event(event_file(1).name);

    files = dir('*.ncs');
    
    for e=1:length(event)
        if isnan(event(1,e).sample)
            stamps  = event(1,e).timestamp;
            hdr = ft_read_header(files(1).name);
            sample4 = (double(stamps)-double(hdr.FirstTimeStamp))/double(hdr.TimeStampPerSample)+1;
            event(1,e).sample=uint64(sample4);
            E(e)=event(1,e).sample;
        end
    end
    
%     save('Events.mat','event');
%     cd ..
%     cd ..
%     save('Events.mat','event');
%     save('File_conversion','elec_name','elec_nbr','Names', 'Names2');
end

%%

% clear subjectdata
% 
% % define the filenames, parameters and other information that is subject specific
% subjectdata.rawdata_dir = d;
% clear d
% padding = 0;
% 
% % define first and last samples in which to look
% % cd('/mnt/data2/5-EPILEPSY/Data/Subj_4/02141')
% cd(path_neuralynx)
% % load('File_conversion');
% % load('Events.mat');
% evt_cell  = struct2cell(event);
% code      = double(cell2mat(evt_cell(2,:)'));
% sample    = double(cell2mat(evt_cell(3,:)')); % samples in 1/4000 s
% sample2   = double(cell2mat(evt_cell(3,:)'))/4; % samples in ms for comparison with .mat
% sample(code~=255)=[];
% sample2(code~=255)=[];
% code(code~=255)=[];
% 
% % find session start lines in sample vector
% sess_lines = find(sample2(2:length(sample2))-sample2(1:length(sample2)-1)<500);
% % 
% cd('/mnt/data2/5-EPILEPSY/Data')
% load('Total_data_sub4')

% align samples to take in cns file on the sample of the first m_pos event
%  = onset of the 1st item to date
% subjectdata.begsample = sample(sess_lines(1)+1)+4*(Age_onset(1,1,1)-Age_cross_onset(1,1,1))-padding*4000; 
% start X = padding seconds before the start of the first trial
% subjectdata.endsample = sample(length(sample))+padding*4000; % look until X = padding seconds after the start of the last trial



%%% identify electrodes
% load('/mnt/data2/5-EPILEPSY/Analysis/Neurons/ActiCap_Elec')
% load('E:\ALIZEE\EPILEPSY\ALIZEE\EPILEPSY\PARIS\ActiCap_Elec.mat')
% 
% for n=1:length(Names)
%     clear a
%     a = find(ismember(elec_name,Names{n}));
%     for m=1:length(a)
%         elec_name1{a(m)}=Names2{n};
%     end
% end
% load elec.dat names/conversion data
% f_id = fopen('elec.dat');
% A = textscan(f_id,'%d%d%d%s');
% fclose(f_id);
% elec(:,1)=A{1}([540:558, 1001:1255]);
% elec(:,2)=A{2}([540:558,  1001:1255]);
% elec(:,3)=A{3}([540:558,  1001:1255]);
% 
% cd(curdir)
% elen=A{4}([540:558,  1001:1255]);