function stat_electrodes=explicit_representation_main(a_path,task,stim,facteur,frequence,smoothing,do_erp,do_gamma,load_stat_electrode)

% try

a_todo = 'tall';
% a_elecdat = '/home/bastinj/MATLAB/tb_grenoble/localizer/elec.dat';
% a_base =['/home/bastinj/Documents/BVS/'];
a_elecdat = 'E:\ALIZEE\EPILEPSY\ALIZEE\elec.dat';
a_base =['E:\ALIZEE\EPILEPSY\ALIZEE\'];
cd(a_base);

% a_path = uigetdir(a_base, 'Select Folder to Analyze');
if task==0
    dir_to_safe=[a_path '\Preproc_explicit'];
end
mkdir(dir_to_safe)
v_f = findstr(a_path,filesep);
s_3 = v_f(length(v_f));
s_2 = v_f(length(v_f)-1);
a_exp_task = a_path(s_3+1:length(a_path));
a_exp = a_path(s_2+1:s_3-1);
a_task = strrep(a_exp_task,[a_exp '_'],'');
FILESEP=filesep;
%%%%%%%%%%%%%%%%%
% XX.  PREF
%%%%%%%%%%%%%%%%
if task==0
    name_task='all';
end

diary([a_exp_task '_' name_task '_smoothing' num2str(smoothing) '_erp' num2str(do_erp) '_gamma' num2str(do_gamma)]);
disp('-----------------------------------------------------------------------------------------')
disp([' #################                    ' a_exp_task '              ##################'])


disp('-----------------------------------------------------------------------------------------')
disp(' #################            Processing pos and log file...           ##################')
oldpos=dir([a_path FILESEP '*BVS_raw.pos']);
a_oldposfile = [a_path FILESEP oldpos.name];
a_newposfile = [dir_to_safe FILESEP oldpos.name([1:end-8 end-3:end])];

loc_chg_pos(a_oldposfile,a_newposfile);
    %% Fusionne log files pour chaque t�che
    a_logfile_part{1}=[a_path FILESEP 'Copy_of_BVS_1b.log'];
    a_logfile_part{2}=[a_path FILESEP 'Copy_of_BVS_2b.log'];
    a_logfile_part{3}=[a_path FILESEP 'Copy_of_BVS_3b.log'];
    a_logfile_part{4}=[a_path FILESEP 'Copy_of_BVS_4b.log'];
    a_logfile_part{5}=[a_path FILESEP 'Copy_of_BVS_5b.log'];
    
    a_logfile=[a_path FILESEP 'Loca_preferences_all.log'];
    Merge_logfiles( a_logfile_part, a_logfile)
    
    chg_pos_main(a_newposfile,a_logfile);
    % OK, I created a new posfile 'named': PAT_CODE_locapref.pos
    a_new_pos_file_for_eeg=a_newposfile;
    a_newposfile = strrep(a_newposfile,'.pos','_locapref_all.pos');
    % create .conf and .par
    %     a_name = [a_path FILESEP a_exp '_' a_task];
    %     a_name = [a_path FILESEP a_task];
    [index_preproc, ind_end]=regexp((a_newposfile),'\Preproc_explicit\');
    
    if exist(strrep(a_new_pos_file_for_eeg([1:index_preproc-1 (ind_end+1):end]),'.pos','_D1.eeg.ent'),'file')
        
        a_name = [a_oldposfile([1:end-8]) '_D1'];
        a_new_name=[a_newposfile([1:end-17]) '_D1'];
    else
        a_name = [a_oldposfile([1:end-8])];
        a_new_name=[a_newposfile([1:end-17])];
    end
    v_code = [100 101]; %these event code need to be created in the posfile !
    

%     loc2_e2cp_pref(a_name,a_elecdat); %VUTILE POUR PRE-PARAMETER DES FICHIERS QUI SERVENT A FAIRE DU
%     TEMPS-FREQUENCE AVEC L'OUTIL "ELAN" développé à LYON
disp('#################         Processing of pos and log file done !         #################')
disp('-----------------------------------------------------------------------------------------')
disp('-----------------------------------------------------------------------------------------')

%PERMET AUSSI DE DETECTER QUAND LES CLINICIENS METTENT DES "zeros" partout
%dans le nom des électrodes : à corriger pour les scripts de JP par
%après...

%STEP 1 : FILTER the signal and then extract gamma band activity :)
%UNCOMMENT HERE (I already computed GAMMA)
% EXTRACT SAMPLING RATE INFO FROM ENT FILE
a_entfile =[a_name '.eeg.ent'];
f_old=fopen([a_entfile ],'r');
for s_i=1:9
    a_line=fgetl(f_old);
end;
s_fs = 1/str2num(a_line);
s_fs = round(s_fs); % sampling frequency
fclose(f_old);

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



smooth_val=smoothing; %0
v_freq=frequence; % (50:10:150)
if do_gamma==1%~exist(new_file,'file')
    
    disp('###############       Start filter and extract gamma band activity       ################')
    
    eeg2env2(a_name,a_newposfile,v_freq,smooth_val);%gamma = done !
    
else
    disp('###############     Filter and extract gamma band activity not asked      ################')
    
end
%     eeg2env2_theta(a_name,(4:8),a_newposfile);%for theta analyses, we modified the way LFP signal is filtered before hilbert transfrom with
%     % Karim Jerbi.
%     eeg2env2_theta(a_name,(8:12),a_newposfile);%for theta analyses, we
% modified the way LFP signal is filtered before hilbert transfrom with
%     eeg2env2(a_name,(15:5:30),a_newposfile);%7BETA = done !

a_code = strvcat('Pref_score (low)','Pref_score (high)');
v_window_ms = [-200 1500];
if exist([a_name 'stats_electrodes.mat'],'file') && load_stat_electrode==1
    load([a_name 'stats_electrodes.mat'])

else
    stat_electrodes=[];
end

if do_erp==1
    disp('-----------------------------------------------------------------------------------------')
    disp('###############             Start evoked response activity               ################')
    stat_electrodes=loc_eeg2erp_main([a_name '.eeg'],a_newposfile,v_code,a_code,v_window_ms,20,stat_electrodes); %evoked response plots
end


save([a_name 'stats_electrodes.mat'],'stat_electrodes')

% display correlations based on gamma enveloppes & pref_scores for different temporal
% smoothing (sorry alizee, it will indeed create many files in the same
% folder...BUT...everything is in the file names !


%NEXT: draw correlation values we just need sf_s to adjust the downsampling
%value (in JP LAchaux lab, we downsample gamma time serie to 64 Hz to
%increase the speed of all computations ; of course we can re-extrcat gamma
%at higher resolution if needed (theoretically grounded); not the case for
%LOCA_PREF...
%before everything, We will need the sampling rate of the EEGfile


%DONE :%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('-----------------------------------------------------------------------------------------')
disp(' #################              Start regressing eeg signal             #################')



a_eeg = [a_new_name,'_f' int2str(min(v_freq)) 'f' int2str(max(v_freq)) '_ds' int2str(s_downsamp) '_sm' num2str(smooth_val) '.eeg'];%try not to smooth to see the results !!
    freq_range_and_smooth=['f' int2str(min(v_freq)) 'f' int2str(max(v_freq)) '_ds' int2str(s_downsamp) '_sm' num2str(smooth_val)];

%     v_code = [51];
% v_window_ms = [-500 2500];
vt_stats=v_window_ms;
if task==0
    
    if strcmp(stim,'all')
        v_code = [31 41 51];
    end
    if strcmp(facteur,'Pref')
        ind_betas=1;
    elseif strcmp(facteur,'Conf')
        ind_betas=2;
    elseif strcmp(facteur,'RT')
        ind_betas=3;
    end
    
    a_code = strvcat(['All_task_on_' stim '_' facteur]);
    if exist(strrep(a_new_pos_file_for_eeg([1:index_preproc-1 (ind_end+1):end]),'.pos','_D1.eeg.ent'),'file')
        
        a_pos =[ a_new_name(1:end-3) '_locapref_all_' 'ds' int2str(s_downsamp) '.pos'];
    else
        a_pos =[ a_new_name(1:end) '_locapref_all_' 'ds' int2str(s_downsamp) '.pos'];
    end
    options.GLM_facteurs={facteur};
    options.GLM_facteurs_index=[ind_betas];
    options.folder_to_save_figures=['\GLM_' facteur '_on_all_' stim];
    mkdir([dir_to_safe options.folder_to_save_figures])
    stat_electrodes=loc_env2plot_rvalues_loca_main_GLM(a_eeg,a_pos,v_code,a_code,v_window_ms,20,options,stat_electrodes,freq_range_and_smooth);
    
    save([a_name 'stats_electrodes.mat'],'stat_electrodes')

end



disp('-----------------------------------------------------------------------------------------')
disp(' #################                       Finish !                       #################')
disp(['Figures have been saved for ' a_exp_task '_' name_task ' smoothing ' num2str(smoothing)])

diary off


% catch ME
%     disp('-----------------------------------------------------------------------------------------')
%     disp(' #################                       Error !                       #################')
%     diary off
%
%     rethrow(ME)
%
% end
