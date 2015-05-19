function LYON_localize_PREF()

a_todo = 'tall';
% a_elecdat = '/home/bastinj/MATLAB/tb_grenoble/localizer/elec.dat';
% a_base =['/home/bastinj/Documents/BVS/'];
a_elecdat = 'E:\ALIZEE\EPILEPSY\ALIZEE\elec.dat';                                
a_base =['E:\ALIZEE\EPILEPSY\ALIZEE\'];     
cd(a_base);

a_path = uigetdir(a_base, 'Select Folder to Analyze');
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

if (strfind(a_task,'PREF'))
    % change pos
%     a_oldposfile = [a_path FILESEP a_exp '_' a_task '_raw.pos'];
%     a_newposfile = [a_path FILESEP a_exp '_' a_task '.pos'];
%     a_oldposfile = [a_path FILESEP a_task '_raw.pos'];
%     a_newposfile = [a_path FILESEP a_task '.pos'];
    oldpos=dir([a_path FILESEP '*_raw.pos']);
    a_oldposfile = [a_path FILESEP oldpos.name];
    a_newposfile = [a_oldposfile([1:end-8 end-3:end])];
   
    loc_chg_pos(a_oldposfile,a_newposfile);
    a_logfile=[a_path FILESEP 'Loca_preferences.log'];
    chg_pos_locapref(a_newposfile,a_logfile)
% OK, I created a new posfile 'named': PAT_CODE_locapref.pos
  a_newposfile = strrep(a_newposfile,'.pos','_locapref.pos');
    % create .conf and .par
%     a_name = [a_path FILESEP a_exp '_' a_task];
%     a_name = [a_path FILESEP a_task];
a_name = a_oldposfile([1:end-8]);
%     loc2_e2cp_pref(a_name,a_elecdat); %VUTILE POUR PRE-PARAMETER DES FICHIERS QUI SERVENT A FAIRE DU
%     TEMPS-FREQUENCE AVEC L'OUTIL "ELAN" développé à LYON 

%PERMET AUSSI DE DETECTER QUAND LES CLINICIENS METTENT DES "zeros" partout
%dans le nom des électrodes : à corriger pour les scripts de JP par
%après...

%STEP 1 : FILTER the signal and then extract gamma band activity :)
%UNCOMMENT HERE (I already computed GAMMA)

eeg2env2(a_name,a_newposfile,(50:10:150),0);%gamma = done !   
eeg2env2_theta(a_name,(4:8),a_newposfile);%for theta analyses, we modified the way LFP signal is filtered before hilbert transfrom with
% Karim Jerbi.
eeg2env2_theta(a_name,(8:12),a_newposfile);%for theta analyses, we
% modified the way LFP signal is filtered before hilbert transfrom with
eeg2env2(a_name,(15:5:30),a_newposfile);%7BETA = done !
    
    v_code = [52 53]; %these event code need to be created in the posfile !
    a_code = strvcat('Pref_score (low)','Pref_score (high)');
    v_window_ms = [-500 3000];
%     loc_eeg2erp([a_name '.eeg'],a_newposfile,v_code,a_code,v_window_ms,20); %evoked response plots
    
   
% display correlations based on gamma enveloppes & pref_scores for different temporal
% smoothing (sorry alizee, it will indeed create many files in the same
% folder...BUT...everything is in the file names !
    

%NEXT: draw correlation values we just need sf_s to adjust the downsampling
%value (in JP LAchaux lab, we downsample gamma time serie to 64 Hz to
%increase the speed of all computations ; of course we can re-extrcat gamma
%at higher resolution if needed (theoretically grounded); not the case for
%LOCA_PREF...
%before everything, We will need the sampling rate of the EEGfile
% EXTRACT SAMPLING RATE INFO FROM ENT FILE
a_entfile =[a_name '.eeg.ent'];
f_old=fopen([a_entfile ],'r');
for s_i=1:9
    a_line=fgetl(f_old);
end;
s_fs = 1/str2num(a_line);
s_fs = round(s_fs); % sampling frequency
fclose(f_old);

switch s_fs
    case 256
        d='ds4';
    case 512
        d='ds8';
    case 1024
        d='ds16';
    case 2048
        d='ds32';
end
%DONE :%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


a_eeg = [a_name,'_f50f150_' d '_sm0.eeg'];%try not to smooth to see the results !!
% a_pos = strrep(a_newposfile,'.pos',['_locapref_' d '.pos']);
a_pos =[ a_name '_locapref_ds16.pos'];  
v_code = [51];% we need to compute correlations between gamma and preference score (which is in column 3)
      
      for j=1:4
          switch j
              case 1
      a_eeg = [a_name,'_f50f150_ds16_sm0.eeg'];%try not to smooth to see the results !!
              case 2
                        a_eeg = [a_name,'_f50f150_ds16_sm250.eeg'];%try not to smooth to see the results !!
              case 3
                        a_eeg = [a_name,'_f50f150_ds16_sm1000.eeg'];%try not to smooth to see the results !!
              case 4
                        a_eeg = [a_name,'_f50f150_ds16_sm500.eeg'];%try not to smooth to see the results !!
          end
      v_code = [51];
a_code = strvcat('PREF FOOD RATING');
v_window_ms = [-500 2500];vt_stats=v_window_ms;
loc_env2plot_rvalues_locaPREF_GLM(a_eeg,a_pos,v_code,a_code,v_window_ms,20); 
      end
end;