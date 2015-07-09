clear
close all
clc


smoothing=0;
freq=50:10:150;
load_old_stat_files=0;
do_erp=1;
compute_gamma=0;

%% LONG VERSIONS
 % DONE 07/07/2015:
% a_path='E:\ALIZEE\EPILEPSY\DATA_SAFE\GRE_PC_18NOV13G_BVS\PC_18NOV13G_BVS\PC_18NOV13G_BVS';
% stat_electrodes=explicit_representation_main(a_path,0,'all','Pref',freq,smoothing,do_erp,compute_gamma,load_old_stat_files);
% save('E:\ALIZEE\EPILEPSY\ALIZEE\EPILEPSY\Subjects_structures\GRE_PC_18NOV13G_BVS_structure_result','stat_electrodes')
% 
% a_path='E:\ALIZEE\EPILEPSY\DATA_SAFE\GRE_2014_DINg_BVS_D1\GRE_2014_DINg_BVS_D1\GRE_2014_DINg_BVS_D1';
% stat_electrodes=explicit_representation_main(a_path,0,'all','Pref',freq,smoothing,do_erp,compute_gamma,load_old_stat_files);
% save('E:\ALIZEE\EPILEPSY\ALIZEE\EPILEPSY\Subjects_structures\GRE_2014_DINg_BVS_D1_structure_result','stat_electrodes')


%% SHORT VERSIONS

 % DONE 07/07/2015:
% a_path='E:\ALIZEE\EPILEPSY\DATA_SAFE\GRE_2015_PETa\GRE_2015_PETa_PREF';
% stat_electrodes=localize_PREF(a_path,1,freq,smoothing,do_erp,compute_gamma,load_old_stat_files)
% save('E:\ALIZEE\EPILEPSY\ALIZEE\EPILEPSY\Subjects_structures\GRE_2015_PETa_PREF_structure_result','stat_electrodes')
% 
% a_path='E:\ALIZEE\EPILEPSY\DATA_SAFE\GRE_2014_LEBf_PREF\GRE_2014_LEBf_PREF\GRE_2014_LEBf_PREF';
% stat_electrodes=localize_PREF(a_path,1,freq,smoothing,do_erp,compute_gamma,load_old_stat_files)
% save('E:\ALIZEE\EPILEPSY\ALIZEE\EPILEPSY\Subjects_structures\GRE_2014_LEBf_PREF_structure_result','stat_electrodes')
% 
% a_path='E:\ALIZEE\EPILEPSY\DATA_SAFE\LYON\LOCA\LYONNEURO_2014_SIEJ_PREF';
% stat_electrodes=localize_PREF(a_path,1,freq,smoothing,do_erp,compute_gamma,load_old_stat_files)
% save('E:\ALIZEE\EPILEPSY\ALIZEE\EPILEPSY\Subjects_structures\LYONNEURO_2014_SIEJ_PREF_structure_result','stat_electrodes')
% 
% a_path='E:\ALIZEE\EPILEPSY\DATA_SAFE\LYON\LOCA\LYONNEURO_2014_LIBI_PREF';
% stat_electrodes=localize_PREF(a_path,1,freq,smoothing,do_erp,compute_gamma,load_old_stat_files)
% save('E:\ALIZEE\EPILEPSY\ALIZEE\EPILEPSY\Subjects_structures\LYONNEURO_2014_LIBI_PREF_structure_result','stat_electrodes')
% 
% 
% a_path='E:\ALIZEE\EPILEPSY\DATA_SAFE\GRE_2014_AUDe_PREF\GRE_2014_AUDe_PREF\GRE_2014_AUDe_PREF';
% stat_electrodes=localize_PREF(a_path,1,freq,smoothing,0,compute_gamma,load_old_stat_files)
% save('E:\ALIZEE\EPILEPSY\ALIZEE\EPILEPSY\Subjects_structures\GRE_2014_AUDe_PREF_structure_result','stat_electrodes')
% 

% a_path='E:\ALIZEE\EPILEPSY\DATA_SAFE\LYON\LOCA\LYONNEURO_2014_FEPk_PREF';
% stat_electrodes=localize_PREF(a_path,1,freq,smoothing,do_erp,1,load_old_stat_files)
% save('E:\ALIZEE\EPILEPSY\ALIZEE\EPILEPSY\Subjects_structures\LYONNEURO_2014_FEPk_PREF_structure_result','stat_electrodes')
% 
% a_path='E:\ALIZEE\EPILEPSY\DATA_SAFE\LYON\LOCA\LYONNEURO_2015_BOUa_PREF';
% stat_electrodes=localize_PREF(a_path,1,freq,smoothing,do_erp,1,load_old_stat_files)
% save('E:\ALIZEE\EPILEPSY\ALIZEE\EPILEPSY\Subjects_structures\LYONNEURO_2015_BOUa_PREF_structure_result','stat_electrodes')
% 

% a_path='E:\ALIZEE\EPILEPSY\DATA_SAFE\LYON\LOCA\LYONNEURO_2015_BOUc1_PREF';
% stat_electrodes=localize_PREF(a_path,1,freq,smoothing,do_erp,1,load_old_stat_files)
% save('E:\ALIZEE\EPILEPSY\ALIZEE\EPILEPSY\Subjects_structures\LYONNEURO_2015_BOUc1_PREF_structure_result','stat_electrodes')
% 

% a_path='E:\ALIZEE\EPILEPSY\DATA_SAFE\LYON\LOCA\LYONNEURO_2015_DUBb_PREF';
% stat_electrodes=localize_PREF(a_path,1,freq,smoothing,do_erp,1,load_old_stat_files)
% save('E:\ALIZEE\EPILEPSY\ALIZEE\EPILEPSY\Subjects_structures\LYONNEURO_2015_DUBb_PREF_structure_result','stat_electrodes')

% a_path='E:\ALIZEE\EPILEPSY\DATA_SAFE\LYON\LOCA\LYONNEURO_2015_PASj_PREF';
% stat_electrodes=localize_PREF(a_path,1,freq,smoothing,do_erp,1,load_old_stat_files)
% save('E:\ALIZEE\EPILEPSY\ALIZEE\EPILEPSY\Subjects_structures\LYONNEURO_2015_PASj_PREF_structure_result','stat_electrodes')
% 
% a_path='E:\ALIZEE\EPILEPSY\DATA_SAFE\LYON\LOCA\LYONNEURO_2015_SCAl_PREF';
% stat_electrodes=localize_PREF(a_path,1,freq,smoothing,do_erp,1,load_old_stat_files)
% save('E:\ALIZEE\EPILEPSY\ALIZEE\EPILEPSY\Subjects_structures\LYONNEURO_2015_SCAl_PREF_structure_result','stat_electrodes')

% a_path='E:\ALIZEE\EPILEPSY\DATA_SAFE\LYON\LOCA\LYONNEURO_2015_SELs_PREF';
% stat_electrodes=localize_PREF(a_path,1,freq,smoothing,do_erp,1,load_old_stat_files)
% save('E:\ALIZEE\EPILEPSY\ALIZEE\EPILEPSY\Subjects_structures\LYONNEURO_2015_SELs_PREF_structure_result','stat_electrodes')


%% =========================================================

% waiting for correct electrode names : 
% a_path='E:\ALIZEE\EPILEPSY\DATA_SAFE\LYON\LOCA\LYONNEURO_2014_DESJ_PREF';
% stat_electrodes=localize_PREF(a_path,1,freq,smoothing,do_erp,compute_gamma,load_old_stat_files)
% save('E:\ALIZEE\EPILEPSY\ALIZEE\EPILEPSY\Subjects_structures\LYONNEURO_2014_DESJ_PREF_structure_result','stat_electrodes')

% a_path='E:\ALIZEE\EPILEPSY\DATA_SAFE\LYON\LOCA\LYONNEURO_2014_MARB_PREF';
% stat_electrodes=localize_PREF(a_path,1,freq,smoothing,do_erp,1,load_old_stat_files)
% save('E:\ALIZEE\EPILEPSY\ALIZEE\EPILEPSY\Subjects_structures\LYONNEURO_2014_MARB_PREF_structure_result','stat_electrodes')

% a_path='E:\ALIZEE\EPILEPSY\DATA_SAFE\LYON\LOCA\LYONNEURO_2014_MONM_PREF';
% stat_electrodes=localize_PREF(a_path,1,freq,smoothing,do_erp,0,load_old_stat_files)
% save('E:\ALIZEE\EPILEPSY\ALIZEE\EPILEPSY\Subjects_structures\LYONNEURO_2014_MONM_PREF_structure_result','stat_electrodes')


%% POS and LOG DO NOT MATCH 

% (86 outliers over 120)
% 
% a_path='E:\ALIZEE\EPILEPSY\DATA_SAFE\LYON\LOCA\LYONNEURO_2014_PERR_PREF';
% stat_electrodes=localize_PREF(a_path,1,freq,smoothing,do_erp,1,load_old_stat_files)
% save('E:\ALIZEE\EPILEPSY\ALIZEE\EPILEPSY\Subjects_structures\LYONNEURO_2014_PERR_PREF_structure_result','stat_electrodes')
% 

% 75 outliers
% a_path='E:\ALIZEE\EPILEPSY\DATA_SAFE\LYON\LOCA\LYONNEURO_2014_REID_PREF';
% stat_electrodes=localize_PREF(a_path,1,freq,smoothing,do_erp,1,load_old_stat_files)
% save('E:\ALIZEE\EPILEPSY\ALIZEE\EPILEPSY\Subjects_structures\LYONNEURO_2014_REID_PREF_structure_result','stat_electrodes')

% a_path='E:\ALIZEE\EPILEPSY\DATA_SAFE\LYON\LOCA\LYONNEURO_2014_CHOk_PREF';
% stat_electrodes=localize_PREF(a_path,1,freq,smoothing,do_erp,1,load_old_stat_files)
% save('E:\ALIZEE\EPILEPSY\ALIZEE\EPILEPSY\Subjects_structures\LYONNEURO_2014_CHOk_PREF_structure_result','stat_electrodes')
% 

%%


