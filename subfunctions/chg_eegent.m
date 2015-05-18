function chg_eegent(a_eegfile,a_elecdat)
 
 % chg_eegent('D:\data_elan\gei2\exp_Gei_7_8_2003_spont','D:\implante\elec.dat');
 % chg_eegent('F:\epi_ab2\ab2_hal\exp_Halici_12_4_2005_ab2_veo','F:\elec.dat');
 % this script changes the eeg.ent to adjust the electrode numbers so that
 % they fit witht the elec.dat file in elan
 
 % first, back up the old eeg.ent
 %
 
 a_oldent=[a_eegfile '.eeg.ent.old'];
 a_bupent=[a_eegfile '.eeg.ent.bup'];
 a_newent=[a_eegfile '.eeg.ent'];
 
 % backup the eeg.ent file
 copyfile(a_newent,a_oldent);
 copyfile(a_newent,a_bupent);
 
 % read the elec.dat
 [v_id v_x v_y a_channame]=textread(a_elecdat,'%d %f %f %s');
 
 a_channame=char(a_channame);
 % then read the eeg.ent file
 f_old=fopen(a_oldent,'r');
 f_new=fopen(a_newent,'w');
 
 for s_i=1:10
     a_line=fgetl(f_old);
     fprintf(f_new,'%s\n',a_line);
 end;
 
 % how many channels ?
 s_nbchannel=str2num(a_line);
 for s_c=1:s_nbchannel-2
     % read the name of the actual data channels
     s_indice=-1;
     a_line=fgetl(f_old);
     % the name is something like Cz.10;
     s_point=find(a_line=='.');
     a_elecname=a_line(1:s_point(1)-1);
     a_elecname=upper(a_elecname);
     % now we have to find it in the list of elec.dat
     v_i=strmatch(a_elecname,a_channame);
     if (~isempty(v_i))
         for s_i=1:length(v_i)
             s_norm=abs(v_x(v_i(s_i)))+abs(v_y(v_i(s_i)));
             if (s_norm==0)
                 s_indice=v_id(v_i(s_i));
                 break;
             end;
         end;
     end;
     if (s_indice==-1)
         disp([a_elecname ' not found in elec.dat']);
     end;
     a_newline=[a_line(1:s_point(1)) int2str(s_indice)];
     fprintf(f_new,'%s\n',a_newline);
 end; % for s_c
 
 while (a_line~=-1)
     a_line=fgetl(f_old);
     if (a_line==-1)
         % do nothing
     else
         fprintf(f_new,'%s\n',a_line);
     end;
 end;
 
 fclose(f_new);
 fclose(f_old);