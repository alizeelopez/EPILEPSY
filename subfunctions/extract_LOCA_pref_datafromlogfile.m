function [a_code v_ttime pref_score stim_id tmp]=extract_LOCA_pref_datafromlogfile(a_logfile)
 
%   New nomenclature:
%     a_base = path;
%     a_patient = pat;
    
%   Creating a list "a_lognames" containing the conditions listed in the condition txt file:


    
%   For all the conditions:

    %   Naming log the name of the current condition:
%         log = char(a_lognames(nummatrix)); 
    %   Naming "a_logfile" the log file corresponding to the current condition and "a_logfile" as opening it:
%         a_logfile = [a_base filesep a_patient filesep 'behavior' filesep log  '.log'];
        f_log=fopen(a_logfile,'r');
    %   Naming "a_trashfile" the trash file corresponding to the current condition and "f_trash" as opening it:
        a_trashfile=strrep((a_logfile),'.log','.trash');
        f_trash=fopen(a_trashfile,'w');
    
    %
        s_do_responses_send_pulses=1;        
     
    %   Effacer les 4 premi�res lignes du fichier log:   
        for n_lines=1:4    
            a_line=fgetl(f_log);
        end
    %   Cr�ation de deux listes vides "a_code_ok" et "v_ttime_ok":    
        a_code_ok=[];
        v_ttime_ok=[]

        if (s_do_responses_send_pulses==0)
            a_repfile=strrep(lower(a_logfile),'.log','.trashrep');
            f_rep=fopen(a_repfile,'w');
        end
        if (~isempty(findstr(lower(a_line),'subject')))
            s_sub=1;
        else
            s_sub=0;
        end

        a_line=fgetl(f_log);
        a_line=fgetl(f_log);

        while ((~isempty(a_line)) & (a_line~=-1))
            if (s_do_responses_send_pulses==0)
                if (isempty(findstr(a_line,'Response')))
                    fprintf(f_trash,'%s\n',a_line);
                else
                %   well, we have to put the response somewhere    
                    fprintf(f_rep,'%s\n',a_line);
                end
                a_line=fgetl(f_log);
            else
                fprintf(f_trash,'%s\n',a_line);
                a_line=fgetl(f_log);
            end
        end
        fclose(f_trash);
        fclose(f_log);
        [a_trial a_evttype a_code v_ttime v_tttime v_unctime v_duration v_uncdur v_reqtime a_reddur]=textread(a_trashfile,'%s %s %s %d %d %d %d %d %d %s');
   
        %correct v_ttime
        v_ttime=(v_ttime-v_ttime(1))/10;%v_ttime starts at 0 and is now in milliseconds

        %on a les codes et les timing associ�s en ms !

        %etape suivante: on ne garde que les codes qui nous interessent +
        %les stim ID
        s_count=1;k=0;tmp=[];
        %in stim_id = we will store the picture number displayed
        %in tmp we will store the index corresponding to a stimulus
        %presentation in the logfile(useful for later)
        %pre_score = rating
        for s_i=1:length(a_evttype)
            uch=char(a_code(s_i));
            try
            if uch(1)=='N'      %stim IF begins by letter 'N'
                stim_id(s_count)=str2double(uch(2:end))
             s_count=s_count+1;
             tmp=[tmp s_i],
            end
            catch
            end
                
        try
            if strcmp(uch(1:4),'CF1=')
                k=k+1
                pref_score(k)=str2num(uch(5:end))
            end
        catch
        end
        end



       
        
        
        
            
       