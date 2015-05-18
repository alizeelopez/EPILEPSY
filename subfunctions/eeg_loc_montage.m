function [v_neworder, v_neworder_rev, m_bipole]=eeg_loc_montage(ss_channel)

%load truc;
v_rt=[];
s_counter=1;
s_other=1;
v_neworder=[];
v_name=[];
for s_c=1:length(ss_channel.values)
    if (s_c == 69)
        a = 69;
    end
    a_string=ss_channel.values(s_c).value;
    a_string=deblank(a_string);
    a_string = strrep(a_string,'-',''); % ADDED TO AVOID PB WITH G-
    a=a_string;
    v_name=strvcat(v_name,a_string);
    v_num=find((a=='0') | (a=='1') | (a=='2') | (a=='3') | (a=='4') | (a=='5') | (a=='6') | (a=='7') | (a=='8') | (a=='9'));
    % to be a real intra site, a_string has to be of the form f'12 or a5 ... that is, all the last chars of a_string are numbers, and there is only one letter with, for the left elecs, a ' sign
    if isempty(v_num)
        s_test=0; % A. we want numbers in the channel name
    else
        s_test=((max(v_num)-min(v_num))==(length(v_num)-1)); % B. we want the numbers to be contiguous
        s_test=s_test*(max(v_num)==length(a_string)); % C. we want the numbers to be the last characters of the string
        %s_test=s_test*((min(v_num)==2)|(min(v_num)==3));
        
        if min(v_num)==2
            s_test=s_test*(isletter(a_string(1)));   % for any string like A + XXX, where A is a letter and XXX = a number string
            a_rt=a_string(1);
            s_id=str2num(a_string(v_num));
            a_name=a_string;
        end
        if min(v_num)==3
            s_test=s_test*(isletter(a_string(1))); % for any string like AM + XXX, where A is a letter and XXX = a number string, and M anything but a digit
            if (strcmp(a_string(1:2),'ET')) % special case "ET" + XXX
                a_rt=a_string(1:2);
                s_id=str2num(a_string(v_num));
                a_name=a_string;
            else
                if (strcmp(upper(a_string(1:2)),'AX'))
                    a_rt=a_string(1:2);
                    s_id=str2num(a_string(v_num));
                    a_name=a_string;
                else
                    s_test=s_test*(a_string(2)==''''); % special case " A' " + XXX, where A is a letter
                    a_rt=a_string(1:2);
                    s_id=str2num(a_string(v_num));
                    a_name=a_string;
                end
            end
        end
        
        if min(v_num)==4
            s_test=s_test*(isletter(a_string(1))); % for any string like AMM + XXX, where A is a letter and XXX = a number string, and M anything but a digit
            if (strcmp(a_string(1:2),'ET'))
                s_test=s_test*(a_string(3)==''''); % special case " ET' " + XXX
                a_rt=a_string(1:3);
                s_id=str2num(a_string(v_num));
                a_name=a_string;
            elseif (strcmp(a_string(1:3),'ECG'))    % special case " ECG " + XXX
                a_rt=a_string(1:3);
                s_id=str2num(a_string(v_num));
                a_name=a_string;
            else
                s_test = 0;
            end
        end
        if min(v_num)>4 % for any string like AMMM + XXX, where A is a letter and XXX = a number string, and M anything but a digit
            s_test = 0;
        end
    end % if isempty
    
    if s_test
        v_find=strmatch(a_rt,v_rt,'exact');
        if isempty(v_find)
            v_rt=strvcat(v_rt,a_rt); % that's a new electrode
            ss_elec(s_counter).a_rt=a_rt;
            ss_elec(s_counter).v_id=[s_id];
            ss_elec(s_counter).v_origid=[s_c]; % where was this elec in the original trc file (in the data)
            s_counter=s_counter+1;
        else
            s_count=v_find(1);
            ss_elec(s_count).v_id=[ss_elec(s_count).v_id s_id];
            ss_elec(s_count).v_origid=[ss_elec(s_count).v_origid s_c]; % where was this elec in the original trc file (in the data)
        end % if isempty
    else
        ss_other(s_other).a_rt=a_string;
        ss_other(s_other).v_id=[];
        ss_other(s_other).v_origid=[s_c]; % where was this elec in the original trc file (in the data)
        s_other=s_other+1;
    end % if (s_test)
end % for s_c


% the objective is to find, for each channel, from 1 to the length of ss_channel, as ordered in ss_channel, in the raw TRC file, the new order in the upcoming ella files
% now we want to reorder everything
s_total=0;
v_c=zeros(1,length(ss_channel.values));
m_bipole=[];
for s_c=1:(s_counter-1) % for each electrode, we want to find which s_c is the first
    if (s_c == 9)
        a = 69;
    end
    a_rt=ss_elec(s_c).a_rt;
    v_id=ss_elec(s_c).v_id;
    v_origid=ss_elec(s_c).v_origid;
    [v_y,v_i]=sort(v_id);
    % question, if elec v_origid(i) is the n-th of this group, what is the total ranking ?
    %  v_origid=v_origid(v_i);
    % what is the ranking of v_origid(i) within this group ?
    % it is the ranking of v_id(i)
    % the smallest value of v_id is v_id(v_i(1)) the second smallest is v_id(v_i(2))
    % then, which one is the smallest ? v_i(1)
    % what is the ranking of v_id(i) ?
    clear v_rank;

    for s_i=1:length(v_id)
        v_rank(s_i)=find(v_y==v_id(s_i))+s_total;
        %                 v_rank(s_i)=find(v_y==v_id(s_i))+s_total;
        
    end
    % v_rank(s_i) is the ranking of channel v_origid(s_i)
    v_c(v_origid)=v_rank;
    %v_c(s_i) is the new indice of the channel originally in position s_i in
    % the trc file.
    % what is the original indice of a channel of new indice s_i?
    
    % we also want to create the bipole txt file
    % that means we want within each electrode, find the site pairs of the type i i+1
    for s_b=2:length(v_y)
        if ((v_y(s_b)-v_y(s_b-1))==1) % then we have a bipole
            m_bipole=[m_bipole;[s_total+s_b s_total+s_b-1]];
        end
    end
    s_total=s_total+length(v_id);
end % for s_c

% at this point, the electrodes that have not beed reordered should have a
% zero in s_c
v_find_zero=find(v_c==0);
for s_c=1:length(v_find_zero)
    s_ii=v_find_zero(s_c);
    s_total=s_total+1;
    v_c(s_ii)=s_total;
end % for s_c

% now, what is v_neworder ? and what do we do with the other channels ?
% we also compute the inverse of v_neworder : which former electrode goes into new position i ? this is v_neworder_rev(i)
v_neworder=v_c; % this means that v_neworder(i) is where I want electrode i (among the visible ones in fff_raw) to be in my new organization.
for s_e=1:max(v_neworder)
    s_j=find(v_neworder==s_e);
    v_neworder_rev(s_e)=s_j;
end;
v_neworder_rev=v_neworder_rev(:);
