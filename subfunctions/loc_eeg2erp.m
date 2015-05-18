function loc_env2plot(a_eegname,a_posname,v_code,a_code,v_window_ms,s_maxsite_per_elec)

% loc_env2bar('D:\data_elan\LT_PENS_TEST\LT_20NOV08_PENS_f50f150_ds8_sm1000.eeg','D:\data_elan\LT_PENS_TEST\lt_20nov08_pens_ds8.pos',[10 20 30 40],strvcat('faible','normal','fort','attente'),[-500 500],14)

% if you want to do it on fdr-corrected files
s_nbsite = s_maxsite_per_elec;
[fff_data_mono,fff_data_bipo,m_pos] = loc_cat2ella_monobipo(a_eegname,a_posname,v_code,v_window_ms);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A. MONOPOLAR   %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% the first thing is to extract the data.
% we will loop over all channels

fff_data = fff_data_mono;

m_bigdata = fff_data.value.getme;

ss_channel = fff_data.dim(2).values;
s_nbchannel = length(ss_channel);

% now, identify the electrode names.
% for each new channel, identify the electrode on which it is, then the
% position on the electrode. and the indice in the m_scf matrix

v_rt=[];
s_counter=1;
s_other=1;
v_neworder=[];
v_name=[];
ss_elec = [];
for s_c=1:(s_nbchannel-2)
    a_string=ss_channel(s_c).value;
    a_string=deblank(a_string);
    a=a_string;
    v_name=strvcat(v_name,a_string);
    v_num=find((a=='0') | (a=='1') | (a=='2') | (a=='3') | (a=='4') | (a=='5') | (a=='6') | (a=='7') | (a=='8') | (a=='9'));
    % to be a real intra site, a_string has to be of the form f'12 or a5 ... that is, all the last chars of a_string are numbers, and there is only one letter with, for the left elecs, a ' sign
    
    a_rt=a_string(1:min(v_num)-1);
    if (length(v_num)>1)
        if (max(v_num)==length(a_string))
            s_tool = length(v_num);
        else
            v_tool = diff(v_num);
            s_tool = min(find(v_tool>1));
        end;
        v_num2 = v_num(1:s_tool);
    else
        v_num2 = v_num;
    end;
    
    s_id=str2num(a_string(v_num2));
    
    
    v_find=strmatch(a_rt,v_rt,'exact');   
    if (isempty(v_find))
        v_rt=strvcat(v_rt,a_rt); % that's a new electrode
        ss_elec(s_counter).a_rt=a_rt;
        ss_elec(s_counter).v_id=[s_id];
        ss_elec(s_counter).v_origid=[s_c]; % where was this elec in the original file (in the data)     
        s_counter=s_counter+1;
    else
        s_count=v_find(1);
        ss_elec(s_count).v_id=[ss_elec(s_count).v_id s_id];
        ss_elec(s_count).v_origid=[ss_elec(s_count).v_origid s_c]; % where was this elec in the original file (in the data)     
    end; % if isempty
end; % for s_c

% end of this difficult exercice



% SECOND PART - MORE SPECIFIC

% the loop is on the electrodes. every time you have finished a series of 3
% electrodes, you start a new figure
s_nbelectrode_per_figure = 3;
s_figure = 0;

m_color = [[0 0 0.6];[0.6 0 0];[0 0.6 0]];
m_color = repmat(m_color,8,1);

for s_c = 1:length(ss_elec)
%    for s_c = 1:3
    
    a_rt = ss_elec(s_c).a_rt;
    v_id = ss_elec(s_c).v_id
    v_origid = ss_elec(s_c).v_origid
    
    s_pos = rem(s_c-1,s_nbelectrode_per_figure);
    if (s_pos == 0)
        if (s_figure>0)
            %print(gcf,'-dtiff','-r600',a_tifname);
            print(gcf,'-djpeg100','-r200',a_tifname);
            close(gcf);
        end;    
        
        s_figure = s_figure+1;
        % start a new figure
        figure('Color','w');
        set(gcf,'PaperUnits','centimeters');
        xSize = 30; ySize = 21;
        xLeft = (30-xSize)/2; yTop = (21-ySize)/2;
        set(gcf,'PaperPosition',[0*xLeft yTop xSize ySize]);
        sr = get (0,'ScreenSize');
        set(gcf,'Position',[sr(1) sr(2) sr(3) sr(4)]);
        set(gcf,'PaperOrientation','Portrait');
        a_tifname = [strrep(a_eegname,'.eeg','') '_erp_mono_f' int2str(s_figure) '.jpg'];
    end;
    % as input, I take up to three electrodes, and for each electrode and site
    % along that electrode, I want to plot a certain number of plots. So I need
    % the data for each condition, with a time-axis.
    
    % 
    % draw electrode, with 20 sites, equally spaced.
    
    
    h_a = axes('position',[0 0 1 1],'visible','off');
    hold on;
    set(h_a,'XLim',[0 1],'YLim',[0 1]);
    
    s_He_pct = 80;
    s_He = s_He_pct/100;
    s_Le_pct = 1;
    s_Le = s_Le_pct/100;
    
    
    s_width_pct = 23; % axis width . SPECIFIC 
    s_inter_plot_pct = 2;
    s_empty_pct = 20; %  % of empty space left of left curve
    s_xelec_pct = 10;
    
    s_x = s_xelec_pct/300 + s_pos/3;
    s_y = 0.1;
    s_eps = s_He/(2*s_nbsite);
    
    h = patch([s_x s_x+s_Le s_x+s_Le s_x],[s_y s_y s_y+s_He s_y+s_He],'w');
    text(s_x,s_y+s_He+s_eps/2,a_rt,'FontWeight','Bold','FontSize',12);
    
    s_width_tool = 0.333 - (s_empty_pct/300);
    s_ncol = ceil(length(v_code)/3);
    
    s_tx = 0;
    s_ty = 8;
    for s_i = 1:length(v_code)
        s_ty = s_ty - 1.5;
        if (s_ty < 2.5 )
            s_ty = 6.5;
            s_tx = s_tx + s_width_tool/s_ncol;
        end;
        text((s_empty_pct/300 + s_pos/3)+s_tx,s_y+s_He+s_ty*s_eps/2,[num2str(s_i) '. ' a_code(s_i,:)],'Color',m_color(s_i,:),'FontWeight','Bold','FontSize',10);
        
    end;
    
    % first we draw the axis
    for s_e = 1:s_nbsite
        s_ye = s_y + s_He -2*s_eps*s_e;
        axes(h_a);
        patch([s_x s_x+s_Le s_x+s_Le s_x],[s_ye s_ye s_ye+s_eps s_ye+s_eps],'k');
        text(s_x-0.02,s_ye+s_eps/2,int2str(s_e));
        
        
        % do we have data to show for that site ?
        v_f = find(v_id==s_e);
        if (~isempty(v_f))
            
            % this site corresponds to id v_f(1) in this elec
            % which should correspond in the original data to 
            s_orig_id = v_origid(v_f(1));
            % then we can extract the data
            m_data = m_bigdata(:,s_orig_id,:);
            m_data= squeeze(m_data);
            
            % then we loop across conditions
            for s_p = 1:length(v_code)
                % for each event code, calculate mean and sem across trials
                % of that condition
                v_f = find(m_pos(:,2)==v_code(s_p)); % the events for that condition
                if (isempty(v_f))
                    disp(['ATTENTION, PLEASE : NO EVENT WITH TYPE = ' int2str(v_code(s_p)) ]);
                    m_erp(:,s_p) = zeros(size(m_erp,1),1);
                    m_lim(:,s_p) = zeros(size(m_erp,1),1);
                else    
                    m_data_se = m_data(:,v_f);
                    %m_data_se = (m_data_se - 1000)/10; % in %
                    v_erp = mean(m_data_se,2);
                    v_std = std(m_data_se,0,2);
                    v_sem = v_std/sqrt(size(m_data_se,2)); % ck
                    v_lim = 1.96*v_sem;
                    
                    
                    m_erp(:,s_p) = v_erp(:);
                    m_lim(:,s_p) = v_lim(:);
                    v_t = linspace(v_window_ms(1),v_window_ms(2),length(v_erp));
                end;
            end;
            
            m_erp = m_erp - mean(m_erp(1,:)); % we remove the same value to all the erp's, so that it is more or less centered on zero at the first sample
            
            s_max1 = max(max(abs(m_erp+m_lim)));
            s_max2 = max(max(abs(m_erp-m_lim)));
            s_max = max(s_max1,s_max2);
            
            text(s_x+0.015,s_ye+s_eps/2,int2str(s_max),'FontWeight','Bold','FontSize',6,'Color','m');
            
            
            axes('Position',[(s_empty_pct/300 + s_pos/3) s_ye-1.5*s_eps s_width_tool 4*s_eps],'Visible','off'); % s_ye-s_eps/2
            hold on;
            
            
            for s_p = 1:length(v_code)
                s_add = (max(v_t) - min(v_t))*floor((s_p-1)/3);
                s_tmin = min(v_t) + s_add;
                s_tmax = max(v_t) + s_add;
                s_tzero = s_add;
                s_t500 = 500 + s_add;
                
                fill([s_tmin s_tzero s_tzero s_tmin],[-1 -1 1 1],[0.8 0.8 0.8],'EdgeColor','none','FaceColor',[0.8 0.8 0.8]);
                h = fill([s_tzero s_t500 s_t500 s_tzero],[-1 -1 1 1],[0.9 0.9 0.9],'EdgeColor','none','FaceColor',[0.9 0.9 0.9]);
                fill([s_t500 s_tmax s_tmax s_t500],[-1 -1 1 1],[1 1 1],'EdgeColor','none','FaceColor',[1 1 1]);
                
                plot(v_t + s_add,0*v_t,'Color','k','LineWidth',0.5,'LineStyle','-');
                
            end;
            
            
            axis([min(v_t) (max(v_t) + s_add) -1.5 1.5]);
            axis off;
            
            
            
            
        end;
        
    end;
    
    % then we draw the curves
    for s_e = 1:s_nbsite
        s_ye = s_y + s_He -2*s_eps*s_e;
        axes(h_a);
        
        
        % do we have data to show for that site ?
        v_f = find(v_id==s_e);
        if (~isempty(v_f))
            
            % this site corresponds to id v_f(1) in this elec
            % which should correspond in the original data to 
            s_orig_id = v_origid(v_f(1));
            % then we can extract the data
            m_data = m_bigdata(:,s_orig_id,:);
            m_data= squeeze(m_data);
            
            % then we loop across conditions
            for s_p = 1:length(v_code)
                % for each event code, calculate mean and sem across trials
                % of that condition
                v_f = find(m_pos(:,2)==v_code(s_p)); % the events for that condition
                m_data_se = m_data(:,v_f);
                m_data_se = (m_data_se - 1000)/10; % in %
                v_erp = mean(m_data_se,2);
                v_std = std(m_data_se,0,2);
                v_sem = v_std/sqrt(size(m_data_se,2)); % ck
                v_lim = 1.96*v_sem;
                
                
                m_erp(:,s_p) = v_erp(:);
                m_lim(:,s_p) = v_lim(:);
                v_t = linspace(v_window_ms(1),v_window_ms(2),length(v_erp));
            end;
            
            m_erp = m_erp - mean(m_erp(1,:)); % we remove the same value to all the erp's, so that it is more or less centered on zero at the first sample
            
            s_max1 = max(max(abs(m_erp+m_lim)));
            s_max2 = max(max(abs(m_erp-m_lim)));
            s_max = max(s_max1,s_max2);
            
            
            
            axes('Position',[(s_empty_pct/300 + s_pos/3) s_ye-1.5*s_eps s_width_tool 4*s_eps],'Visible','off'); % s_ye-s_eps/2
            hold on;
            
            
            for s_p = 1:length(v_code)
                v_erp = m_erp(:,s_p)/s_max;
                v_lim = m_lim(:,s_p)/s_max;
                
                s_add = (max(v_t) - min(v_t))*floor((s_p-1)/3);
                
                
                plot(v_t + s_add,v_erp,'Color',m_color(s_p,:),'LineWidth',1.5,'LineStyle','-');
                plot(v_t + s_add,v_erp+v_lim,'Color',m_color(s_p,:),'LineWidth',1,'LineStyle','-');
                plot(v_t + s_add,v_erp-v_lim,'Color',m_color(s_p,:),'LineWidth',1,'LineStyle','-');
                
                plot(v_t + s_add,0*v_t,'Color','k','LineWidth',0.5,'LineStyle','-');
                
            end;
            
            axis([min(v_t) (max(v_t) + s_add) -1.5 1.5]);
            axis off;
            
            
            
            
        end;
        
    end;
    
end;




% to end with, print last figure
%print(gcf,'-dtiff','-r200',a_tifname);
print(gcf,'-djpeg100','-r200',a_tifname);
close(gcf);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A. BIPOPOLAR   %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% the first thing is to extract the data.
% we will loop over all channels

fff_data = fff_data_bipo;

m_bigdata = fff_data.value.getme;

ss_channel = fff_data.dim(2).values;
s_nbchannel = length(ss_channel);

% now, identify the electrode names.
% for each new channel, identify the electrode on which it is, then the
% position on the electrode. and the indice in the m_scf matrix

v_rt=[];
s_counter=1;
s_other=1;
v_neworder=[];
v_name=[];
ss_elec = [];
for s_c=1:(s_nbchannel-2)
    a_string=ss_channel(s_c).value;
    a_string=deblank(a_string);
    a=a_string;
    v_name=strvcat(v_name,a_string);
    v_num=find((a=='0') | (a=='1') | (a=='2') | (a=='3') | (a=='4') | (a=='5') | (a=='6') | (a=='7') | (a=='8') | (a=='9'));
    % to be a real intra site, a_string has to be of the form f'12 or a5 ... that is, all the last chars of a_string are numbers, and there is only one letter with, for the left elecs, a ' sign
    
    a_rt=a_string(1:min(v_num)-1);
    if (length(v_num)>1)
        if (max(v_num)==length(a_string))
            s_tool = length(v_num);
        else
            v_tool = diff(v_num);
            s_tool = min(find(v_tool>1));
        end;
        v_num2 = v_num(1:s_tool);
    else
        v_num2 = v_num;
    end;
    
    s_id=str2num(a_string(v_num2));
    
    
    v_find=strmatch(a_rt,v_rt,'exact');   
    if (isempty(v_find))
        v_rt=strvcat(v_rt,a_rt); % that's a new electrode
        ss_elec(s_counter).a_rt=a_rt;
        ss_elec(s_counter).v_id=[s_id];
        ss_elec(s_counter).v_origid=[s_c]; % where was this elec in the original file (in the data)     
        s_counter=s_counter+1;
    else
        s_count=v_find(1);
        ss_elec(s_count).v_id=[ss_elec(s_count).v_id s_id];
        ss_elec(s_count).v_origid=[ss_elec(s_count).v_origid s_c]; % where was this elec in the original file (in the data)     
    end; % if isempty
end; % for s_c

% end of this difficult exercice



% SECOND PART - MORE SPECIFIC

% the loop is on the electrodes. every time you have finished a series of 3
% electrodes, you start a new figure
s_nbelectrode_per_figure = 3;
s_figure = 0;

m_color = [[0 0 0.6];[0.6 0 0];[0 0.6 0]];
m_color = repmat(m_color,8,1);

for s_c = 1:length(ss_elec)
%    for s_c = 1:3
    
    a_rt = ss_elec(s_c).a_rt;
    v_id = ss_elec(s_c).v_id
    v_origid = ss_elec(s_c).v_origid
    
    s_pos = rem(s_c-1,s_nbelectrode_per_figure);
    if (s_pos == 0)
        if (s_figure>0)
            print(gcf,'-dill','-r300',a_tifname2);
            print(gcf,'-djpeg100','-r200',a_tifname);
            close(gcf);
        end;    
        
        s_figure = s_figure+1;
        % start a new figure
        figure('Color','w');
        set(gcf,'PaperUnits','centimeters');
        xSize = 30; ySize = 21;
        xLeft = (30-xSize)/2; yTop = (21-ySize)/2;
        set(gcf,'PaperPosition',[0*xLeft yTop xSize ySize]);
        sr = get (0,'ScreenSize');
        set(gcf,'Position',[sr(1) sr(2) sr(3) sr(4)]);
        set(gcf,'PaperOrientation','Portrait');
        a_tifname = [strrep(a_eegname,'.eeg','') '_erp_bipo_f' int2str(s_figure) '.jpg'];
                a_tifname2 = [strrep(a_eegname,'.eeg','') '_erp_bipo_f' int2str(s_figure) '.ai'];
    end;
    % as input, I take up to three electrodes, and for each electrode and site
    % along that electrode, I want to plot a certain number of plots. So I need
    % the data for each condition, with a time-axis.
    
    % 
    % draw electrode, with 20 sites, equally spaced.
    
    
    h_a = axes('position',[0 0 1 1],'visible','off');
    hold on;
    set(h_a,'XLim',[0 1],'YLim',[0 1]);
    
    s_He_pct = 80;
    s_He = s_He_pct/100;
    s_Le_pct = 1;
    s_Le = s_Le_pct/100;
    
    
    s_width_pct = 23; % axis width . SPECIFIC 
    s_inter_plot_pct = 2;
    s_empty_pct = 20; %  % of empty space left of left curve
    s_xelec_pct = 10;
    
    s_x = s_xelec_pct/300 + s_pos/3;
    s_y = 0.1;
    s_eps = s_He/(2*s_nbsite);
    
    h = patch([s_x s_x+s_Le s_x+s_Le s_x],[s_y s_y s_y+s_He s_y+s_He],'w');
    text(s_x,s_y+s_He+s_eps/2,a_rt,'FontWeight','Bold','FontSize',12);
    
    s_width_tool = 0.2 - (s_empty_pct/300);%set it back to 0.333 !!
    s_ncol = ceil(length(v_code)/3);
    
    s_tx = 0;
    s_ty = 8;
    for s_i = 1:length(v_code)
        s_ty = s_ty - 1.5;
        if (s_ty < 2.5 )
            s_ty = 6.5;
            s_tx = s_tx + s_width_tool/s_ncol;
        end;
        text((s_empty_pct/300 + s_pos/3)+s_tx,s_y+s_He+s_ty*s_eps/2,[num2str(s_i) '. ' a_code(s_i,:)],'Color',m_color(s_i,:),'FontWeight','Bold','FontSize',10);
        
    end;
    
    % first we draw the axis
    for s_e = 1:s_nbsite
        s_ye = s_y + s_He -2*s_eps*s_e;
        axes(h_a);
        patch([s_x s_x+s_Le s_x+s_Le s_x],[s_ye s_ye s_ye+s_eps s_ye+s_eps],'k');
        text(s_x-0.02,s_ye+s_eps/2,int2str(s_e));
        
        
        % do we have data to show for that site ?
        v_f = find(v_id==s_e);
        if (~isempty(v_f))
            
            % this site corresponds to id v_f(1) in this elec
            % which should correspond in the original data to 
            s_orig_id = v_origid(v_f(1));
            % then we can extract the data
            m_data = m_bigdata(:,s_orig_id,:);
            m_data= squeeze(m_data);
            
            % then we loop across conditions
            for s_p = 1:length(v_code)
                % for each event code, calculate mean and sem across trials
                % of that condition
                v_f = find(m_pos(:,2)==v_code(s_p)); % the events for that condition
                if (isempty(v_f))
                    disp(['ATTENTION, PLEASE : NO EVENT WITH TYPE = ' int2str(v_code(s_p)) ]);
                    m_erp(:,s_p) = zeros(size(m_erp,1),1);
                    m_lim(:,s_p) = zeros(size(m_erp,1),1);
                else    
                    m_data_se = m_data(:,v_f);
                    %m_data_se = (m_data_se - 1000)/10; % in %
                    v_erp = mean(m_data_se,2);
                    v_std = std(m_data_se,0,2);
                    v_sem = v_std/sqrt(size(m_data_se,2)); % ck
                    v_lim = 1.96*v_sem;
                    
                    
                    m_erp(:,s_p) = v_erp(:);
                    m_lim(:,s_p) = v_lim(:);
                    v_t = linspace(v_window_ms(1),v_window_ms(2),length(v_erp));
                end;
            end;
            
            m_erp = m_erp - mean(m_erp(1,:)); % we remove the same value to all the erp's, so that it is more or less centered on zero at the first sample
            
            s_max1 = max(max(abs(m_erp+m_lim)));
            s_max2 = max(max(abs(m_erp-m_lim)));
            s_max = max(s_max1,s_max2);
            
            text(s_x+0.015,s_ye+s_eps/2,int2str(s_max),'FontWeight','Bold','FontSize',6,'Color','m');
            
            
            axes('Position',[(s_empty_pct/300 + s_pos/3) s_ye-1.5*s_eps s_width_tool 4*s_eps],'Visible','off'); % s_ye-s_eps/2
            hold on;
            
            
            for s_p = 1:length(v_code)
                s_add = (max(v_t) - min(v_t))*floor((s_p-1)/3);
                s_tmin = min(v_t) + s_add;
                s_tmax = max(v_t) + s_add;
                s_tzero = s_add;
                s_t500 = 500 + s_add;
                
                fill([s_tmin s_tzero s_tzero s_tmin],[-1 -1 1 1],[0.8 0.8 0.8],'EdgeColor','none','FaceColor',[0.8 0.8 0.8]);
                h = fill([s_tzero s_t500 s_t500 s_tzero],[-1 -1 1 1],[0.9 0.9 0.9],'EdgeColor','none','FaceColor',[0.9 0.9 0.9]);
                fill([s_t500 s_tmax s_tmax s_t500],[-1 -1 1 1],[1 1 1],'EdgeColor','none','FaceColor',[1 1 1]);
                
                plot(v_t + s_add,0*v_t,'Color','k','LineWidth',0.5,'LineStyle','-');
                
            end;
            
            
            axis([min(v_t) (max(v_t) + s_add) -1.5 1.5]);
            axis off;
            
            
            
            
        end;
        
    end;
    
    % then we draw the curves
    for s_e = 1:s_nbsite
        s_ye = s_y + s_He -2*s_eps*s_e;
        axes(h_a);
        
        
        % do we have data to show for that site ?
        v_f = find(v_id==s_e);
        if (~isempty(v_f))
            
            % this site corresponds to id v_f(1) in this elec
            % which should correspond in the original data to 
            s_orig_id = v_origid(v_f(1));
            % then we can extract the data
            m_data = m_bigdata(:,s_orig_id,:);
            m_data= squeeze(m_data);
            
            % then we loop across conditions
            for s_p = 1:length(v_code)
                % for each event code, calculate mean and sem across trials
                % of that condition
                v_f = find(m_pos(:,2)==v_code(s_p)); % the events for that condition
                m_data_se = m_data(:,v_f);
                m_data_se = (m_data_se - 1000)/10; % in %
                v_erp = mean(m_data_se,2);
                v_std = std(m_data_se,0,2);
                v_sem = v_std/sqrt(size(m_data_se,2)); % ck
                v_lim = 1.96*v_sem;
                
                
                m_erp(:,s_p) = v_erp(:);
                m_lim(:,s_p) = v_lim(:);
                v_t = linspace(v_window_ms(1),v_window_ms(2),length(v_erp));
            end;
            
            m_erp = m_erp - mean(m_erp(1,:)); % we remove the same value to all the erp's, so that it is more or less centered on zero at the first sample
            
            s_max1 = max(max(abs(m_erp+m_lim)));
            s_max2 = max(max(abs(m_erp-m_lim)));
            s_max = max(s_max1,s_max2);
            
            
            
            axes('Position',[(s_empty_pct/300 + s_pos/3) s_ye-1.5*s_eps s_width_tool 4*s_eps],'Visible','off'); % s_ye-s_eps/2
            hold on;
            
            
            for s_p = 1:length(v_code)
                v_erp = m_erp(:,s_p)/s_max;
                v_lim = m_lim(:,s_p)/s_max;
                
                s_add = (max(v_t) - min(v_t))*floor((s_p-1)/3);
                
                
                plot(v_t + s_add,v_erp,'Color',m_color(s_p,:),'LineWidth',1.5,'LineStyle','-');
                plot(v_t + s_add,v_erp+v_lim,'Color',m_color(s_p,:),'LineWidth',1,'LineStyle','-');
                plot(v_t + s_add,v_erp-v_lim,'Color',m_color(s_p,:),'LineWidth',1,'LineStyle','-');
                
                plot(v_t + s_add,0*v_t,'Color','k','LineWidth',0.5,'LineStyle','-');
                
            end;
            
            axis([min(v_t) (max(v_t) + s_add) -1.5 1.5]);
            axis off;
            
            
            
            
        end;
        
    end;
    
end;




% to end with, print last figure
print(gcf,'-dill','-r300',a_tifname2);
print(gcf,'-djpeg100','-r200',a_tifname);
close(gcf);







function [fff_data_mono,fff_data_bipo,m_pos] = loc_cat2ella_monobipo(a_eegname,a_posname,v_code,v_window_ms)

% this is directly from cat2ella

a_eegfile = a_eegname;
a_entfile=[a_eegname '.ent'];

% Step 1 : read the header file .ent
f_ent=fopen(a_entfile,'r'); % this is the txt file that contains the hdr in the ent file

a_version=fgetl(f_ent);
a_history=fgetl(f_ent);
a_patient=fgetl(f_ent);

a_line=fgetl(f_ent);
a_line=fgetl(f_ent);
a_line=fgetl(f_ent);
a_line=fgetl(f_ent);
a_line=fgetl(f_ent);

a_sampleperiod_s=fgetl(f_ent);
s_sampleperiod_s=str2num(a_sampleperiod_s);
s_fs=round(1/s_sampleperiod_s); % sampling frequency in Hz

v_window_sam = round(s_fs*v_window_ms/1000);

%s_nbsample_bef,s_nbsample_aft;

a_nbchannel=fgetl(f_ent); % number of channels including techincal stuff
s_nbchannel=str2num(a_nbchannel);

a_elecname=[];
for s_c=1:s_nbchannel-2
    % read the name of the actual data channels
    a_line=fgetl(f_ent);
    % the name is something like Cz.10;
    s_point=find(a_line=='.');
    a_elecname=strvcat(a_elecname,a_line(1:s_point-1));
    v_id(s_c)=str2num(a_line(s_point+1:length(a_line)));
end; % for s_c
a_line=fgetl(f_ent);
a_line=fgetl(f_ent);

a_elecname=strvcat(a_elecname,'Num1');
a_elecname=strvcat(a_elecname,'Num2');

% what's the type of each channel
a_electype=[];
for s_c=1:s_nbchannel-2
    % read the name of the actual data channels
    a_line=fgetl(f_ent);
    a_electype=strvcat(a_electype,a_line);
end; % for s_c
a_electype=strvcat(a_electype,'dateur');
a_electype=strvcat(a_electype,'event');
a_line=fgetl(f_ent);
a_line=fgetl(f_ent);



% general information 
fff_data_mono.name=['data from Elan, version ' a_version];
fff_data_mono.exam=a_patient;
fff_data_mono.history=a_history;
fff_data_mono.type.name='sample/channel/event';

% general information 
fff_data_bipo.name=['data from Elan, version ' a_version];
fff_data_bipo.exam=a_patient;
fff_data_bipo.history=a_history;
fff_data_bipo.type.name='sample/channel/event';


% first dimension = sample
clear ss_h
ss_h.id='time';
ss_h.unit='ms';
ss_h.type='sample';
ss_h.filter=[1 1 1 1];

s_sampling_freq=s_fs; % 

ss_h.values(1).value=v_window_sam(1); % this means that the first sample in this file is 400 samples before the event of reference
ss_h.values(1).authorized=1; % always 1 for samples
ss_h.values(1).description.id=ss_h.values(1).value*(1000/s_fs); % the equivalent latency in ms.
ss_h.values(1).description.sampling_freq=s_fs; % the equivalence sample-ms in Hz

ss_h.values(2).value=v_window_sam(2); % this means that the last sample in this file is 600 samples before the event of reference
ss_h.values(2).authorized=1; % always 1 for samples
ss_h.values(2).description.id=ss_h.values(2).value*(1000/s_fs); % the equivalent latency in ms.
ss_h.values(2).description.sampling_freq=s_fs; % the equivalence sample-ms in Hz

ss_h.values(3).value=v_window_sam(1); % this means that the first selected sample in this file is 200 samples before the event of reference
ss_h.values(3).authorized=1; % always 1 for samples
ss_h.values(3).description.id=ss_h.values(3).value*(1000/s_fs); % the equivalent latency in ms.
ss_h.values(3).description.sampling_freq=s_fs; % the equivalence sample-ms in Hz

ss_h.values(4).value=v_window_sam(2); % this means that the last selected sample is 200 samples after the event of reference
ss_h.values(4).authorized=1; % always 1 for samples
ss_h.values(4).description.id=ss_h.values(4).value*(1000/s_fs); % the equivalent latency in ms.
ss_h.values(4).description.sampling_freq=s_fs; % the equivalence sample-ms in Hz

fff_data_mono.dim(1)=ss_h;
fff_data_bipo.dim(1)=ss_h;


% dim 2: channels - monopolar
clear ss_h
ss_h.id='channels';
ss_h.unit='none';
ss_h.type='catego';
ss_h.filter=[ones(1,s_nbchannel-2) 0 0];

for s_i=1:s_nbchannel
    ss_h.values(s_i).value=a_elecname(s_i,:);
    ss_h.values(s_i).authorized=1;
    ss_h.values(s_i).description.id=a_elecname(s_i,:);
    ss_h.values(s_i).description.positive_input_label=a_elecname(s_i,:);
    ss_h.values(s_i).description.negative_input_label='G2';
end; % for s_i
ss_h.values(s_nbchannel-1).authorized=0;
ss_h.values(s_nbchannel).authorized=0;
fff_data_mono.dim(2)=ss_h;

ss_channel = ss_h;
[v_neworder, v_neworder_rev, m_bipole]=loc2_montage(ss_channel);
m_bipole_origid = v_neworder_rev(m_bipole); % for each bipole, specifies the indice of the big and litle brothers according to indexing in ss_channel and eeg file
% CAREFUL !!
v_bigsister = m_bipole(:,1);
v_litsister = m_bipole(:,2);
v_bigsister_o = m_bipole_origid(:,1);
v_litsister_o = m_bipole_origid(:,2);





% dim 2: channels - bipolar
% we want one field par bipole
s_nbbipole=length(v_bigsister);

clear ss_h
ss_h.id='channels';
ss_h.unit='none';
ss_h.type='catego';
ss_h.filter=[ones(1,s_nbbipole)];

for s_i=1:s_nbbipole
    ss_h.values(s_i).value=a_elecname(v_bigsister_o(s_i),:);
    ss_h.values(s_i).authorized=1;
    ss_h.values(s_i).description.id=a_elecname(v_bigsister_o(s_i),:);
    ss_h.values(s_i).description.positive_input_label=a_elecname(v_bigsister_o(s_i),:);
    ss_h.values(s_i).description.negative_input_label='G2';
end; % for s_i
fff_data_bipo.dim(2)=ss_h;



% dim 3: events
% we have to read the pos file first
m_pos = load(a_posname);
v_sub = zeros(size(m_pos,1),1); % ck
for s_i = 1:length(v_code)
    s_id = v_code(s_i);
    v_sub = v_sub + (m_pos(:,2)==s_id);
end;
v_f = find(v_sub); % only the rows with an event code as in v_code
m_pos = m_pos(v_f,:);

s_nbevent=size(m_pos,1);

clear ss_h
ss_h.id='([ALL])';
ss_h.unit='none';
ss_h.type='events';
ss_h.filter=ones(1,s_nbevent);

for s_i=1:s_nbevent
    a_fill='*******';
    ss_h.values(s_i).value=m_pos(s_i,2);
    ss_h.values(s_i).authorized=1;
    ss_h.values(s_i).description.id=m_pos(s_i,2);
    ss_h.values(s_i).description.sample_latency=m_pos(s_i,1);
end; % for s_i
fff_data_mono.dim(3)=ss_h;
fff_data_bipo.dim(3)=ss_h;



% now creates the data

disp(['reading data mono']);
f_in=fopen(a_eegfile,'r','ieee-be');
for s_e=1:s_nbevent
    % get the data for this latency
    s_adress_of_data=fff_data_mono.dim(3).values(s_e).description.sample_latency; % the latency of the event for the first channel in this file
    fseek(f_in,(s_adress_of_data+v_window_sam(1))*s_nbchannel*2,-1);
    m_data=fread(f_in,[s_nbchannel,v_window_sam(2)-v_window_sam(1)+1],'int16');
    if (s_e == 1)
        m_bigdata_mono = zeros(size(m_data,2),size(m_data,1),s_nbevent);
        m_bigdata_mono(:,:,1) = m_data';
    else
        m_bigdata_mono(:,:,s_e) = m_data';
    end;        
end; % for s_e
fclose(f_in);




disp(['reading data bipo']);
f_in=fopen(a_eegfile,'r','ieee-be');
for s_e=1:s_nbevent
    % get the data for this latency
    s_adress_of_data=fff_data_bipo.dim(3).values(s_e).description.sample_latency; % the latency of the event for the first channel in this file
    fseek(f_in,(s_adress_of_data+v_window_sam(1))*s_nbchannel*2,-1);
    m_data=fread(f_in,[s_nbchannel,v_window_sam(2)-v_window_sam(1)+1],'int16');
    m_databig=m_data(v_bigsister_o,:);
    m_datalit = 0*m_databig; % in case there is a 0 as little sister : it means monopolar
    v_f = find(v_litsister_o);
    m_datalit(v_f,:)=m_data(v_litsister_o(v_f),:); 
    m_data=m_databig-m_datalit;
    if (s_e == 1)
        m_bigdata_bipo = zeros(size(m_data,2),size(m_data,1),s_nbevent);
        m_bigdata_bipo(:,:,1) = m_data';
    else
        m_bigdata_bipo(:,:,s_e) = m_data';
    end;        
end; % for s_e
fclose(f_in);


% values
% We stop here, ....
clear ss_h
ss_h.id='voltage';
ss_h.unit='default';
ss_h.type='linear';
ss_h.logic_min=0;
ss_h.logic_max=0;
ss_h.phys_min=0;
ss_h.phys_max=0;
ss_h.comment='none';
ss_h.where='here';
ss_h.getme=m_bigdata_mono;
ss_h.data_format='float32';

fff_data_mono.value=ss_h;



% values
% We stop here, ....
clear ss_h
ss_h.id='voltage';
ss_h.unit='default';
ss_h.type='linear';
ss_h.logic_min=0;
ss_h.logic_max=0;
ss_h.phys_min=0;
ss_h.phys_max=0;
ss_h.comment='none';
ss_h.where='here';
ss_h.getme=m_bigdata_bipo;
ss_h.data_format='float32';

fff_data_bipo.value=ss_h;
