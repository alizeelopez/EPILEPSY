function v_sigout=bpfilter(v_sigin,s_lofreq,s_hifreq,s_samfreq)

% v_sigout=bpfilter(v_sigin,s_lofreq,s_hifreq,s_samfreq)
%
% this function band-pass filter signal v_sig between s_lowfreq and s_hifreq without changing its phase
%
% INPUTS
%
% v_sigin     : signal to be filtered. 1-d vector
% s_lofreq : lower frequency, in Hz
% s_hifreq  : higher frequency, in Hz
% s_samfreq : sampling frequency of v_sig, in Hz
%
% OUTPUTS
%
% v_sigout : filtered signal , same size as v_sigin
%
% J.Ph.Lachaux + MARIO - 6-Dec-01. please report any bug to jean-philippe.lachaux@chups.jussieu.fr

s_mario=1; % set it to 1 if you want to use Mario's version

if (s_mario)
   
   
   
   
   s_test=0;
   if (size(v_sigin,1)>size(v_sigin,2))
      v_sigin=v_sigin';
   	s_test=1;   
   end;
   
	v_sigout = bandpassFilter(v_sigin,s_samfreq,s_lofreq,s_hifreq);   
   
   % this is to make sure that v_sigin and v_sigout have the same dimension
   if (s_test)
      % if the samples used to be the first dimension of v_sigin in the user input
      if (size(v_sigout,1)<=size(v_sigout,2))
   	   v_sigout=v_sigout';
   	end;
   else
      if (size(v_sigout,1)>size(v_sigout,2))
   	   v_sigout=v_sigout';
   	end;
   end;
   
   
   
   
   
else   % if s_mario

if (size(v_sigin,1)>size(v_sigin,2))
   v_sigin=v_sigin';
end;

s_wn_lo=s_lofreq/(0.5*s_samfreq);
v_b_lo=n_fir1(350,s_wn_lo);
v_siglo=conv2(v_sigin,v_b_lo,'same');

v_siglo=v_sigin-v_siglo;

s_wn_hi=s_hifreq/(0.5*s_samfreq);
v_b_hi=n_fir1(350,s_wn_hi);
v_sighi=conv2(v_sigin,v_b_hi,'same');

v_sighi=v_sigin-v_sighi;

v_sigout=v_siglo-v_sighi;

end; % if s_mario

%%%% Mario's Function
