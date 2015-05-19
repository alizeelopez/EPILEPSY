function eeg2eeg_chg_pos(a_oldpos,a_newpos,s_downsamp)

m_pos=load(a_oldpos);

v_sam = m_pos(:,1);
v_sam = round(v_sam/s_downsamp);

m_pos(:,1) = v_sam;

f=fopen(a_newpos,'w');

for s_i=1:size(m_pos,1)
    s_id=m_pos(s_i,2);
    s_latency=m_pos(s_i,1);
    fprintf(f,'%d\t%d\t%d\n',s_latency,s_id,0);
end;

fclose(f);
