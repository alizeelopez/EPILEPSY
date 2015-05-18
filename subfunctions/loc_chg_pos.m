function loc_chg_pos(a_oldposfile,a_newposfile)

% this function changes the codes of the pos file in the clock protocol

m_pos=load(a_oldposfile);

v_f = find(m_pos(:,2)==99);
if (~isempty(v_f))
    v_select = (min(v_f)+3:size(m_pos,1));
else
    v_select = (1:size(m_pos,1));
end;

m_pos = m_pos(v_select,:);

f=fopen(a_newposfile,'w');

for s_i=1:size(m_pos,1)
    s_id=m_pos(s_i,2);
    s_latency=m_pos(s_i,1);
    fprintf(f,'%d\t%d\t%d\n',s_latency,s_id,0);
end;

fclose(f);