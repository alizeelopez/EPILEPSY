function [rapport, rap]=check_data(id_fixpos,id_cue,id_resp)

if length(id_fixpos)~=120
    rapport{1,:}=['Number of id fixation is not 120 : ' num2str(length(id_fixpos)) ', we try to correct that'];
    rap(1)=0;
else
    rapport{1,:}=['Number of id fixation is 120 : everything''s ok'];
    rap(1)=1;
end
if length(id_cue)~=120
    rapport{2,:}=['Number of id cue is not 120 : ' num2str(length(id_cue)) ', we try to correct that'];
    rap(2)=0;
else
    rapport{2,:}=['Number of id cue is 120 : everything''s ok'];
    rap(2)=1;
end

if length(id_resp)~=120
    rapport{3,:}=['Number of id response is not 120 : ' num2str(length(id_resp)) ', we try to correct that'];
    rap(3)=0;
else
    rapport{3,:}=['Number of id response is 120 : everything''s ok'];
    rap(3)=1;
end
% display(rapport);