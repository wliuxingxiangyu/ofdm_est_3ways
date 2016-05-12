%原來的块狀插入导頻函數,不完善,已棄用
function [output_with_pilot,pilot_num,pilot_code]=insert_pilot(pilot_inter,pilot_symbol_bit,output);
[c,s]=size(output);
pilot_num= floor(s/pilot_inter)+1;
pilot_map=qam16(pilot_symbol_bit);
pilot_code= pilot_map;
%output=zeros(c,s+pilot_num);
turn=reshape(conj(output'),pilot_inter,c*s/pilot_inter);
pilot_sequence=zeros(1, c*s/pilot_inter);
for a=1:c*s/pilot_inter
    pilot_sequence(1,a)= pilot_map;
end
with_pilot1=[pilot_sequence;turn];
[c2,s2]=size(with_pilot1);
with_pilot2=reshape(with_pilot1,c2*s2/c,c);
output_with_pilot1= conj(with_pilot2');
pilot_sequence_last=zeros(c,1);
for a=1:c
    pilot_sequence_last(a,1)= pilot_map;
end
output_with_pilot=zeros(c,s+pilot_num);
output_with_pilot=[output_with_pilot1,pilot_sequence_last];
end