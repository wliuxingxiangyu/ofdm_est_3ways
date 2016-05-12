%块狀插入导頻函數,測試用,亦已棄用
function [output_with_pilot,pilot_num,pilot_code]=insert_pilot(pilot_inter,pilot_symbol_bit,output);
[c,s]=size(output);
pilot_num= c*(s/ pilot_inter);
turn=reshape(output',pilot_inter,pilot_num);
pilot_map=qam16(pilot_symbol_bit);
pilot_code= pilot_map;
pilot_sequence=zeros(1, pilot_num);
for a=1:pilot_num
    pilot_sequence(1,a)= pilot_map;
end
with_pilot1=[pilot_sequence;turn];
[c2,s2]=size(with_pilot1);
with_pilot2=reshape(with_pilot1,c2* pilot_num,c);
output_with_pilot= with_pilot2';