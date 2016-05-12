% ·†Ó≤Â»Î(æQ–‘É»≤Â)∫Øîµ
function [output_with_pilot,pilot_num,pilot_code]=insert_pilot2(pilot_inter,pilot_symbol_bit,output);
output=conj(output');
[c,s]=size(output);
pilot_num= floor(s/pilot_inter)+1;
pilot_map=qam16(pilot_symbol_bit);
pilot_code= pilot_map;
new=zeros(c,s+pilot_num);
for a=1:c
    pilot_sequence(1,a)= pilot_map;
end
a=1;
gg=0;
hh=pilot_inter*floor(s/pilot_inter);
while a<hh
    new(:,a)=pilot_sequence;
    new(:,(a+1):(a+pilot_inter))=output(:,(a-gg):(a+4-gg));
    gg=gg+1;
a=a+pilot_inter+1;
end
    (s+pilot_num-1-(s-hh-1));
    (s+pilot_num-1);
    hh+1;
    s;
new(:,(s+pilot_num-1-(s-hh-1)):(s+pilot_num-1))=output(:,(hh+1:s));
new(:,s+pilot_num)=pilot_sequence;
output_with_pilot=new;
output_with_pilot=conj(output_with_pilot');
end