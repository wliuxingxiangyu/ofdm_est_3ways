%¶à½Ë¥ÂäÕ{ÓÃº¯”µ
function output=multi_chan(input,m,n)
delay=[0,2e-6,4e-6,8e-6,12e-6];
trm=4e-6;
copy=zeros(size(input));
for k=1:length(delay)
	for i=1+k:length(input)
    		copy(i)=(exp(-delay(k)/trm))*input(i-k);
	end
input=input+copy;
end
output=reshape(input,m,n);