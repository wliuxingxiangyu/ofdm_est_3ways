%Jake模型的函数  
function  [H] = jakes(Datalen,chanpf,fm,fs,No,startp)
%%%%%%%  调用Jake模型相关的信道参数   %%%%%%%
if (~exist('startp')),
   startp = 0;
end
if (~exist('No')),
   startp = 0;
end

if startp == 0;
   flag = 0;
else
   flag = 1;
end

pathnum = chanpf;
Dopplershift=fm;
wm  = 2*pi*fm;
Ts = 1/fs;           
t = 0:Ts:(Datalen-1)*Ts; 
t = t + startp*Ts;
No = 8;
N = (2*No+1)* 2;
alpha = 0;
belta = pi*(1:No)/(No+1);
w = wm * cos(2*pi*(1:No)/N);   

if flag == 1;
    rand('state',0);  
end
phase_init =rand(No,1)*2*pi;
for k = 1:No
	resin = zeros(No,length(t));
	resqu = zeros(No,length(t));
	for i = 1:No     
        resin(i,:) = 2*cos(belta(i))*cos(w(i)*t + phase_init(k));
        resqu(i,:) = 2*sin(belta(i))*cos(w(i)*t + phase_init(k));
	end
	res2in = sum(resin) + sqrt(2) * cos(alpha) * cos(wm*t);
	res2qu = sum(resqu) + sqrt(2) * sin(alpha) * cos(wm*t);	
	y(k,:) = sqrt(1/(2*No+1)).*(res2in + sqrt(-1)* res2qu) ;
end

for i = 1:No
   y(i,:) = y(i,:)/norm(y(i,:))*sqrt(Datalen);
end

H = y(1:pathnum,1:end);

return;
