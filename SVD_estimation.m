%SVD∑÷Ω‚À„∑®∫Øîµ
function output=SVD_estimation(input,pilot_inter,pilot_sequence,pilot_num,trms,t_max,snr,cp);
beta=17/9;
[N,NL]=size(input);
Rhh=zeros(N,N);
for k=1:N
    for l=1:N
        Rhh(k,l)=(1-exp((-1)*t_max*((1/trms)+j*2*pi*(k-l)/N)))./(trms*(1-exp((-1)*t_max/trms))*((1/trms)+j*2*pi*(k-l)/N));
    end
end
[U,D]=eig(Rhh);
dlamda=diag(D);
[dlamda_sort,IN]=sort(dlamda);
for k=1:N
    lamda_sort(k)=dlamda_sort(N-k+1);
    IN_new(k)=IN(N-k+1);
end
U_new=zeros(N,N);
for k=1:N
    U_new(:,k)=U(:,IN_new(k));
end
delta=zeros(N,1);
for k=1:N
    if k<=cp
       delta(k)=lamda_sort(k)/(lamda_sort(k)+beta/snr);
    else
        delta(k)=0;
    end
end
D_new=diag(delta);
output=zeros(N,NL-pilot_num);
i=1;
count=0;
while i<=NL
    Hi=input(:,i)./pilot_sequence;
    Hlr=U_new*D_new*(U_new')*Hi;
    count=count+1;
    if  count*pilot_inter<=(NL-pilot_num)
        for p=((count-1)*pilot_inter+1):count*pilot_inter
            output(:,p)=input(:,(i+p-(count-1)*pilot_inter))./Hlr;
        end
    else
        for p=((count-1)*pilot_inter+1):(NL-pilot_num)
            output(:,p)=input(:,(i+p-(count-1)*pilot_inter))./Hlr;
        end
    end
    i=i+pilot_inter+1;
  end