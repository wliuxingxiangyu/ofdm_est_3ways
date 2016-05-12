%Ofdmϵ�y�������, LS,LMMSE,SVD�ֽ����ֹ����㷨�ıȽ�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���x
clear all;
close all;
pilot_inter=5;%���l�g��
pilot_symbol_bit=[0 0 0 1];%���l��̖
carrier_count=256;%���d����
bits_per_symbol_16QAM=4;%����ԓ�{��ʱһ��̖ռλ��
cp_length=16;%cp�L��
SNR_dB=[0 4 8 12 16 20];%��ͬ��SNR
ls_err_ber=zeros(1,length(SNR_dB));
lmmse_err_ber=zeros(1,length(SNR_dB));
SVD_err_ber=zeros(1,length(SNR_dB));%��ͬ�㷨��Ӌ���BER
no_est_error_ber=zeros(1,length(SNR_dB));%������Ӌ
symbols_per_carrier=100;
%ʱ�g����,lmmse��svd�㷨�r��
trms=4e-6;
t_interval=1e-6;%�������Ϊ1us
trms_1=trms/t_interval;
t_max=16e-6/t_interval;
%Jake's Model��setting
Datalen=carrier_count*symbols_per_carrier;
N=carrier_count;
fm=20;   %���������l����
fs=1e6;      
startp=0;
No=9;
L=8;
[Hf] =jakes(Datalen,L,fm,fs,No,startp);
H_power=([exp(0) exp(-0.2) exp(-0.4) exp(-0.6) exp(-0.8) exp(-1) exp(-1.2) exp(-1.4)])';
Hf=H_power*ones(1,L).*eye(L)*Hf;
%Jake's Model��setting OK...
for i=1:length(SNR_dB) %ÿ��SNR���Ϸ������ɴ�
    ls_error_bit=0;%�`�a����,��ͬ
    lmmse_error_bit=0;
    SVD_error_bit=0;
    no_est_error_bit=0;
    total_bit_num=0;
%��̖��ʼ��
baseband_out_num=carrier_count*(symbols_per_carrier)*bits_per_symbol_16QAM;
baseband_out=round(rand(1, baseband_out_num));   
%����������ӳ��Ϊ16������
convert_matrix=reshape(baseband_out, bits_per_symbol_16QAM*carrier_count,length(baseband_out)/ (bits_per_symbol_16QAM*carrier_count));
total_bit_num=baseband_out_num;
%Ӱ��
map_out=map_16_QAM(convert_matrix);
%���뵼�l
[insert_pilot_out,pilot_num,pilot_code]=insert_pilot(pilot_inter,pilot_symbol_bit,map_out);
%Jake's ModelӰ��
Hb=zeros(size(insert_pilot_out));
for er=1:symbols_per_carrier+pilot_num
    ggc=conj(Hf(:,er)');
    H=fft(ggc,carrier_count);
    Hb(:,er)=H;
end
insert_pilot_out=insert_pilot_out.*Hb;
%FFT׃�Q
ofdm_mode_out=ifft(insert_pilot_out,128);
%����ѭ��ǰ׺
ofdm_cp_out=insert_cp(ofdm_mode_out,cp_length);
[ofdm_cp_out_m,ofdm_cp_out_n]=size(ofdm_cp_out);
%�ͳ�
ofdm_modulation=reshape(ofdm_cp_out,1,ofdm_cp_out_m*ofdm_cp_out_n);
Tx_data=ofdm_modulation;
%�ྶ�ŵ�
multipath_signal=multi_chan(Tx_data,ofdm_cp_out_m,ofdm_cp_out_n);
%multipath_signal=ofdm_cp_out;
passchan_ofdm_symbol=reshape(multipath_signal,ofdm_cp_out_m,ofdm_cp_out_n);
%����
snr=10^(SNR_dB(i)/10);
    [nnl1,mml1]=size(passchan_ofdm_symbol);
    spow=0;
    for k=1:nnl1
      for b=1:mml1
        spow=spow+real(passchan_ofdm_symbol(k,b))^2+imag(passchan_ofdm_symbol(k,b))^2;
      end
    end
    spow1=spow/(nnl1*mml1);        
    sgma=sqrt(spow1/(2*snr));%sgma��μ��㣬�뵱ǰSNR���ź�ƽ�������й�ϵ
    receive_ofdm_symbol=add_noise(sgma,passchan_ofdm_symbol);%���������˹������,receive_ofdm_symbolΪ���ս��ջ��յ���ofdm���ſ�
%ȥѭ�hǰ�Y
cut_cp_symbol=del_cp(receive_ofdm_symbol,cp_length);
%FFT
ofdm_demodulation_out=fft(cut_cp_symbol,128);%��128��FFT���㣬���ofdm���
%��Ӌ
ls_zf_detect_sig=ls_estimation(ofdm_demodulation_out,pilot_inter,pilot_code,pilot_num);%����LS�����㷨��õ��Ľ����ź�
lmmse_zf_detect_sig=lmmse_estimation(ofdm_demodulation_out,pilot_inter,pilot_code,pilot_num,trms_1,t_max,snr);%����LMMSE�����㷨��õ��Ľ����ź�
SVD_sig=SVD_estimation(ofdm_demodulation_out,pilot_inter,pilot_code,pilot_num,trms_1,t_max,snr,cp_length);%����SVD�ֽ�����㷨��õ��Ľ����ź�
no_detect_sig=de_p(ofdm_demodulation_out,pilot_inter,pilot_code,pilot_num);%�������κι����㷨��õ��Ľ����ź�
%��ӳ��
    ls_receive_bit_sig=de_map(ls_zf_detect_sig);
    lmmse_receive_bit_sig=de_map(lmmse_zf_detect_sig);
    SVD_receive_bit_sig=de_map(SVD_sig);
    no_est_receive_bit_sig=de_map(no_detect_sig);
%���¹���ͳ�Ƹ��ֹ����㷨�õ��Ľ����ź��еĴ��������
    ls_err_num=error_count(convert_matrix,ls_receive_bit_sig);
    lmmse_err_num=error_count(convert_matrix,lmmse_receive_bit_sig);
    SVD_err_num=error_count(convert_matrix,SVD_receive_bit_sig);
    no_est_err_num=error_count(convert_matrix,no_est_receive_bit_sig);
    
    ls_error_bit=ls_error_bit+ls_err_num;
    lmmse_error_bit=lmmse_error_bit+lmmse_err_num;
    SVD_error_bit=SVD_error_bit+SVD_err_num;
    no_est_error_bit=no_est_error_bit+no_est_err_num;
%end
%������ֹ����㷨���������
ls_err_ber(i)=ls_error_bit/total_bit_num;
lmmse_err_ber(i)=lmmse_error_bit/total_bit_num;
SVD_err_ber(i)=SVD_error_bit/total_bit_num;
no_est_error_ber(i)=no_est_error_bit/total_bit_num;
end
plot(SNR_dB,ls_err_ber,'r-o',SNR_dB,lmmse_err_ber,'b-*',SNR_dB,SVD_err_ber,'g-+',SNR_dB,no_est_error_ber,'y--')
legend('LS׼��','LMMSE׼��','SVD�ֽ�','������Ӌ')
title('���ֹ����㷨�ıȽ�');
xlabel('SNR')
ylabel('BER')
grid on
hold on