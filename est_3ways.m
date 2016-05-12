%Ofdm系統仿真程序, LS,LMMSE,SVD分解三种估计算法的比较
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%定義
clear all;
close all;
pilot_inter=5;%导頻間隔
pilot_symbol_bit=[0 0 0 1];%导頻符號
carrier_count=256;%子載波數
bits_per_symbol_16QAM=4;%采用該調制时一符號占位數
cp_length=16;%cp長度
SNR_dB=[0 4 8 12 16 20];%不同的SNR
ls_err_ber=zeros(1,length(SNR_dB));
lmmse_err_ber=zeros(1,length(SNR_dB));
SVD_err_ber=zeros(1,length(SNR_dB));%不同算法估計后的BER
no_est_error_ber=zeros(1,length(SNR_dB));%不作估計
symbols_per_carrier=100;
%时間參數,lmmse和svd算法時用
trms=4e-6;
t_interval=1e-6;%采样间隔为1us
trms_1=trms/t_interval;
t_max=16e-6/t_interval;
%Jake's Model　setting
Datalen=carrier_count*symbols_per_carrier;
N=carrier_count;
fm=20;   %作多普勒頻率用
fs=1e6;      
startp=0;
No=9;
L=8;
[Hf] =jakes(Datalen,L,fm,fs,No,startp);
H_power=([exp(0) exp(-0.2) exp(-0.4) exp(-0.6) exp(-0.8) exp(-1) exp(-1.2) exp(-1.4)])';
Hf=H_power*ones(1,L).*eye(L)*Hf;
%Jake's Model　setting OK...
for i=1:length(SNR_dB) %每个SNR点上仿真若干次
    ls_error_bit=0;%誤碼个數,下同
    lmmse_error_bit=0;
    SVD_error_bit=0;
    no_est_error_bit=0;
    total_bit_num=0;
%信號初始化
baseband_out_num=carrier_count*(symbols_per_carrier)*bits_per_symbol_16QAM;
baseband_out=round(rand(1, baseband_out_num));   
%二进制数据映射为16进制数
convert_matrix=reshape(baseband_out, bits_per_symbol_16QAM*carrier_count,length(baseband_out)/ (bits_per_symbol_16QAM*carrier_count));
total_bit_num=baseband_out_num;
%影射
map_out=map_16_QAM(convert_matrix);
%插入导頻
[insert_pilot_out,pilot_num,pilot_code]=insert_pilot(pilot_inter,pilot_symbol_bit,map_out);
%Jake's Model影响
Hb=zeros(size(insert_pilot_out));
for er=1:symbols_per_carrier+pilot_num
    ggc=conj(Hf(:,er)');
    H=fft(ggc,carrier_count);
    Hb(:,er)=H;
end
insert_pilot_out=insert_pilot_out.*Hb;
%FFT變換
ofdm_mode_out=ifft(insert_pilot_out,128);
%插入循环前缀
ofdm_cp_out=insert_cp(ofdm_mode_out,cp_length);
[ofdm_cp_out_m,ofdm_cp_out_n]=size(ofdm_cp_out);
%送出
ofdm_modulation=reshape(ofdm_cp_out,1,ofdm_cp_out_m*ofdm_cp_out_n);
Tx_data=ofdm_modulation;
%多径信道
multipath_signal=multi_chan(Tx_data,ofdm_cp_out_m,ofdm_cp_out_n);
%multipath_signal=ofdm_cp_out;
passchan_ofdm_symbol=reshape(multipath_signal,ofdm_cp_out_m,ofdm_cp_out_n);
%加噪聲
snr=10^(SNR_dB(i)/10);
    [nnl1,mml1]=size(passchan_ofdm_symbol);
    spow=0;
    for k=1:nnl1
      for b=1:mml1
        spow=spow+real(passchan_ofdm_symbol(k,b))^2+imag(passchan_ofdm_symbol(k,b))^2;
      end
    end
    spow1=spow/(nnl1*mml1);        
    sgma=sqrt(spow1/(2*snr));%sgma如何计算，与当前SNR和信号平均能量有关系
    receive_ofdm_symbol=add_noise(sgma,passchan_ofdm_symbol);%加入随机高斯白噪声,receive_ofdm_symbol为最终接收机收到的ofdm符号块
%去循環前綴
cut_cp_symbol=del_cp(receive_ofdm_symbol,cp_length);
%FFT
ofdm_demodulation_out=fft(cut_cp_symbol,128);%作128点FFT运算，完成ofdm解调
%估計
ls_zf_detect_sig=ls_estimation(ofdm_demodulation_out,pilot_inter,pilot_code,pilot_num);%采用LS估计算法测得到的接收信号
lmmse_zf_detect_sig=lmmse_estimation(ofdm_demodulation_out,pilot_inter,pilot_code,pilot_num,trms_1,t_max,snr);%采用LMMSE估计算法测得到的接收信号
SVD_sig=SVD_estimation(ofdm_demodulation_out,pilot_inter,pilot_code,pilot_num,trms_1,t_max,snr,cp_length);%采用SVD分解估计算法测得到的接收信号
no_detect_sig=de_p(ofdm_demodulation_out,pilot_inter,pilot_code,pilot_num);%不采用任何估计算法测得到的接收信号
%解映射
    ls_receive_bit_sig=de_map(ls_zf_detect_sig);
    lmmse_receive_bit_sig=de_map(lmmse_zf_detect_sig);
    SVD_receive_bit_sig=de_map(SVD_sig);
    no_est_receive_bit_sig=de_map(no_detect_sig);
%以下过程统计各种估计算法得到的接收信号中的错误比特数
    ls_err_num=error_count(convert_matrix,ls_receive_bit_sig);
    lmmse_err_num=error_count(convert_matrix,lmmse_receive_bit_sig);
    SVD_err_num=error_count(convert_matrix,SVD_receive_bit_sig);
    no_est_err_num=error_count(convert_matrix,no_est_receive_bit_sig);
    
    ls_error_bit=ls_error_bit+ls_err_num;
    lmmse_error_bit=lmmse_error_bit+lmmse_err_num;
    SVD_error_bit=SVD_error_bit+SVD_err_num;
    no_est_error_bit=no_est_error_bit+no_est_err_num;
%end
%计算各种估计算法的误比特率
ls_err_ber(i)=ls_error_bit/total_bit_num;
lmmse_err_ber(i)=lmmse_error_bit/total_bit_num;
SVD_err_ber(i)=SVD_error_bit/total_bit_num;
no_est_error_ber(i)=no_est_error_bit/total_bit_num;
end
plot(SNR_dB,ls_err_ber,'r-o',SNR_dB,lmmse_err_ber,'b-*',SNR_dB,SVD_err_ber,'g-+',SNR_dB,no_est_error_ber,'y--')
legend('LS准则','LMMSE准则','SVD分解','不作估計')
title('三种估计算法的比较');
xlabel('SNR')
ylabel('BER')
grid on
hold on