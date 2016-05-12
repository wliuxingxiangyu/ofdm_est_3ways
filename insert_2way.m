%Ofdm系統仿真程序  块狀和梳狀插入成式在不用情況下之比較
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%梳狀插入
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
%定義
pilot_inter=5;%导頻間隔
pilot_symbol_bit=[0 0 0 1];%导頻符號
carrier_count=256;%子載波數
bits_per_symbol_16QAM=4;%采用該調制时一符號占位數
cp_length=16;
SNR_dB=[0 4 8 12 16 20];
ls_err_ber=zeros(1,length(SNR_dB));
symbols_per_carrier=100;
%Jake's Model setting
Datalen=carrier_count*symbols_per_carrier;
N=carrier_count;
fm=200;   %这里設定doppler頻率,0為慢衰落,200為快衰落
fs=1e6;      
startp=0;
No=9;
L=8;
[Hf] =jakes(Datalen,L,fm,fs,No,startp);
H_power=([exp(0) exp(-0.2) exp(-0.4) exp(-0.6) exp(-0.8) exp(-1) exp(-1.2) exp(-1.4)])';
Hf=H_power*ones(1,L).*eye(L)*Hf;
%Jake's Model setting OK...
%每个SNR点上仿真若干次
for i=1:length(SNR_dB) 
    ls_error_bit=0;%誤碼个數,下同
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
[insert_pilot_out,pilot_num,pilot_code]=insert_pilot2(pilot_inter,pilot_symbol_bit,map_out);%采用insert_pilot2函數,為梳狀插入方式
%Jake's Model影响
Hb=zeros(size(insert_pilot_out+pilot_num));
for er=1:symbols_per_carrier
ggb=conj(Hf(:,er)');
H=fft(ggb,carrier_count+pilot_num);
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
multipath_signal=ofdm_cp_out;
multipath_signal=multi_chan(Tx_data,ofdm_cp_out_m,ofdm_cp_out_n);%若测定非多徑信道下则情況则在本行头加上"%"號
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
ofdm_demodulation_out=conj(ofdm_demodulation_out');
ls_zf_detect_sig=ls_estimation(ofdm_demodulation_out,pilot_inter,pilot_code,pilot_num);%采用LS估计算法测得到的接收信号
ls_zf_detect_sig=conj(ls_zf_detect_sig');
%解映射
    ls_receive_bit_sig=de_map(ls_zf_detect_sig);    
    %以下过程统计各种估计算法得到的接收信号中的错误比特数
    ls_err_num=error_count(convert_matrix,ls_receive_bit_sig);
    ls_error_bit=ls_error_bit+ls_err_num;
%计算各种估计算法的误比特率
ls_err_ber(i)=ls_error_bit/total_bit_num;
end
plot(SNR_dB,ls_err_ber,'b-*')
title('快衰落情况下块状和梳状导频插入方式比较(LS准则)');
xlabel('SNR')
ylabel('BER')
grid on
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%块狀插入
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%定義
pilot_inter=5;
pilot_symbol_bit=[0 0 0 1];
carrier_count=256;
bits_per_symbol_16QAM=4;
cp_length=16;
SNR_dB=[0 4 8 12 16 20];
ls_err_ber=zeros(1,length(SNR_dB));
no_est_error_ber=zeros(1,length(SNR_dB));
symbols_per_carrier=100;
%Jake's Model setting
Datalen=carrier_count*symbols_per_carrier;
N=carrier_count;
fm=200;   %这里設定doppler頻率,0為慢衰落,200為快衰落
fs=1e6;      
startp=0;
No=9;
L=8;
[Hf] =jakes(Datalen,L,fm,fs,No,startp);
H_power=([exp(0) exp(-0.2) exp(-0.4) exp(-0.6) exp(-0.8) exp(-1) exp(-1.2) exp(-1.4)])';
Hf=H_power*ones(1,L).*eye(L)*Hf;
%Jake's Model setting ok...
%每个SNR点上仿真若干次
for i=1:length(SNR_dB) 
    ls_error_bit=0;
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
[insert_pilot_out,pilot_num,pilot_code]=insert_pilot(pilot_inter,pilot_symbol_bit,map_out);%采用insert_pilot函數,為块狀插入方式
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
%送出(考慮和多徑合在一起)
ofdm_modulation=reshape(ofdm_cp_out,1,ofdm_cp_out_m*ofdm_cp_out_n);
Tx_data=ofdm_modulation;
%多径信道
multipath_signal=ofdm_cp_out;
multipath_signal=multi_chan(Tx_data,ofdm_cp_out_m,ofdm_cp_out_n);%若测定非多徑信道下则情況则在本行头加上"%"號
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
no_detect_sig=de_p(ofdm_demodulation_out,pilot_inter,pilot_code,pilot_num);%不作估計時测得到的接收信号,直接把导頻删去
%解映射
ls_receive_bit_sig=de_map(ls_zf_detect_sig);%16QAM解映射
no_est_receive_bit_sig=de_map(no_detect_sig);%16QAM解映射
    %以下过程统计各种估计算法得到的接收信号中的错误比特数
    ls_err_num=error_count(convert_matrix,ls_receive_bit_sig);
    no_est_err_num=error_count(convert_matrix,no_est_receive_bit_sig);
    
    ls_error_bit=ls_error_bit+ls_err_num;
    no_est_error_bit=no_est_error_bit+no_est_err_num;
%计算各种估计算法的误比特率
ls_err_ber(i)=ls_error_bit/total_bit_num;
no_est_error_ber(i)=no_est_error_bit/total_bit_num;
end
plot(SNR_dB,ls_err_ber,'r-o',SNR_dB,no_est_error_ber,'g--')
legend('梳状','块状','不作估計')
grid on
hold on