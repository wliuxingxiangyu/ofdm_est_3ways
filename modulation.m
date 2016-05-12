%Ofdm系統仿真程序   (非)多径情况下BPSK,4QAM和16QAM调制方式的比较
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  BPSK调制方式下
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
%定義
pilot_inter=5;%导頻間隔
pilot_symbol_bit=1;%导頻符號
carrier_count=256;%子載波數
%bits_per_symbol_16QAM=4;%采用該調制时一符號占位數(BPSK不用)
cp_length=16;%cp長度
SNR_dB=[0 2 4 6 8 10 12 14 16];
ls_err_ber=zeros(1,length(SNR_dB));
for i=1:length(SNR_dB) %每个SNR点上仿真若干次
    ls_error_bit=0;
    total_bit_num=0;
symbols_per_carrier=50*4;
loop=5;
%for l=1:loop
baseband_out_num=carrier_count*symbols_per_carrier;
baseband_out=round(rand(1, baseband_out_num));   
for k=1:baseband_out_num
    if baseband_out(k)>=0.5
        baseband_out(k)=1;
    else
        baseband_out(k)=-1;
    end
end
%二进制数据映射为2进制数
convert_matrix=reshape(baseband_out, carrier_count,symbols_per_carrier);
total_bit_num=baseband_out_num;
%影射
%map_out=map_16_QAM(convert_matrix);
map_out=convert_matrix;
%插入导頻
[insert_pilot_out,pilot_num,pilot_code]=insert_pilot(pilot_inter,pilot_symbol_bit,map_out);
%FFT變換
ofdm_mode_out=ifft(insert_pilot_out,128);
%插入循环前缀
ofdm_cp_out=insert_cp(ofdm_mode_out,cp_length);
[ofdm_cp_out_m,ofdm_cp_out_n]=size(ofdm_cp_out);
%送出(考慮和多徑合在一起)
ofdm_modulation=reshape(ofdm_cp_out,1,ofdm_cp_out_m*ofdm_cp_out_n);
Tx_data=ofdm_modulation;
%静态信道（多径）
multipath_signal=ofdm_cp_out;
%multipath_signal=multi_chan(Tx_data,ofdm_cp_out_m,ofdm_cp_out_n);%若测定多徑信道下则情況则解除本行头的"%"號
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
%********************** 以下就是对接收ofdm信号进行信道估计和信号检测的过程************************
ls_zf_detect_sig=ls_estimation(ofdm_demodulation_out,pilot_inter,pilot_code,pilot_num);%采用LS估计算法测得到的接收信号
%解映射
    ls_receive_bit_sig=de_map1(ls_zf_detect_sig);%BPSK解映射
    %以下过程统计各种估计算法得到的接收信号中的错误比特数
    ls_err_num=error_count(convert_matrix,ls_receive_bit_sig);    
    ls_error_bit=ls_error_bit+ls_err_num;
%end
%计算各种估计算法的误比特率
ls_err_ber(i)=ls_error_bit/total_bit_num;
end
plot(SNR_dB,ls_err_ber,'b-*')
title('非多径情况下BPSK,4QAM和16QAM调制方式的比较');
legend('BPSK')
xlabel('SNR')
ylabel('BER')
grid on
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  16QAM调制方式下
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%定義
pilot_inter=5;
pilot_symbol_bit=[0 0 0 1];
carrier_count=256;
bits_per_symbol_16QAM=4;%%采用該調制时一符號占位數(4)
cp_length=16;
SNR_dB=[0 2 4 6 8 10 12 14 16];
ls_err_ber=zeros(1,length(SNR_dB));
lmmse_err_ber=zeros(1,length(SNR_dB));
lr_lmmse_err_ber=zeros(1,length(SNR_dB));
for i=1:length(SNR_dB) %每个SNR点上仿真若干次
    ls_error_bit=0;
    lmmse_error_bit=0;
    lr_lmmse_error_bit=0;
    total_bit_num=0;
symbols_per_carrier=50;
loop=5;
%for l=1:loop
baseband_out_num=carrier_count*(symbols_per_carrier)*bits_per_symbol_16QAM;
baseband_out=round(rand(1, baseband_out_num));   
%二进制数据映射为16进制数
convert_matrix=reshape(baseband_out, bits_per_symbol_16QAM*carrier_count,length(baseband_out)/ (bits_per_symbol_16QAM*carrier_count));
total_bit_num=baseband_out_num;
%影射
map_out=map_16_QAM(convert_matrix);
%插入导頻
[insert_pilot_out,pilot_num,pilot_code]=insert_pilot(pilot_inter,pilot_symbol_bit,map_out);
%FFT變換
ofdm_mode_out=ifft(insert_pilot_out,128);
%插入循环前缀
ofdm_cp_out=insert_cp(ofdm_mode_out,cp_length);
[ofdm_cp_out_m,ofdm_cp_out_n]=size(ofdm_cp_out);
%送出(考慮和多徑合在一起)
ofdm_modulation=reshape(ofdm_cp_out,1,ofdm_cp_out_m*ofdm_cp_out_n);
Tx_data=ofdm_modulation;
%静态信道（多径）
multipath_signal=ofdm_cp_out;
%multipath_signal=multi_chan(Tx_data,ofdm_cp_out_m,ofdm_cp_out_n);%若测定多徑信道下则情況则解除本行头的"%"號
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
%解映射
    ls_receive_bit_sig=de_map(ls_zf_detect_sig);%16QAM解映射
    %以下过程统计各种估计算法得到的接收信号中的错误比特数
    ls_err_num=error_count(convert_matrix,ls_receive_bit_sig);
    ls_error_bit=ls_error_bit+ls_err_num;
%计算各种估计算法的误比特率
ls_err_ber(i)=ls_error_bit/total_bit_num;
end
plot(SNR_dB,ls_err_ber,'r-o')
title('非多径情况下BPSK,4QAM和16QAM调制方式的比较')
legend('插入导頻(ls准測)')
xlabel('SNR')
ylabel('BER')
grid on
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  4QAM调制方式下
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%定義
pilot_inter=5;
pilot_symbol_bit=[0 0 0 1];
carrier_count=256;
bits_per_symbol_4QAM=2;%采用該調制时一符號占位數(2)
cp_length=16;
SNR_dB=[0 2 4 6 8 10 12 14 16];
ls_err_ber=zeros(1,length(SNR_dB));
lmmse_err_ber=zeros(1,length(SNR_dB));
lr_lmmse_err_ber=zeros(1,length(SNR_dB));
for i=1:length(SNR_dB) %每个SNR点上仿真若干次
    ls_error_bit=0;
    total_bit_num=0;
symbols_per_carrier=50;
loop=5;
%for l=1:loop
baseband_out_num=carrier_count*(symbols_per_carrier+pilot_inter)*bits_per_symbol_4QAM;
baseband_out=round(rand(1, baseband_out_num));   
%二进制数据映射为16进制数
convert_matrix=reshape(baseband_out, bits_per_symbol_4QAM*carrier_count,length(baseband_out)/ (bits_per_symbol_4QAM*carrier_count));
total_bit_num=baseband_out_num;
%影射
map_out=map_4_QAM(convert_matrix);
%插入导頻
[insert_pilot_out,pilot_num,pilot_code]=insert_pilot(pilot_inter,pilot_symbol_bit,map_out);
%FFT變換
ofdm_mode_out=ifft(insert_pilot_out,128);
%插入循环前缀
ofdm_cp_out=insert_cp(ofdm_mode_out,cp_length);
[ofdm_cp_out_m,ofdm_cp_out_n]=size(ofdm_cp_out);
%送出(考慮和多徑合在一起)
ofdm_modulation=reshape(ofdm_cp_out,1,ofdm_cp_out_m*ofdm_cp_out_n);
Tx_data=ofdm_modulation;
%静态信道（多径）
multipath_signal=ofdm_cp_out;
%multipath_signal=multi_chan(Tx_data,ofdm_cp_out_m,ofdm_cp_out_n);%若测定多徑信道下
%则情況则解除本行头的"%"號
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
%解映射
    ls_receive_bit_sig=de_map4(ls_zf_detect_sig);%16QAM解映射
    %以下过程统计各种估计算法得到的接收信号中的错误比特数
    ls_err_num=error_count(convert_matrix,ls_receive_bit_sig);    
    ls_error_bit=ls_error_bit+ls_err_num;
%计算各种估计算法的误比特率
ls_err_ber(i)=ls_error_bit/total_bit_num;
lmmse_err_ber(i)=lmmse_error_bit/total_bit_num;
lr_lmmse_err_ber(i)=lr_lmmse_error_bit/total_bit_num;
end
plot(SNR_dB,ls_err_ber,'g-+')
title('非多径情况下BPSK,4QAM和16QAM调制方式的比较')
legend('BPSK','16QAM','4QAM')
hold on
