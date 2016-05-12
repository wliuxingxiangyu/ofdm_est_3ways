%Ofdm系統仿真程序  比較不同估計方法的誤碼性能
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%測量不作估計時和采用LS估计算法時的誤碼性能
clear all;
close all;
%定義
pilot_inter=5;%导頻間隔
pilot_symbol_bit=[0 0 0 1];%导頻符號
carrier_count=128;%子載波數
bits_per_symbol_16QAM=4;%采用該調制时一符號占位數
cp_length=16;%cp長度
SNR_dB=[0 2 4 6 8 10 12 14 16];%不同的SNR
ls_err_ber=zeros(1,length(SNR_dB));
no_est_error_ber=zeros(1,length(SNR_dB));
for i=1:length(SNR_dB) %每个SNR点上仿真若干次
    ls_error_bit=0;
    no_est_error_bit=0;
    total_bit_num=0;
symbols_per_carrier=50;
baseband_out_num=carrier_count*(symbols_per_carrier+pilot_inter)*bits_per_symbol_16QAM;
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
%送出
ofdm_modulation=reshape(ofdm_cp_out,1,ofdm_cp_out_m*ofdm_cp_out_n);
Tx_data=ofdm_modulation;
%静态信道（多径）
multipath_signal=multi_chan(Tx_data,ofdm_cp_out_m,ofdm_cp_out_n);
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
    sgma=sqrt(spow1/(2*snr));%sgma计算方式，与当前SNR和信号平均能量有关系
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
plot(SNR_dB,no_est_error_ber,'b-*',SNR_dB,ls_err_ber,'r-o')
title('多径情况下插入导频,训练符号及不作估计方式比较')
xlabel('SNR')
ylabel('BER')
grid on
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ofdm系統仿真程序  訓練序列 有估計 16qam  有多徑 
%定義
%pilot_inter=5;
%pilot_symbol_bit=[0 0 0 1];
carrier_count=128;
bits_per_symbol_16QAM=4;
cp_length=16;
SNR_dB=[0 2 4 6 8 10 12 14 16];
training_err_ber=zeros(1,length(SNR_dB));
for i=1:length(SNR_dB) %每个SNR点上仿真若干次
    training_err_bit=0;
    total_bit_num=0;
symbols_per_carrier=50;
baseband_out_num=carrier_count*(symbols_per_carrier)*bits_per_symbol_16QAM;
baseband_out=round(rand(1, baseband_out_num));   
%二进制数据映射为16进制数
convert_matrix=reshape(baseband_out, bits_per_symbol_16QAM*carrier_count,length(baseband_out)/ (bits_per_symbol_16QAM*carrier_count));
total_bit_num=baseband_out_num;
%影射
map_out=map_16_QAM(convert_matrix);
%不插入导頻
%[insert_pilot_out,pilot_num,pilot_code]=insert_pilot(pilot_inter,pilot_symbol_bit,map_out);
%插入訓練序列
training_symbols=[1+j;1-j;-1+j;-1-j;3+j;1+3*j;3-j;1-3*j;-1+3*j;-3+j;-3-j;-1-3*j;3+3*j;-3+3*j;-3-3*j;3-3*j];
training_symbols=[training_symbols;training_symbols];
training_symbols=[training_symbols;training_symbols];
training_symbols=[training_symbols;training_symbols];
training_symbols=[training_symbols,training_symbols];
training_symbols=[training_symbols,training_symbols];   %造出一個128*4的訓練符號矩陣
insert_training_out=[training_symbols,map_out];    %合井到要傳的數据上
train_length=size(training_symbols,2);
total_length=size(insert_training_out,2);
%FFT變換
ofdm_mode_out=ifft(insert_training_out,128);
%插入循环前缀
ofdm_cp_out=insert_cp(ofdm_mode_out,cp_length);
[ofdm_cp_out_m,ofdm_cp_out_n]=size(ofdm_cp_out);
%送出
ofdm_modulation=reshape(ofdm_cp_out,1,ofdm_cp_out_m*ofdm_cp_out_n);
Tx_data=ofdm_modulation;
%静态信道（多径）
multipath_signal=multi_chan(Tx_data,ofdm_cp_out_m,ofdm_cp_out_n);
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
    sgma=sqrt(spow1/(2*snr));%sgma计算方式，与当前SNR和信号平均能量有关系
    receive_ofdm_symbol=add_noise(sgma,passchan_ofdm_symbol);%加入随机高斯白噪声,receive_ofdm_symbol为最终接收机收到的ofdm符号块
%去循環前綴
cut_cp_symbol=del_cp(receive_ofdm_symbol,cp_length);
%FFT
ofdm_demodulation_out=fft(cut_cp_symbol,128);%作128点FFT运算，完成ofdm解调

Rx_carriers = ofdm_demodulation_out(:,(train_length+1):total_length );
Rx_training_symbols = ofdm_demodulation_out(:,(1: train_length)  );
Rx_training_symbols1=Rx_training_symbols;
%信道估计
Rx_training_symbols = Rx_training_symbols./ training_symbols;
Rx_training_symbols_deno = Rx_training_symbols.^2;
Rx_training_symbols_deno = Rx_training_symbols_deno(:,1)+Rx_training_symbols_deno(:,2)+Rx_training_symbols_deno(:,3)+Rx_training_symbols_deno(:,4) ;
Rx_training_symbols_nume = Rx_training_symbols(:,1) +Rx_training_symbols(:,2) + Rx_training_symbols(:,3) +Rx_training_symbols(:,4) ;
Rx_training_symbols_nume = conj(Rx_training_symbols_nume) ;
% 取4个向量的导频符号是为了进行平均优化
% 都是针对 “行向量”即单个的OFDM符号 进行操作
% 原理：寻求1/H，对FFT之后的数据进行频域补偿
% 1/H = conj(H)/H^2 because H^2 = H * conj(H)
Rx_training_symbols = Rx_training_symbols_nume./Rx_training_symbols_deno;
Rx_training_symbols_2 =[Rx_training_symbols,Rx_training_symbols];
Rx_training_symbols_4 =[Rx_training_symbols_2,Rx_training_symbols_2];
Rx_training_symbols_8 =[Rx_training_symbols_4,Rx_training_symbols_4];
Rx_training_symbols_16 =[Rx_training_symbols_8,Rx_training_symbols_8];
Rx_training_symbols_32 =[Rx_training_symbols_16,Rx_training_symbols_16];
Rx_training_symbols_48 =[Rx_training_symbols_32,Rx_training_symbols_16];
Rx_training_symbols_50 =[Rx_training_symbols_48,Rx_training_symbols_2];
Rx_carriers = Rx_training_symbols_50.*Rx_carriers;
%解映射
    receive_bit_sig=de_map(Rx_carriers);
    %以下过程统计各种估计算法得到的接收信号中的错误比特数
    err_num=error_count(convert_matrix,receive_bit_sig);
    trainning_error_bit=training_err_bit+err_num;
%计算各种估计算法的误比特率
trainning_err_ber(i)=trainning_error_bit/total_bit_num;
end
plot(SNR_dB,trainning_err_ber,'g-+')
legend('不作估计','插入导频(ls准测)估计','利用训练符号的估计')
grid on
hold on
