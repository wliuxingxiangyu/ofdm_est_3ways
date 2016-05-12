%Ofdmϵ�y�������  ���^��ͬ��Ӌ�������`�a����
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�y��������Ӌ�r�Ͳ���LS�����㷨�r���`�a����
clear all;
close all;
%���x
pilot_inter=5;%���l�g��
pilot_symbol_bit=[0 0 0 1];%���l��̖
carrier_count=128;%���d����
bits_per_symbol_16QAM=4;%����ԓ�{��ʱһ��̖ռλ��
cp_length=16;%cp�L��
SNR_dB=[0 2 4 6 8 10 12 14 16];%��ͬ��SNR
ls_err_ber=zeros(1,length(SNR_dB));
no_est_error_ber=zeros(1,length(SNR_dB));
for i=1:length(SNR_dB) %ÿ��SNR���Ϸ������ɴ�
    ls_error_bit=0;
    no_est_error_bit=0;
    total_bit_num=0;
symbols_per_carrier=50;
baseband_out_num=carrier_count*(symbols_per_carrier+pilot_inter)*bits_per_symbol_16QAM;
baseband_out=round(rand(1, baseband_out_num));   
%����������ӳ��Ϊ16������
convert_matrix=reshape(baseband_out, bits_per_symbol_16QAM*carrier_count,length(baseband_out)/ (bits_per_symbol_16QAM*carrier_count));
total_bit_num=baseband_out_num;
%Ӱ��
map_out=map_16_QAM(convert_matrix);
%���뵼�l
[insert_pilot_out,pilot_num,pilot_code]=insert_pilot(pilot_inter,pilot_symbol_bit,map_out);
%FFT׃�Q
ofdm_mode_out=ifft(insert_pilot_out,128);
%����ѭ��ǰ׺
ofdm_cp_out=insert_cp(ofdm_mode_out,cp_length);
[ofdm_cp_out_m,ofdm_cp_out_n]=size(ofdm_cp_out);
%�ͳ�
ofdm_modulation=reshape(ofdm_cp_out,1,ofdm_cp_out_m*ofdm_cp_out_n);
Tx_data=ofdm_modulation;
%��̬�ŵ����ྶ��
multipath_signal=multi_chan(Tx_data,ofdm_cp_out_m,ofdm_cp_out_n);
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
    sgma=sqrt(spow1/(2*snr));%sgma���㷽ʽ���뵱ǰSNR���ź�ƽ�������й�ϵ
    receive_ofdm_symbol=add_noise(sgma,passchan_ofdm_symbol);%���������˹������,receive_ofdm_symbolΪ���ս��ջ��յ���ofdm���ſ�
%ȥѭ�hǰ�Y
cut_cp_symbol=del_cp(receive_ofdm_symbol,cp_length);
%FFT
ofdm_demodulation_out=fft(cut_cp_symbol,128);%��128��FFT���㣬���ofdm���
%��Ӌ
ls_zf_detect_sig=ls_estimation(ofdm_demodulation_out,pilot_inter,pilot_code,pilot_num);%����LS�����㷨��õ��Ľ����ź�
no_detect_sig=de_p(ofdm_demodulation_out,pilot_inter,pilot_code,pilot_num);%������Ӌ�r��õ��Ľ����ź�,ֱ�Ӱѵ��lɾȥ
%��ӳ��
ls_receive_bit_sig=de_map(ls_zf_detect_sig);%16QAM��ӳ��
no_est_receive_bit_sig=de_map(no_detect_sig);%16QAM��ӳ��
    %���¹���ͳ�Ƹ��ֹ����㷨�õ��Ľ����ź��еĴ��������
    ls_err_num=error_count(convert_matrix,ls_receive_bit_sig);
    no_est_err_num=error_count(convert_matrix,no_est_receive_bit_sig);
    
    ls_error_bit=ls_error_bit+ls_err_num;
    no_est_error_bit=no_est_error_bit+no_est_err_num;
%������ֹ����㷨���������
ls_err_ber(i)=ls_error_bit/total_bit_num;
no_est_error_ber(i)=no_est_error_bit/total_bit_num;
end
plot(SNR_dB,no_est_error_ber,'b-*',SNR_dB,ls_err_ber,'r-o')
title('�ྶ����²��뵼Ƶ,ѵ�����ż��������Ʒ�ʽ�Ƚ�')
xlabel('SNR')
ylabel('BER')
grid on
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ofdmϵ�y�������  Ӗ������ �й�Ӌ 16qam  �ж��� 
%���x
%pilot_inter=5;
%pilot_symbol_bit=[0 0 0 1];
carrier_count=128;
bits_per_symbol_16QAM=4;
cp_length=16;
SNR_dB=[0 2 4 6 8 10 12 14 16];
training_err_ber=zeros(1,length(SNR_dB));
for i=1:length(SNR_dB) %ÿ��SNR���Ϸ������ɴ�
    training_err_bit=0;
    total_bit_num=0;
symbols_per_carrier=50;
baseband_out_num=carrier_count*(symbols_per_carrier)*bits_per_symbol_16QAM;
baseband_out=round(rand(1, baseband_out_num));   
%����������ӳ��Ϊ16������
convert_matrix=reshape(baseband_out, bits_per_symbol_16QAM*carrier_count,length(baseband_out)/ (bits_per_symbol_16QAM*carrier_count));
total_bit_num=baseband_out_num;
%Ӱ��
map_out=map_16_QAM(convert_matrix);
%�����뵼�l
%[insert_pilot_out,pilot_num,pilot_code]=insert_pilot(pilot_inter,pilot_symbol_bit,map_out);
%����Ӗ������
training_symbols=[1+j;1-j;-1+j;-1-j;3+j;1+3*j;3-j;1-3*j;-1+3*j;-3+j;-3-j;-1-3*j;3+3*j;-3+3*j;-3-3*j;3-3*j];
training_symbols=[training_symbols;training_symbols];
training_symbols=[training_symbols;training_symbols];
training_symbols=[training_symbols;training_symbols];
training_symbols=[training_symbols,training_symbols];
training_symbols=[training_symbols,training_symbols];   %���һ��128*4��Ӗ����̖���
insert_training_out=[training_symbols,map_out];    %�Ͼ���Ҫ���Ĕ�����
train_length=size(training_symbols,2);
total_length=size(insert_training_out,2);
%FFT׃�Q
ofdm_mode_out=ifft(insert_training_out,128);
%����ѭ��ǰ׺
ofdm_cp_out=insert_cp(ofdm_mode_out,cp_length);
[ofdm_cp_out_m,ofdm_cp_out_n]=size(ofdm_cp_out);
%�ͳ�
ofdm_modulation=reshape(ofdm_cp_out,1,ofdm_cp_out_m*ofdm_cp_out_n);
Tx_data=ofdm_modulation;
%��̬�ŵ����ྶ��
multipath_signal=multi_chan(Tx_data,ofdm_cp_out_m,ofdm_cp_out_n);
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
    sgma=sqrt(spow1/(2*snr));%sgma���㷽ʽ���뵱ǰSNR���ź�ƽ�������й�ϵ
    receive_ofdm_symbol=add_noise(sgma,passchan_ofdm_symbol);%���������˹������,receive_ofdm_symbolΪ���ս��ջ��յ���ofdm���ſ�
%ȥѭ�hǰ�Y
cut_cp_symbol=del_cp(receive_ofdm_symbol,cp_length);
%FFT
ofdm_demodulation_out=fft(cut_cp_symbol,128);%��128��FFT���㣬���ofdm���

Rx_carriers = ofdm_demodulation_out(:,(train_length+1):total_length );
Rx_training_symbols = ofdm_demodulation_out(:,(1: train_length)  );
Rx_training_symbols1=Rx_training_symbols;
%�ŵ�����
Rx_training_symbols = Rx_training_symbols./ training_symbols;
Rx_training_symbols_deno = Rx_training_symbols.^2;
Rx_training_symbols_deno = Rx_training_symbols_deno(:,1)+Rx_training_symbols_deno(:,2)+Rx_training_symbols_deno(:,3)+Rx_training_symbols_deno(:,4) ;
Rx_training_symbols_nume = Rx_training_symbols(:,1) +Rx_training_symbols(:,2) + Rx_training_symbols(:,3) +Rx_training_symbols(:,4) ;
Rx_training_symbols_nume = conj(Rx_training_symbols_nume) ;
% ȡ4�������ĵ�Ƶ������Ϊ�˽���ƽ���Ż�
% ������� ������������������OFDM���� ���в���
% ԭ��Ѱ��1/H����FFT֮������ݽ���Ƶ�򲹳�
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
%��ӳ��
    receive_bit_sig=de_map(Rx_carriers);
    %���¹���ͳ�Ƹ��ֹ����㷨�õ��Ľ����ź��еĴ��������
    err_num=error_count(convert_matrix,receive_bit_sig);
    trainning_error_bit=training_err_bit+err_num;
%������ֹ����㷨���������
trainning_err_ber(i)=trainning_error_bit/total_bit_num;
end
plot(SNR_dB,trainning_err_ber,'g-+')
legend('��������','���뵼Ƶ(ls׼��)����','����ѵ�����ŵĹ���')
grid on
hold on
