%Ofdmϵ�y�������   (��)�ྶ�����BPSK,4QAM��16QAM���Ʒ�ʽ�ıȽ�
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  BPSK���Ʒ�ʽ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
%���x
pilot_inter=5;%���l�g��
pilot_symbol_bit=1;%���l��̖
carrier_count=256;%���d����
%bits_per_symbol_16QAM=4;%����ԓ�{��ʱһ��̖ռλ��(BPSK����)
cp_length=16;%cp�L��
SNR_dB=[0 2 4 6 8 10 12 14 16];
ls_err_ber=zeros(1,length(SNR_dB));
for i=1:length(SNR_dB) %ÿ��SNR���Ϸ������ɴ�
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
%����������ӳ��Ϊ2������
convert_matrix=reshape(baseband_out, carrier_count,symbols_per_carrier);
total_bit_num=baseband_out_num;
%Ӱ��
%map_out=map_16_QAM(convert_matrix);
map_out=convert_matrix;
%���뵼�l
[insert_pilot_out,pilot_num,pilot_code]=insert_pilot(pilot_inter,pilot_symbol_bit,map_out);
%FFT׃�Q
ofdm_mode_out=ifft(insert_pilot_out,128);
%����ѭ��ǰ׺
ofdm_cp_out=insert_cp(ofdm_mode_out,cp_length);
[ofdm_cp_out_m,ofdm_cp_out_n]=size(ofdm_cp_out);
%�ͳ�(���]�Ͷ�������һ��)
ofdm_modulation=reshape(ofdm_cp_out,1,ofdm_cp_out_m*ofdm_cp_out_n);
Tx_data=ofdm_modulation;
%��̬�ŵ����ྶ��
multipath_signal=ofdm_cp_out;
%multipath_signal=multi_chan(Tx_data,ofdm_cp_out_m,ofdm_cp_out_n);%���ⶨ�����ŵ�������r��������ͷ��"%"̖
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
%********************** ���¾��ǶԽ���ofdm�źŽ����ŵ����ƺ��źż��Ĺ���************************
ls_zf_detect_sig=ls_estimation(ofdm_demodulation_out,pilot_inter,pilot_code,pilot_num);%����LS�����㷨��õ��Ľ����ź�
%��ӳ��
    ls_receive_bit_sig=de_map1(ls_zf_detect_sig);%BPSK��ӳ��
    %���¹���ͳ�Ƹ��ֹ����㷨�õ��Ľ����ź��еĴ��������
    ls_err_num=error_count(convert_matrix,ls_receive_bit_sig);    
    ls_error_bit=ls_error_bit+ls_err_num;
%end
%������ֹ����㷨���������
ls_err_ber(i)=ls_error_bit/total_bit_num;
end
plot(SNR_dB,ls_err_ber,'b-*')
title('�Ƕྶ�����BPSK,4QAM��16QAM���Ʒ�ʽ�ıȽ�');
legend('BPSK')
xlabel('SNR')
ylabel('BER')
grid on
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  16QAM���Ʒ�ʽ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���x
pilot_inter=5;
pilot_symbol_bit=[0 0 0 1];
carrier_count=256;
bits_per_symbol_16QAM=4;%%����ԓ�{��ʱһ��̖ռλ��(4)
cp_length=16;
SNR_dB=[0 2 4 6 8 10 12 14 16];
ls_err_ber=zeros(1,length(SNR_dB));
lmmse_err_ber=zeros(1,length(SNR_dB));
lr_lmmse_err_ber=zeros(1,length(SNR_dB));
for i=1:length(SNR_dB) %ÿ��SNR���Ϸ������ɴ�
    ls_error_bit=0;
    lmmse_error_bit=0;
    lr_lmmse_error_bit=0;
    total_bit_num=0;
symbols_per_carrier=50;
loop=5;
%for l=1:loop
baseband_out_num=carrier_count*(symbols_per_carrier)*bits_per_symbol_16QAM;
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
%�ͳ�(���]�Ͷ�������һ��)
ofdm_modulation=reshape(ofdm_cp_out,1,ofdm_cp_out_m*ofdm_cp_out_n);
Tx_data=ofdm_modulation;
%��̬�ŵ����ྶ��
multipath_signal=ofdm_cp_out;
%multipath_signal=multi_chan(Tx_data,ofdm_cp_out_m,ofdm_cp_out_n);%���ⶨ�����ŵ�������r��������ͷ��"%"̖
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
%��ӳ��
    ls_receive_bit_sig=de_map(ls_zf_detect_sig);%16QAM��ӳ��
    %���¹���ͳ�Ƹ��ֹ����㷨�õ��Ľ����ź��еĴ��������
    ls_err_num=error_count(convert_matrix,ls_receive_bit_sig);
    ls_error_bit=ls_error_bit+ls_err_num;
%������ֹ����㷨���������
ls_err_ber(i)=ls_error_bit/total_bit_num;
end
plot(SNR_dB,ls_err_ber,'r-o')
title('�Ƕྶ�����BPSK,4QAM��16QAM���Ʒ�ʽ�ıȽ�')
legend('���뵼�l(ls׼�y)')
xlabel('SNR')
ylabel('BER')
grid on
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  4QAM���Ʒ�ʽ��
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���x
pilot_inter=5;
pilot_symbol_bit=[0 0 0 1];
carrier_count=256;
bits_per_symbol_4QAM=2;%����ԓ�{��ʱһ��̖ռλ��(2)
cp_length=16;
SNR_dB=[0 2 4 6 8 10 12 14 16];
ls_err_ber=zeros(1,length(SNR_dB));
lmmse_err_ber=zeros(1,length(SNR_dB));
lr_lmmse_err_ber=zeros(1,length(SNR_dB));
for i=1:length(SNR_dB) %ÿ��SNR���Ϸ������ɴ�
    ls_error_bit=0;
    total_bit_num=0;
symbols_per_carrier=50;
loop=5;
%for l=1:loop
baseband_out_num=carrier_count*(symbols_per_carrier+pilot_inter)*bits_per_symbol_4QAM;
baseband_out=round(rand(1, baseband_out_num));   
%����������ӳ��Ϊ16������
convert_matrix=reshape(baseband_out, bits_per_symbol_4QAM*carrier_count,length(baseband_out)/ (bits_per_symbol_4QAM*carrier_count));
total_bit_num=baseband_out_num;
%Ӱ��
map_out=map_4_QAM(convert_matrix);
%���뵼�l
[insert_pilot_out,pilot_num,pilot_code]=insert_pilot(pilot_inter,pilot_symbol_bit,map_out);
%FFT׃�Q
ofdm_mode_out=ifft(insert_pilot_out,128);
%����ѭ��ǰ׺
ofdm_cp_out=insert_cp(ofdm_mode_out,cp_length);
[ofdm_cp_out_m,ofdm_cp_out_n]=size(ofdm_cp_out);
%�ͳ�(���]�Ͷ�������һ��)
ofdm_modulation=reshape(ofdm_cp_out,1,ofdm_cp_out_m*ofdm_cp_out_n);
Tx_data=ofdm_modulation;
%��̬�ŵ����ྶ��
multipath_signal=ofdm_cp_out;
%multipath_signal=multi_chan(Tx_data,ofdm_cp_out_m,ofdm_cp_out_n);%���ⶨ�����ŵ���
%����r��������ͷ��"%"̖
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
%��ӳ��
    ls_receive_bit_sig=de_map4(ls_zf_detect_sig);%16QAM��ӳ��
    %���¹���ͳ�Ƹ��ֹ����㷨�õ��Ľ����ź��еĴ��������
    ls_err_num=error_count(convert_matrix,ls_receive_bit_sig);    
    ls_error_bit=ls_error_bit+ls_err_num;
%������ֹ����㷨���������
ls_err_ber(i)=ls_error_bit/total_bit_num;
lmmse_err_ber(i)=lmmse_error_bit/total_bit_num;
lr_lmmse_err_ber(i)=lr_lmmse_error_bit/total_bit_num;
end
plot(SNR_dB,ls_err_ber,'g-+')
title('�Ƕྶ�����BPSK,4QAM��16QAM���Ʒ�ʽ�ıȽ�')
legend('BPSK','16QAM','4QAM')
hold on
