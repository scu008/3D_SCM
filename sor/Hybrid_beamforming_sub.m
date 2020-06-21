% Hybrid beamforming test (all-array)
clear, clc


% ���� ������ ����
fft_len = 64;
cp_len = fft_len / 4;
mod_type = 4;
data_len = fft_len * mod_type;


% �۽Ŵ� �� ���Ŵܿ����� �� ������ ��Ʈ���� ��
N_s = 2;
N_d = 2;
A = [1] ;  % BD ��� �� ���ű� �� ���׳� ��

% �ۼ��� �� �迭 ���׳��� ���׳� ��
ant_num = [4 1 0.5; 1 1 0.5];
N_tx = ant_num(1,1) * ant_num(1,2);
N_rx = ant_num(2,1) * ant_num(2,2);


% ber ��� ��� ���͵�
snr = 0:2:30;
ber_result = zeros(1, length(snr) );
iter = 100;


% ä�� ���� �Ķ���� �ʱ�ȭ
path = 7;
H = zeros(path, fft_len+cp_len, N_rx * N_d, N_tx * N_s);
rx_angle = zeros(4, path, N_d, N_s);
t_He = zeros(path, N_d, N_s);
He = zeros(fft_len, N_d, N_s); 
m_path = zeros(1,N_d);

for i = 1:iter
    
%         tic
    
    for j = 1:length(snr)
        
        %============= ä�� ��� ============================================
        
        % ä�� ��� ����
        for d = 1:N_d
            for s = 1:N_s
                [temp, rx_angle(:,:,d,s)] = FD_channel(fft_len + cp_len, path, ant_num);
                H(:,:,1 + (d - 1) * N_rx : d * N_rx, 1 + (s - 1) * N_tx : s * N_tx) = temp;
                if d == s, [~, m_path(d)] = max( abs( temp(:,1,1,1) ) ); end
            end
        end
        
        % Beamforming ��� ��� (�۽�, ���� ���)
        Wt = beamformer(ant_num(1,:), rx_angle(1:2,:,:,:), N_s, 2, m_path);
        Wr = beamformer(ant_num(2,:), rx_angle(3:4,:,:,:), N_d, 2, m_path);
        
        % �ð� ���� Effective ä�� ��� ���
        temp = zeros(N_rx * N_d, N_tx * N_s);
        for k = 1:path
            temp(:,:) = H(k,1,:,:);
            t_He(k,:,:)= Wr.' * temp * Wt;
        end
        
        % ���ļ� ���� ��ȯ
        for d = 1:N_d
            for s = 1:N_s
                He(:,d,s) = fft(t_He( :,d,s), fft_len);
            end
        end
        
        %============== ��ȣ �۽� ===========================================
        
        % ���� �ɺ� ���� �� precoding ����
        bit = randi([0 1], N_s, data_len);
        Dsym = base_mod(bit, mod_type);
                  [Dsym, ~, Wd] = ZF_precoding(Dsym, He);
%                 [Dsym, ~, Wd] = MMSE_precoding(Dsym, He);
%         [Dsym, ~, bd_H, Wd] = BD_precoding(Dsym, He, A);
        
        % ����ȭ ��� ���
         factor = zeros(1,fft_len);
        for k = 1:fft_len
            t_Wd(:,:) = Wd(k,:,:);
            temp = Wt * t_Wd;
            factor(k) = sqrt( trace( temp * temp' ) );
        end 
        
        % OFDM �ɺ� ����
        Dsym = Dsym ./ factor;
        Isym = ifft(Dsym, fft_len, 2) * sqrt(fft_len);
        tx_ofdm = [ Isym(:, fft_len - cp_len + 1 : end) Isym ];
        
        % �۽� ����� ����
        tx_ofdm = Wt * tx_ofdm;
        
        %============== ��ȣ ���� ===========================================
        
        
        % ä�� ���
        rx_ofdm = awgn_noise( FD_fading( tx_ofdm, H ), snr(j) );
        
        % ���� ������ ��� ����
        rx_ofdm = Wr.' * rx_ofdm;
        
        % OFDM ����
        rx_Isym = rx_ofdm(:, cp_len + 1 : fft_len + cp_len);
        rx_Dsym = fft(rx_Isym, fft_len, 2) / sqrt(fft_len);
        rx_Dsym = rx_Dsym .* factor;
        
        
        % BD p recoding�� ���Ǿ��� ��� ���ű⺰ MIMO ���� ����
        if length(A) > 1
            N = 1;   % ä�� ��Ŀ��� �� ���ű��� ù ��° �ε���
            for n = 1:length(A)
                idx = N : N + A(n) - 1;  % �ش� ���ű��� �밢 �ε��� ����
%                   rx_Dsym(idx,:) = OSIC_detect( rx_Dsym(idx,:), bd_H(:,idx,idx),mod_type,1 );
                 rx_Dsym(idx,:) = ZF_detect( rx_Dsym(idx,:), bd_H(:,idx,idx));
%                   rx_Dsym(idx,:) = MMSE_detect( rx_Dsym(idx,:), bd_H(:,idx,idx),1);
%                   rx_Dsym(idx,:) = ML_detect( rx_Dsym(idx,:), bd_H(:,idx,idx),mod_type);
%                   rx_Dsym(idx,:) = DFE_detect( rx_Dsym(idx,:), bd_H(:,idx,idx),mod_type,1,1); % mode - 1: ZF, 2: MMSE
                N = N + A(n);
            end
        end
        
        % base demodulation
        rx_bit = base_demod(rx_Dsym, mod_type);
        ber_result(j) = ber_result(j) + sum( sum( bit ~= rx_bit, 2 ), 1 ) / (data_len * N_s);
    end
    
%         toc
end

% ��� ���ȭ
ber_result = ber_result / iter;

% ��� ���
semilogy(snr, ber_result, 'r-');
title('Hybrid Beamforming')
legend('Detection result')
ylabel('BER')
xlabel('SNR (dB)')
grid on
axis([0 max(snr) 10^-5 1])
