% Hybrid beamforming test (all-array)
clear, clc


% 전송 데이터 스펙
fft_len = 64;
cp_len = fft_len / 4;
mod_type = 4;
data_len = fft_len * mod_type;


% 송신단 및 수신단에서의 총 데이터 스트림의 수
N_s = 2;
N_d = 2;
A = [1] ;  % BD 사용 시 수신기 당 안테나 수

% 송수신 각 배열 안테나의 안테나 수
ant_num = [4 1 0.5; 1 1 0.5];
N_tx = ant_num(1,1) * ant_num(1,2);
N_rx = ant_num(2,1) * ant_num(2,2);


% ber 결과 출력 벡터들
snr = 0:2:30;
ber_result = zeros(1, length(snr) );
iter = 100;


% 채널 관련 파라미터 초기화
path = 7;
H = zeros(path, fft_len+cp_len, N_rx * N_d, N_tx * N_s);
rx_angle = zeros(4, path, N_d, N_s);
t_He = zeros(path, N_d, N_s);
He = zeros(fft_len, N_d, N_s); 
m_path = zeros(1,N_d);

for i = 1:iter
    
%         tic
    
    for j = 1:length(snr)
        
        %============= 채널 계산 ============================================
        
        % 채널 계수 생성
        for d = 1:N_d
            for s = 1:N_s
                [temp, rx_angle(:,:,d,s)] = FD_channel(fft_len + cp_len, path, ant_num);
                H(:,:,1 + (d - 1) * N_rx : d * N_rx, 1 + (s - 1) * N_tx : s * N_tx) = temp;
                if d == s, [~, m_path(d)] = max( abs( temp(:,1,1,1) ) ); end
            end
        end
        
        % Beamforming 계수 계산 (송신, 수신 계수)
        Wt = beamformer(ant_num(1,:), rx_angle(1:2,:,:,:), N_s, 2, m_path);
        Wr = beamformer(ant_num(2,:), rx_angle(3:4,:,:,:), N_d, 2, m_path);
        
        % 시간 영역 Effective 채널 계수 계산
        temp = zeros(N_rx * N_d, N_tx * N_s);
        for k = 1:path
            temp(:,:) = H(k,1,:,:);
            t_He(k,:,:)= Wr.' * temp * Wt;
        end
        
        % 주파수 영역 변환
        for d = 1:N_d
            for s = 1:N_s
                He(:,d,s) = fft(t_He( :,d,s), fft_len);
            end
        end
        
        %============== 신호 송신 ===========================================
        
        % 전송 심볼 생성 및 precoding 수행
        bit = randi([0 1], N_s, data_len);
        Dsym = base_mod(bit, mod_type);
                  [Dsym, ~, Wd] = ZF_precoding(Dsym, He);
%                 [Dsym, ~, Wd] = MMSE_precoding(Dsym, He);
%         [Dsym, ~, bd_H, Wd] = BD_precoding(Dsym, He, A);
        
        % 정규화 상수 계산
         factor = zeros(1,fft_len);
        for k = 1:fft_len
            t_Wd(:,:) = Wd(k,:,:);
            temp = Wt * t_Wd;
            factor(k) = sqrt( trace( temp * temp' ) );
        end 
        
        % OFDM 심볼 생성
        Dsym = Dsym ./ factor;
        Isym = ifft(Dsym, fft_len, 2) * sqrt(fft_len);
        tx_ofdm = [ Isym(:, fft_len - cp_len + 1 : end) Isym ];
        
        % 송신 빔계수 적용
        tx_ofdm = Wt * tx_ofdm;
        
        %============== 신호 수신 ===========================================
        
        
        % 채널 통과
        rx_ofdm = awgn_noise( FD_fading( tx_ofdm, H ), snr(j) );
        
        % 수신 빔포밍 계수 적용
        rx_ofdm = Wr.' * rx_ofdm;
        
        % OFDM 복조
        rx_Isym = rx_ofdm(:, cp_len + 1 : fft_len + cp_len);
        rx_Dsym = fft(rx_Isym, fft_len, 2) / sqrt(fft_len);
        rx_Dsym = rx_Dsym .* factor;
        
        
        % BD p recoding이 사용되었을 경우 수신기별 MIMO 검출 수행
        if length(A) > 1
            N = 1;   % 채널 행렬에서 각 수신기의 첫 번째 인덱스
            for n = 1:length(A)
                idx = N : N + A(n) - 1;  % 해당 수신기의 대각 인덱스 범위
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

% 결과 평균화
ber_result = ber_result / iter;

% 결과 출력
semilogy(snr, ber_result, 'r-');
title('Hybrid Beamforming')
legend('Detection result')
ylabel('BER')
xlabel('SNR (dB)')
grid on
axis([0 max(snr) 10^-5 1])
