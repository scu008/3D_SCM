% Hybrid beamforming test (all-array)
clear, clc
% 전송 데이터 스펙
fft_len = 128;
cp_len = fft_len / 4;
mod_type = 2;
data_len = fft_len * mod_type;
cl = 3;
code_rate = 1/2;
TRELLIS = poly2trellis(cl,codegenerator(code_rate,3));

% 송신단 및 수신단에서의 총 데이터 스트림의 수
N_s = 4;
N_d = 4;
N_rf = 4;
A = [1];  % BD 사용 시 수신기 당 안테나 수

% 송수신 각 배열 안테나의 안테나 수
ant_num = [32 1 0.5; 1 1 0.5];
N_tx = ant_num(1,1) * ant_num(1,2);
N_rx = ant_num(2,1) * ant_num(2,2);

% ber 결과 출력 벡터들
snr = -10:5:15;
ber_result = zeros(1, length(snr) );
C_simulation = zeros(1,length(snr));

iter = 30;

% 채널 관련 파라미터 초기화
path = 7;
H = zeros(path, fft_len+cp_len, N_rx * N_d, N_tx);
rx_angle = zeros(4, path, N_d, N_s);
t_He = zeros(path, N_d, N_rf);
He = zeros(fft_len, N_d, N_s);
m_path = zeros(1,N_d);

for j = 1:length(snr)    
        tic    
    for  i = 1:iter
        
        %============= 채널 계산 ============================================
        
        % 채널 계수 생성
        for d = 1:N_d
            [temp, rx_angle(:,:,d,d)] = FD_channel(fft_len + cp_len, path, ant_num);
            H(:,:,1+(d-1)*N_rx:d*N_rx,:) = temp;
            [~, m_path(d)] = max( abs( temp(:,1,1,1) ) );
        end
        
        % Beamforming 계수 계산 (송신, 수신 계수)
%         Wt = beamformer(ant_num(1,:), rx_angle(1:2,:,:,:), N_s, 1, m_path);
        [Wt, Wd] = sor_beamformer(ant_num,N_rf, N_d, H, fft_len, snr(j));
%         Wt = sor_beamformer2(ant_num, N_s, N_rf, H, fft_len);
        Wr = beamformer(ant_num(2,:), rx_angle(3:4,:,:,:), N_d, 2);
        
        % 시간 영역 Effective 채널 계수 계산
        temp = zeros(N_rx * N_d, N_tx);
        for k = 1:path
            temp(:,:) = H(k,1,:,:);
%             t_He(k,:,:)= Wr.' * temp * Wt;
            t_He(k,:,:)= temp * Wt;
        end
        % 주파수 영역 변환
        for d = 1 : N_d
            for s = 1 : N_rf
                He(:,d,s) = fft(t_He(:,d,s), fft_len,1);                
            end
        end
        
        
        %============== 신호 송신 ===========================================
        
        % 전송 심볼 생성 및 precoding 수행
        bit = randi([0 1], N_s, data_len * code_rate);
        bit_coded = zeros(N_s, data_len);
        bit_inter = zeros(N_s, data_len);
        for n = 1 : N_s
            bit_coded(n,:) = convenc(bit(n,:), TRELLIS); % coded
            bit_inter(n,:) = interleaver(bit_coded(n,:), 1);
        end
        Dsym = base_mod(bit_inter, mod_type);
        Dsym_ = zeros(N_rf, fft_len);
        Heff = zeros(fft_len, N_s, N_s);
        for k = 1 : fft_len
            % Digital Precoding 수행
            Dsym_(:,k) = squeeze(Wd(k,:,:)) * Dsym(:,k);
            
            % effective channel
            Heff(k,:,:) = squeeze(He(k,:,:)) * squeeze(Wd(k,:,:));
        end
        Heff = squeeze(mean(Heff));
        
%         Dsym = Wd * Dsym;
%           [Dsym, ~, Wd] = ZF_precoding(Dsym, He);
%         [Dsym, ~, bd_H, Wd] = BD_precoding(Dsym, He, A);
        
        % 정규화 상수 계산
        factor = zeros(1,fft_len);
        for k = 1:fft_len
            t_Wd(:,:) = Wd(k,:,:);
            temp = Wt * t_Wd;
            factor(k) = sqrt( trace( temp * temp' ) );
        end
        
        % OFDM 심볼 생성
        Dsym = Dsym_ ./ factor;
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
        
        
        % BD precoding이 사용되었을 경우 수신기별 MIMO 검출 수행
        if length(A) > 1
            N = 1;   % 채널 행렬에서 각 수신기의 첫 번째 인덱스
            for n = 1:length(A)
                idx = N : N + A(n) - 1;  % 해당 수 신기의 대각 인덱스 범위
                rx_Dsym(idx,:) = ZF_detect( rx_Dsym(idx,:), bd_H(:,idx,idx) );
                N = N + A(n);
            end
        end
        
        
        % base demodulation
        rx_bit = base_demod(rx_Dsym, mod_type);
        bit_deinter = zeros(N_s, data_len);
        bit_decoded = zeros(N_s, data_len * code_rate);
        for n = 1 : N_s
            bit_deinter(n,:) = interleaver(rx_bit(n,:), 2);
            bit_decoded(n,:) = vitdec(bit_deinter(n,:), TRELLIS, cl*5,'trunc','hard');
        end
        ber_result(j) = ber_result(j) + sum( sum( bit ~= bit_decoded, 2 ), 1 ) / (data_len * code_rate * N_s);
        C_simulation(j) = C_simulation(j) + log2(det(eye(N_s) + 10^(snr(j)/10)*(abs(Heff*Heff'))))/iter;
    end
    
        toc
end

% 결과 평균화
ber_result = ber_result / iter;

% 결과 출력
figure(1)
semilogy(snr, ber_result, 'r-d');
title('Hybrid Beamforming')
legend('Detection result')
ylabel('BER')
xlabel('SNR (dB)')
grid on
axis([min(snr) max(snr) 10^-5 1])

figure(2)
plot(snr, C_simulation,'-o');
movegui('northeast');
title('SumRate'), xlabel('SNR'), ylabel('SumRate')
legend('.'), grid on
