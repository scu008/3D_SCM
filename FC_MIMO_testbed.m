% MU-MIMO test

% 기본 파라미터 설정
model = SCM();
model.n_path = path;
model.n_mray = scatter;
cp_len = fft_len / 4;
data_len = fft_len * mod_type;

% 채널 환경
N_d = rx_node; 
N_s = rx_node;
A = [1];

% 송수신 각 배열 안테나의 안테나 수 (수신, 송신)
model.ant(1, tx_ant);
N_tx = model.Ntx;
N_rx = model.Nrx;

% 기본 파라미터 설정
r_H = zeros(path, fft_len + cp_len, N_rx * N_d, N_tx);
t_H = zeros(path, N_rx * N_d, N_tx);
eq_sym = zeros(N_d, fft_len);
result_mimo = zeros(1, length(snr) );
result_cap = zeros(1, length(snr) );

% 반복 시작
tic
for i = 1 : iter    
    % 채널 계수 생성
    for d = 1:N_d
        temp = model.FD_channel(fft_len + cp_len);
        r_H(:,:,1+(d-1)*N_rx:d*N_rx,:) = temp;
    end
    t_H(:,:,:) = r_H(:,1,:,:);
    H = fft(t_H, fft_len, 1);
    
    % 데이터 생성
    bit = randi([0 1], N_s, data_len);
    
    % modulation
    sym = base_mod(bit, mod_type);    
    
    % SNR에 따라 반복
    for j = 1 : length(snr)
        
        % 수신 벡터 초기화
        [~, No] = awgn_noise(zeros(1), snr(j));
        
        
        % Precoding 수행
%         [sym_hat, K, rx_H] = BD_precoding(sym, H, A);
        [sym_hat, K] = ZF_precoding(sym, H);
%         [sym_hat, W, K, U] = SVD_precoding(sym, H, No);
        sym_hat = sym_hat ./ K;
        
        
        % 채널 통과
        rx_sym = zeros(N_d, fft_len);
        for r = 1 : N_d
           for t = 1 : N_tx
               rx_sym(r,:) = rx_sym(r,:) + sym_hat(t,:) .* H(:,r,t).';
           end
        end
        
        % 신호 수신
        [eq_sym, No] = awgn_noise( rx_sym, snr(j) );
        eq_sym = eq_sym .* K;    
        
        
        % SVD 및 water filling 사용시
%         for k = 1:fft_len
%            t_W(:,:) = U(k,:,:);
%            rx_sym(:,k) = t_W' * rx_sym(:,k);     
%         end

        
        result_cap(j) = result_cap(j) + mean( sum( log2( 1 + abs(rx_sym).^2 / No ) ) );
        
        
        % BD precoding이 사용되었을 경우 수신기별 MIMO 검출 수행
        if length(A) > 1
           N = 1;   % 채널 행렬에서 각 수신기의 첫 번째 인덱스
           for n = 1:length(A)
               idx = N : N + A(n) - 1;  % 해당 수신기의 대각 인덱스 범위
               eq_sym(idx,:) = ML_detect( eq_sym(idx,:), rx_H(:,idx,idx), mod_type );
               N = N + A(n);
           end
        end
       
        
        % demodulation
        rx_bit = base_demod(eq_sym, mod_type);
        
        % 결과 저장
        result_mimo(j) = result_mimo(j) + sum( sum( bit ~= rx_bit ) )/ (data_len * N_s);
        
    end
end
toc

% 결과 평균
result_mimo = result_mimo / iter;
result_cap = result_cap / iter;

if mode == 1
    % Capacity 결과 출력
    plot(snr, result_cap, plot_format);
    title('Sum Rate Performance')
    legend('Sum Capacity')
    ylabel('Average Spectral Efficiency (bps/Hz)')
    xlabel('SNR (dB)')
    grid on
    
elseif mode == 2
    % BER 결과 출력
    semilogy(snr, result_mimo, plot_format);
    title('BER Performance')
    legend('Detection result')
    ylabel('BER')
    xlabel('SNR (dB)')
    grid on
    axis([snr(1) max(snr) 10^-5 1])
end
