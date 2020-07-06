% MU-MIMO test

% �⺻ �Ķ���� ����
model = SCM();
model.n_path = path;
model.n_mray = scatter;
cp_len = fft_len / 4;
data_len = fft_len * mod_type;

% ä�� ȯ��
N_d = rx_node; 
N_s = rx_node;
A = [1];

% �ۼ��� �� �迭 ���׳��� ���׳� �� (����, �۽�)
model.ant(1, tx_ant);
N_tx = model.Ntx;
N_rx = model.Nrx;

% �⺻ �Ķ���� ����
r_H = zeros(path, fft_len + cp_len, N_rx * N_d, N_tx);
t_H = zeros(path, N_rx * N_d, N_tx);
eq_sym = zeros(N_d, fft_len);
result_mimo = zeros(1, length(snr) );
result_cap = zeros(1, length(snr) );

% �ݺ� ����
tic
for i = 1 : iter    
    % ä�� ��� ����
    for d = 1:N_d
        temp = model.FD_channel(fft_len + cp_len);
        r_H(:,:,1+(d-1)*N_rx:d*N_rx,:) = temp;
    end
    t_H(:,:,:) = r_H(:,1,:,:);
    H = fft(t_H, fft_len, 1);
    
    % ������ ����
    bit = randi([0 1], N_s, data_len);
    
    % modulation
    sym = base_mod(bit, mod_type);    
    
    % SNR�� ���� �ݺ�
    for j = 1 : length(snr)
        
        % ���� ���� �ʱ�ȭ
        [~, No] = awgn_noise(zeros(1), snr(j));
        
        
        % Precoding ����
%         [sym_hat, K, rx_H] = BD_precoding(sym, H, A);
        [sym_hat, K] = ZF_precoding(sym, H);
%         [sym_hat, W, K, U] = SVD_precoding(sym, H, No);
        sym_hat = sym_hat ./ K;
        
        
        % ä�� ���
        rx_sym = zeros(N_d, fft_len);
        for r = 1 : N_d
           for t = 1 : N_tx
               rx_sym(r,:) = rx_sym(r,:) + sym_hat(t,:) .* H(:,r,t).';
           end
        end
        
        % ��ȣ ����
        [eq_sym, No] = awgn_noise( rx_sym, snr(j) );
        eq_sym = eq_sym .* K;    
        
        
        % SVD �� water filling ����
%         for k = 1:fft_len
%            t_W(:,:) = U(k,:,:);
%            rx_sym(:,k) = t_W' * rx_sym(:,k);     
%         end

        
        result_cap(j) = result_cap(j) + mean( sum( log2( 1 + abs(rx_sym).^2 / No ) ) );
        
        
        % BD precoding�� ���Ǿ��� ��� ���ű⺰ MIMO ���� ����
        if length(A) > 1
           N = 1;   % ä�� ��Ŀ��� �� ���ű��� ù ��° �ε���
           for n = 1:length(A)
               idx = N : N + A(n) - 1;  % �ش� ���ű��� �밢 �ε��� ����
               eq_sym(idx,:) = ML_detect( eq_sym(idx,:), rx_H(:,idx,idx), mod_type );
               N = N + A(n);
           end
        end
       
        
        % demodulation
        rx_bit = base_demod(eq_sym, mod_type);
        
        % ��� ����
        result_mimo(j) = result_mimo(j) + sum( sum( bit ~= rx_bit ) )/ (data_len * N_s);
        
    end
end
toc

% ��� ���
result_mimo = result_mimo / iter;
result_cap = result_cap / iter;

if mode == 1
    % Capacity ��� ���
    plot(snr, result_cap, plot_format);
    title('Sum Rate Performance')
    legend('Sum Capacity')
    ylabel('Average Spectral Efficiency (bps/Hz)')
    xlabel('SNR (dB)')
    grid on
    
elseif mode == 2
    % BER ��� ���
    semilogy(snr, result_mimo, plot_format);
    title('BER Performance')
    legend('Detection result')
    ylabel('BER')
    xlabel('SNR (dB)')
    grid on
    axis([snr(1) max(snr) 10^-5 1])
end
