% Sepideh algorithm test (fully-connected)

% ���� ������ ���� ����
model = SCM();
model.n_mray = scatter;
model.n_path = path;
cp_len = fft_len / 4;
data_len = fft_len * mod_type;

% �۽Ŵ� �� ���Ŵܿ����� �� ������ ��Ʈ�� ����
N_s = rx_node;
N_d = rx_node;
N_rf = num_rf;

% �ۼ��� �� �迭 ���׳��� ���׳� ��
model.ant(rx_ant, tx_ant);
N_tx = model.Ntx;
N_rx = model.Nrx;

% ber ��� ��� ���͵�
result_ber = zeros(1, length(snr) );
result_cap = zeros(1, length(snr) );

% ä�� ���� �Ķ���� �ʱ�ȭ
H = zeros(path, fft_len+cp_len, N_rx * N_d, N_tx);
t_He = zeros(path, N_d, N_rf);
He = zeros(fft_len, N_d, N_rf);


for i = 1:iter

%     tic
    for j = 1:length(snr)
        
        %============= ä�� ��� ============================================
        
        % ä�� ��� ����
        for d = 1:N_d
            [temp, rx_angle] = model.FD_channel(fft_len + cp_len);
            H(:,:,1+(d-1)*N_rx:d*N_rx,:) = temp;
        end
        
        % �ð� ä�� ���ļ� ��ȯ
        tmp_H(:,:,:) = H(:,1,:,:);
        H_f = fft(tmp_H, fft_len, 1);        
        
        % Beamforming ��� ��� (�۽�, ���� ���)
        Wt = sep_precoding(H_f, N_rf);
        
        %============== ��ȣ �۽� ===========================================
        
        % ���� �ɺ� ����
        bit = randi([0 1], N_s, data_len);
        Dsym = base_mod(bit, mod_type);
        
        % precoding ����
        for k = 1:fft_len
            tmp_H_f(:,:) = H_f(k,:,:);
            He(k,:,:) = tmp_H_f * Wt;
        end
        [Dsym, ~, Wd] = ZF_precoding(Dsym, He);
        
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
        % ========================= BER ��� ==============================
        % ä�� ���
        [rx_ofdm, No] = awgn_noise( model.FD_fading( tx_ofdm, H ), snr(j) );
        
        % ���� ������ ��� ����
        rx_ofdm = rx_ofdm;
        
        % OFDM ����
        rx_Isym = rx_ofdm(:, cp_len + 1 : fft_len + cp_len);
        rx_Dsym = fft(rx_Isym, fft_len, 2) / sqrt(fft_len);
        rx_Dsym = rx_Dsym .* factor;
        
        % base demodulation
        rx_bit = base_demod(rx_Dsym, mod_type);
        result_ber(j) = result_ber(j) + sum( sum( bit ~= rx_bit, 2 ), 1 ) / (data_len * N_s);
        
        
        % ======================== Capacity ��� ==========================
        
        % Capacity ���
        rx_ofdm = model.FD_fading( tx_ofdm, H );
        rx_Isym = rx_ofdm(:, cp_len + 1 : fft_len + cp_len);
        rx_Dsym = fft(rx_Isym, fft_len, 2) / sqrt(fft_len);
        result_cap(j) = result_cap(j) + mean( sum( log2( 1 + abs(rx_Dsym).^2 / (No * rx_ant) ) ) );
        
        
    end
%     toc
    
end


% ��� ���ȭ
result_ber = result_ber / iter;
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
    semilogy(snr, result_ber, plot_format);
    title('BER Performance')
    legend('Detection result')
    ylabel('BER')
    xlabel('SNR (dB)')
    grid on
    axis([snr(1) max(snr) 10^-5 1])
    
elseif mode == 3
    % Repeatition ��� ���
    plot(snr, result_rep, plot_format);
    title('Repeatition Performance')
    ylabel('Repeatition Number')
    xlabel('SNR (dB)')
    grid on

end
