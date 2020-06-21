function [Vrf, Vd, sum_rate] = sor_beamformer(ant_num, Nrf, N_d, H, fft_len, SNR)
% %% 전송 데이터 스펙
% clear
% fft_len = 128;
% cp_len = fft_len / 4;
% mod_type = 2;
% data_len = fft_len * mod_type;
% path = 7;
% % 송신단 및 수신단에서의 총 데이터 스트림의 수
% N_s = 4;
% N_d = 4;
% % 송수신 각 배열 안테나의 안테나 수
% ant_num = [32 1 0.5; 1 1 0.5];
% N_tx = ant_num(1,1) * ant_num(1,2);
% N_rx = ant_num(2,1) * ant_num(2,2);
% 
% % 채널 계수 생성
% for d = 1:N_d
%     [temp, ~ ] = FD_channel(fft_len + cp_len, path, ant_num);
%     H(:,:,1+(d-1)*N_rx:d*N_rx,:) = temp;
% end
% Nrf = 4;
% SNR = 15;
%% Analog beamforming
sum_rate = 0;
SNR = 10^(SNR/10);
sigma = 1;
N = ant_num(1,1) * ant_num(1,2);
H = squeeze(H(:,1,:,:));
H_f = fft(H,fft_len,1);
F1 = zeros(fft_len, N, N);
for k = 1 : fft_len
    temp = squeeze(H_f(k,:,:));
    F1(k,:,:) =  temp'*temp;
end
F1 = squeeze(sum(F1,1)/fft_len);
 
Vrf = ones(N, Nrf);
gamma_squ = SNR/(N*Nrf);
iter = 1;
count1 = 0;

while (iter)
    temp = Vrf;
    for j = 1 : Nrf
        Vj = Vrf;
        Vj(:,j) = [];
        C = eye(Nrf - 1) + gamma_squ * Vj' * F1 * Vj;
        G = gamma_squ * F1 - gamma_squ^2 * F1 * Vj / C * Vj' * F1;
        for i = 1 : N
            eta = G(i,:)*Vrf(:,j) - G(i,i)*Vrf(i,j);
            if abs(eta) < 0.0001
                Vrf(i,j) = 1;
            else
                Vrf(i,j) = eta/abs(eta);
            end
        end
    end
    if mean(abs(temp - Vrf),'all') < 0.01
        iter = 0;
    end
    count1 = count1 + 1;
end
temp = zeros(N_d, N);
t_He = zeros(size(H,1), N_d, Nrf);
for k = 1 : size(H,1)
    temp(:,:) = H(k,:,:);
    t_He(k,:,:)= temp * Vrf;
end
%% digital beamforming
Vd = zeros(fft_len,Nrf, N_d);
for k = 1 : fft_len
    count2 = 0;
    iter = 1;
    h = squeeze(H_f(k,:,:));
    g_ = h*Vrf;
    Vd_k = g_' /( g_ * g_' );
    Vd_k = Vd_k/sqrt(trace(Vrf*Vd_k*Vd_k'*Vrf')) * sqrt(SNR);
    w = zeros(N_d, 1);
    e = zeros(N_d, 1);
    g_ = h * Vrf;
    while(iter)
        for n = 1 : N_d
            temp = Vd_k;
            w(n) = g_(n,:)*temp(:,n) / (sigma + sum(abs(g_(n,:)*temp).^2));
            e(n) = sum( (abs(w(n)*g_(n,:)*temp).^2 ) - 2 * real((w(n))*g_(n,:)*temp) + abs(w(n)).^2 * sigma + 1);
        end
        t = 1./e;
        temp2 = 0;
        for n = 1 : N_d
            temp2 = temp2 +  t(n) * abs(w(n)).^2 * g_(n,:)' * g_(n,:);  
        end
        for n = 1 : N_d
            lambda = sum(t(n) * abs(w.^2)) / SNR ;
            J = temp2 + lambda * (Vrf'*Vrf);
            Vd_k(:,n) = t(n) * abs(w(n)) * inv(J) * g_(n,:)';  
        end
        
        a = abs(temp - Vd_k);
        count2 = count2 + 1;
        if (a < 0.001)
            iter = 0;
        end
    end
%     SNR
%     trace((Vrf*(Vd_k*Vd_k')*Vrf'))
%     g_*Vd_k
    
    Vd(k,:,:) = Vd_k;
    
    for n = 1 : N_d
        sum_rate = sum_rate + log2(1 + trace(abs(h*Vrf*Vd_k).^2) / (sigma + sum(abs(h*Vrf*Vd_k).^2,'all') - trace(abs(h*Vrf*Vd_k).^2)));
    end
end
sum_rate = sum_rate / fft_len;