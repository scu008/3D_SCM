% function Vrf = sor_beamformer2(ant_num, K, Nrf, H, fft_len, snr)
clear
% 전송 데이터 스펙
fft_len = 64;
cp_len = fft_len / 4;
mod_type = 2;
data_len = fft_len * mod_type;
path = 7;
% 송신단 및 수신단에서의 총 데이터 스트림의 수
N_s = 2;
N_d = 2;
A = [1];  % BD 사용 시 수신기 당 안테나 수
% 송수신 각 배열 안테나의 안테나 수
ant_num = [8 1 0.5; 1 1 0.5];
N_tx = ant_num(1,1) * ant_num(1,2);
N_rx = ant_num(2,1) * ant_num(2,2);

% 채널 계수 생성
for d = 1:N_d
    [temp, ~ ] = FD_channel(fft_len + cp_len, path, ant_num);
    H(:,:,1+(d-1)*N_rx:d*N_rx,:) = temp;
end
Nrf = 4;
K = N_s;
SNR = 30;
%% beamforming start
SNR = 10^(SNR/10);
N = ant_num(1,1) * ant_num(1,2);
H = squeeze(H(:,1,:,:));
H = squeeze(mean(fft(H,fft_len,1),1));
% Vrf = randi([0 1], [N, Nrf])*2-1;
Vrf = orth(exp(-1j*2*pi*rand(N, Nrf)))*sqrt(N);
% Vrf = Vrf./abs(Vrf);
% Vrf = ones(N, Nrf);
P = eye(K);
iter = 1;
count = 0;
while( iter )% = 1 : 100
    temp = Vrf;
    for j = 1 : Nrf
        H_ = P^(-1/2) * H;
        Vj = Vrf; Vj(:,j) = [];
        A = H_ * (Vj*Vj') * H_';
        B = H_'/(A*A)*H_;
        D = H_'/A*H_;
        for i = 1 : N
            [zeta_B, zeta_D, eta_B, eta_D] = appendA(N,i,j,B,D,Vrf);
            c = (1 + zeta_D)*eta_B - zeta_B*eta_D;
            z = imag(2*conj(eta_B)*eta_D);
            if (real(c) >= 0)
                PI = asin(imag(c)/abs(c));
            else
                PI = pi - asin(imag(c)/abs(c));
            end
            theta_1 = - PI + asin(z/abs(c));
            theta_2 = pi - PI - asin(z/abs(c));
            Vrf(i,j) = appendB(theta_1, theta_2, zeta_B, zeta_D, eta_B, eta_D, A, N);
        end
    end
    if (sum(abs(temp - Vrf), 'all') < 0.5)
        iter = 0;
    end
    count= count + 1;
    if (rem(count,100) == 99)
%         iter = 0;
%         Vrf = randi([0 1], [N, Nrf])*2-1;
        Vrf = orth(exp(-1j*2*pi*rand(N, Nrf)));
        Vrf = Vrf./abs(Vrf);
    end    
end
Vd = Vrf'*H'/(H*(Vrf*Vrf')*H');
Q = Vd'*Vrf'*Vrf*Vd;