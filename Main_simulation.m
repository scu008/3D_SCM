%% Full-digital Beamforming (ZF)
clc, clear

% 1: Capacity, 2: BER
mode = 1;

% 전송 파라미터 설정
fft_len = 64;               % OFDM 부반송파의 수
mod_type = 2;               % 변조 차수 ex) 1 - BPSK, 2 - QPSK, 4 - 16QAM, 6 - 64QAM, 8 - 256QAM
rx_node = 4;                % 수신기의 수 (수신기의 안테나는 1개)
tx_ant = 16;                % 기지국의 안테나 수
snr = -20:5:10;               % 전송 채널 SNR 범위
path = 7;
scatter = 10;
iter = 300;                 % 전송 반복 횟수

% 설정 파라미터 기반의 프로그램 실행
plot_format = 's-b';
FC_MIMO_testbed             % 실행 프로그램
hold on

%% Beam steering method
clc, clear

% 1: Capacity, 2: BER
mode = 1;

% 전송 파라미터 설정
fft_len = 64;               % OFDM 부반송파의 수
mod_type = 2;               % 변조 차수 ex) 1 - BPSK, 2 - QPSK, 4 - 16QAM, 6 - 64QAM, 8 - 256QAM
rx_node = 4;                % 수신기의 수 (수신기의 안테나는 1개)
tx_ant = 16;                % 기지국의 안테나 수
snr = -20:5:10;               % 전송 채널 SNR 범위
path = 7;
scatter = 10;
iter = 300;                 % 전송 반복 횟수

% 설정 파라미터 기반의 프로그램 실행
rx_ant = 1;                 % 수신기의 개별 안테나 수
plot_format = 'r-*';
Hybrid_beamforming          % 실행 프로그램
hold on

%% Sepideh beamforming
clc, clear

% 1: Capacity, 2: BER
mode = 1;

% 전송 파라미터 설정
fft_len = 64;               % OFDM 부반송파의 수
mod_type = 2;               % 변조 차수 ex) 1 - BPSK, 2 - QPSK, 4 - 16QAM, 6 - 64QAM, 8 - 256QAM
rx_node = 4;                % 수신기의 수 (수신기의 안테나는 1개)
tx_ant = 16;                % 기지국의 안테나 수
snr = -20:5:10;               % 전송 채널 SNR 범위
path = 7;
scatter = 10;
iter = 10;                 % 전송 반복 횟수

% 설정 파라미터 기반의 프로그램 실행
num_rf = 8;                 % 기지국의 RF chain의 수 *(필수 조건: rx_node <= num_rf <= tx_ant)*
rx_ant = 1;                 % 수신기의 개별 안테나 수
plot_format = 'm-*';
Hybrid_beamforming_sep      % 실행 프로그램
hold on
