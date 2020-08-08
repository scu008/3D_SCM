%% Full-digital Beamforming (ZF)
clc, clear

% 1: Capacity, 2: BER
mode = 1;

% ���� �Ķ���� ����
fft_len = 64;               % OFDM �ιݼ����� ��
mod_type = 2;               % ���� ���� ex) 1 - BPSK, 2 - QPSK, 4 - 16QAM, 6 - 64QAM, 8 - 256QAM
rx_node = 4;                % ���ű��� �� (���ű��� ���׳��� 1��)
tx_ant = 16;                % �������� ���׳� ��
snr = -20:5:10;               % ���� ä�� SNR ����
path = 7;
scatter = 10;
iter = 300;                 % ���� �ݺ� Ƚ��

% ���� �Ķ���� ����� ���α׷� ����
plot_format = 's-b';
FC_MIMO_testbed             % ���� ���α׷�
hold on

%% Beam steering method
clc, clear

% 1: Capacity, 2: BER
mode = 1;

% ���� �Ķ���� ����
fft_len = 64;               % OFDM �ιݼ����� ��
mod_type = 2;               % ���� ���� ex) 1 - BPSK, 2 - QPSK, 4 - 16QAM, 6 - 64QAM, 8 - 256QAM
rx_node = 4;                % ���ű��� �� (���ű��� ���׳��� 1��)
tx_ant = 16;                % �������� ���׳� ��
snr = -20:5:10;               % ���� ä�� SNR ����
path = 7;
scatter = 10;
iter = 300;                 % ���� �ݺ� Ƚ��

% ���� �Ķ���� ����� ���α׷� ����
rx_ant = 1;                 % ���ű��� ���� ���׳� ��
plot_format = 'r-*';
Hybrid_beamforming          % ���� ���α׷�
hold on

%% Sepideh beamforming
clc, clear

% 1: Capacity, 2: BER
mode = 1;

% ���� �Ķ���� ����
fft_len = 64;               % OFDM �ιݼ����� ��
mod_type = 2;               % ���� ���� ex) 1 - BPSK, 2 - QPSK, 4 - 16QAM, 6 - 64QAM, 8 - 256QAM
rx_node = 4;                % ���ű��� �� (���ű��� ���׳��� 1��)
tx_ant = 16;                % �������� ���׳� ��
snr = -20:5:10;               % ���� ä�� SNR ����
path = 7;
scatter = 10;
iter = 10;                 % ���� �ݺ� Ƚ��

% ���� �Ķ���� ����� ���α׷� ����
num_rf = 8;                 % �������� RF chain�� �� *(�ʼ� ����: rx_node <= num_rf <= tx_ant)*
rx_ant = 1;                 % ���ű��� ���� ���׳� ��
plot_format = 'm-*';
Hybrid_beamforming_sep      % ���� ���α׷�
hold on
