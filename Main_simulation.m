%% Full-digital Beamforming (ZF)
clc, clear

% 1: Capacity, 2: BER
mode = 1;

% ���� �Ķ���� ����
fft_len = 64;               % OFDM �ιݼ����� ��
mod_type = 2;               % ���� ���� ex) 1 - BPSK, 2 - QPSK, 4 - 16QAM, 6 - 64QAM, 8 - 256QAM
rx_node = 4;                % ���ű��� �� (���ű��� ���׳��� 1��)
tx_ant = 16;                % �������� ���׳� ��

% ���� �Ķ���� ����� ���α׷� ����
path = 1;
snr = -20:5:10;             % ���� ä�� SNR ����
iter = 100;                 % ���� �ݺ� Ƚ��
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
rx_ant = 1;                 % ���ű��� ���� ���׳� ��

% ���� �Ķ���� ����� ���α׷� ����
path = 1;
snr = -20:5:10;             % ���� ä�� SNR ����
iter = 300;                 % ���� �ݺ� Ƚ��
plot_format = 'r-*';
Hybrid_beamforming          % ���� ���α׷�
hold on
