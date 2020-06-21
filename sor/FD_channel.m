% 2018 02
function [coeff, result_angle] = FD_channel(sym_len, PDP, ant_num, normalized_power, distance_rate, t_angle, fs, dopp_f, mode, rad_pattern)
% 3���� MIMO ä���� ����� �����ϴ� �Լ�

% sym_len - �۽� �ɺ��� ��
% PDP [W] - ��Į�� ���� ��� multipath�� �ν�, ������ ��� Power delay profile�� ���
% ant_num - �ۼ��ű��� ���׳� �迭 [ Tx_����, Tx_����, �۽� ���׳� ����(lamda�� ����) ; Rx_����, Rx_����, ���� ���׳� ����(lamda�� ����) ]
% normalized_power - ����ȭ�� �۽� ���� ('1'�� ������� ����)
% distance_rate - ����ȭ�� �Ÿ� ('1'�� ������� ����)
% t_angle - 1~2 �࿡�� �ǿ� ���� AoD�� 3~4 �࿡�� �ǿ� ���� AoA�� ����
% fs - sampling frequency
% dopp_f - 1 by n or 2 by n ���÷� ���ļ� ��� (ù ���� ������Ʈ�� ���� ���÷� ���ļ�, �� ° ���� ������Ʈ�� ���� ���÷� ����Ʈ�� )
% dopp_f�� ��� ���÷� ���ļ��� �����ϴ� ���ͷ� ��� ���� ( default spectrum - Jake model)
% dopp_f�� ��Į�� ���� ���� �� ��� �ǿ� ������ ���÷� ȿ�� ����
% mode - 1: NLOS, mode - 2: LOS
% coeff - ������ ä�ΰ���� ����� ( 4����: ä�� �Ǽ�, �۽� �ɺ� ����, ���� ���׳� index, ���ž��׳� index )
% result_angle - �Է� �� Ȥ�� ���Ƿ� ������ AoD, AoA�� ����

% *********** dopp_filter.m, Doppler_rnd.m �ʿ� *************** 
% ( �ú� ä�� ������ ���� �⺻ ���� )


% �Ű������� ���� �ʱ�ȭ
%=========================================================
if nargin < 10, rad_pattern = 1; end
if nargin < 9, mode = 1; end
if nargin < 8, dopp_f = 0; end
if nargin < 7, fs = 10e6; end
if nargin < 6 || length(t_angle) < 4, t_angle = zeros(4,1); end
if nargin < 5, distance_rate = 1; end
if nargin < 4, normalized_power = 1; end
if nargin < 3 || length(ant_num) < 3, ant_num = [1 1 2 2; 1 1 2 2]; end
if nargin < 2, PDP = 7; end


% �⺻ ���� ����
%==========================================================
% �迭 ���� ����
c = 3e8;                % ���� �ӵ�
f = 2.4e9;              % ĳ���� ���ļ�
lamda = c/f;            % ��ȣ�� ����
k = 2*pi / lamda;       % �ļ�
tdy = ant_num(1,3)*lamda;          % �۽� ���׳� �� �Ÿ� (����)
tdz = ant_num(1,3)*lamda;          % �۽� ���׳� �� �Ÿ� (����)
rdy = ant_num(2,3)*lamda;          % ���� ���׳� �� �Ÿ� (����)
rdz = ant_num(2,3)*lamda;          % ���� ���׳� �� �Ÿ� (����)

% �ۼ��� ���׳� ��
Ntx = ant_num(1,1) * ant_num(1,2);
Nrx = ant_num(2,1) * ant_num(2,2);

% ���׳� ��ġ ���
temp1 = repmat(0:ant_num(1,1)-1, ant_num(1,2), 1);
temp2 = repmat(0:ant_num(1,2)-1, 1, ant_num(1,1));
ant_mat1 = [ zeros(Ntx,1) temp1(:)*tdy (temp2.')*tdz];

temp3 = repmat(0:ant_num(2,1)-1, ant_num(2,2), 1);
temp4 = repmat(0:ant_num(2,2)-1, 1, ant_num(2,1));
ant_mat2 = [ zeros(Nrx,1) temp3(:)*rdy (temp4.')*rdz];

% ��ǥ�� ��ȯ �Լ�
trans_f = @(t_phi, t_theta) [ sin(t_theta).*cos(t_phi); sin(t_theta).*sin(t_phi); cos(t_theta) ];


% Path loss ���
%==========================================================
% �Ϲ��� 3 ~ 4, �̻��� 2
exp_beta = 3;

% �Ϲ��� 5 ~ 12, LTE: 10
sigma = 0;
shadowing = randn(1) * sigma;

% + log10(k) : �ۼ��� ���׳� gain �� �߰� ȸ�ο� ���� �ս� = 1 �� ����ȭ
path_loss = -10 * exp_beta * log10(distance_rate);

% path loss ����
loss = 10^( (path_loss + shadowing) / 10 );


% ä�� �ǿ� PDP �� AoD, AoA ���� 
%==========================================================
% PDP ����
if length(PDP) == 1
    % �⺻ ä�� ��������
    Multipath = PDP;
    tap = exp( -(1:Multipath) / 5 );
    tap = sqrt( tap / sum(tap) );
else
    tap = sqrt( PDP / sum(PDP) );
end

% ���� ��� �� angluar spread (AS) ����
tap_len = length( tap );


% AoD ���� (phi, theta)
result_angle = zeros(4,tap_len);

if sum( sum( t_angle(1:2,:), 1 ), 2 ) == 0
    result_angle(1,:) = -pi/2 + rand(1,tap_len)*pi;
    result_angle(2,:) = rand(1,tap_len)*pi;
else
    % AoD, AoA ����� ä�� �� ���� ��ġ ���� ���� ��� ó��
    [~,angle_num] = size(t_angle);
    
    if angle_num < tap_len
        t_angle = [t_angle repmat( t_angle(end), 1, tap_len-angle_num )];
    elseif angle_num > tap_len
        t_angle = t_angle(:,tap_len);
    end

    result_angle(1:2,:) = t_angle(1:2,:);
end

% AoA ���� (phi, theta)
if sum( sum( t_angle(3:4,:), 1 ), 2 ) == 0
    result_angle(3,:) = -pi/2 + rand(1,tap_len)*pi;
    result_angle(4,:) = rand(1,tap_len)*pi;
else
    result_angle(3:4,:) = t_angle(3:4,:);
end


% tap ����
tap = tap.' * ones(1, sym_len);
tap = tap * sqrt(normalized_power) * sqrt(loss);

% ä�� ��� ����
%===============================================================

% ��ȯ�� �ۼ��� ��ġ ��� (ä�� �� �� by ���׳� ��)
if rad_pattern == 1, rad_pattern = @(x,y) ones(tap_len,1); end
temp_tx_v = exp( 1j*k* ant_mat1 * trans_f( result_angle(1,:), result_angle(2,:) ) ).';
tx_v = temp_tx_v .* repmat( rad_pattern(result_angle(1,:).', result_angle(2,:).'), 1, Ntx);

temp_rx_v = exp( 1j*k* ant_mat2 * trans_f( result_angle(3,:), result_angle(4,:) ) ).';
rx_v = temp_rx_v .* repmat( rad_pattern(result_angle(3,:).', result_angle(4,:).'), 1, Nrx);


% ���� ������ ����
if (ant_num(1,3) > 0.5) && (ant_num(2,3) > 0.5), corr = 0;
elseif ant_num(1,3) > 0.5, corr = Ntx;
elseif ant_num(2,3) > 0.5, corr = Nrx;
else corr = 1; 
end

temp_rnd = zeros(tap_len, sym_len, corr);
for i = 1:corr, temp_rnd(:,:,i) = Doppler_rnd(sym_len, tap_len, fs, dopp_f, mode); end


% ä�� ��� �迭 ����
coeff = zeros(tap_len, sym_len, Nrx, Ntx);

% ä�� ��� ����
for i = 1:Nrx
    for j = 1:Ntx
        temp = tx_v(:,j) .* rx_v(:,i) * ones(1,sym_len);
        
        if corr == 0, rnd = Doppler_rnd(sym_len, tap_len, fs, dopp_f, mode);
        elseif corr == Ntx, rnd = temp_rnd(:,:,j);
        elseif corr == Nrx, rnd = temp_rnd(:,:,i);
        else rnd = temp_rnd;
        end
        
        coeff(:,:,i,j) = tap .* rnd .* temp;
    end
end


end