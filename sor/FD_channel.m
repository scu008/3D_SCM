% 2018 02
function [coeff, result_angle] = FD_channel(sym_len, PDP, ant_num, normalized_power, distance_rate, t_angle, fs, dopp_f, mode, rad_pattern)
% 3차원 MIMO 채널의 계수를 생성하는 함수

% sym_len - 송신 심볼의 수
% PDP [W] - 스칼라 값일 경우 multipath로 인식, 벡터의 경우 Power delay profile로 사용
% ant_num - 송수신기의 안테나 배열 [ Tx_수평, Tx_수직, 송신 안테나 간격(lamda의 비율) ; Rx_수평, Rx_수직, 수신 안테나 간격(lamda의 비율) ]
% normalized_power - 정규화된 송신 전력 ('1'의 상대적인 비율)
% distance_rate - 정규화된 거리 ('1'의 상대적인 비율)
% t_angle - 1~2 행에는 탭에 따른 AoD가 3~4 행에는 탭에 따른 AoA를 저장
% fs - sampling frequency
% dopp_f - 1 by n or 2 by n 도플러 주파수 행렬 (첫 행의 엘리먼트는 탭의 도플러 주파수, 둘 째 행의 엘리먼트는 탭의 도플러 스펙트럼 )
% dopp_f의 경우 도플러 주파수만 존재하는 벡터로 사용 가능 ( default spectrum - Jake model)
% dopp_f에 스칼라 값을 넣을 시 모든 탭에 동일한 도플러 효과 적용
% mode - 1: NLOS, mode - 2: LOS
% coeff - 생성된 채널계수가 저장됨 ( 4차원: 채널 탭수, 송신 심볼 길이, 수신 안테나 index, 수신안테나 index )
% result_angle - 입력 값 혹은 임의로 생성된 AoD, AoA를 저장

% *********** dopp_filter.m, Doppler_rnd.m 필요 *************** 
% ( 시변 채널 설정을 위한 기본 변수 )


% 매개변수에 따른 초기화
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


% 기본 변수 설정
%==========================================================
% 배열 팩터 변수
c = 3e8;                % 빛의 속도
f = 2.4e9;              % 캐리어 주파수
lamda = c/f;            % 신호의 파장
k = 2*pi / lamda;       % 파수
tdy = ant_num(1,3)*lamda;          % 송신 안테나 간 거리 (가로)
tdz = ant_num(1,3)*lamda;          % 송신 안테나 간 거리 (세로)
rdy = ant_num(2,3)*lamda;          % 수신 안테나 간 거리 (가로)
rdz = ant_num(2,3)*lamda;          % 수신 안테나 간 거리 (세로)

% 송수신 안테나 수
Ntx = ant_num(1,1) * ant_num(1,2);
Nrx = ant_num(2,1) * ant_num(2,2);

% 안테나 위치 행렬
temp1 = repmat(0:ant_num(1,1)-1, ant_num(1,2), 1);
temp2 = repmat(0:ant_num(1,2)-1, 1, ant_num(1,1));
ant_mat1 = [ zeros(Ntx,1) temp1(:)*tdy (temp2.')*tdz];

temp3 = repmat(0:ant_num(2,1)-1, ant_num(2,2), 1);
temp4 = repmat(0:ant_num(2,2)-1, 1, ant_num(2,1));
ant_mat2 = [ zeros(Nrx,1) temp3(:)*rdy (temp4.')*rdz];

% 좌표계 변환 함수
trans_f = @(t_phi, t_theta) [ sin(t_theta).*cos(t_phi); sin(t_theta).*sin(t_phi); cos(t_theta) ];


% Path loss 계산
%==========================================================
% 일반적 3 ~ 4, 이상적 2
exp_beta = 3;

% 일반적 5 ~ 12, LTE: 10
sigma = 0;
shadowing = randn(1) * sigma;

% + log10(k) : 송수신 안테나 gain 및 중간 회로에 의한 손실 = 1 로 정규화
path_loss = -10 * exp_beta * log10(distance_rate);

% path loss 적용
loss = 10^( (path_loss + shadowing) / 10 );


% 채널 탭에 PDP 및 AoD, AoA 적용 
%==========================================================
% PDP 저장
if length(PDP) == 1
    % 기본 채널 프로파일
    Multipath = PDP;
    tap = exp( -(1:Multipath) / 5 );
    tap = sqrt( tap / sum(tap) );
else
    tap = sqrt( PDP / sum(PDP) );
end

% 다중 경로 및 angluar spread (AS) 정의
tap_len = length( tap );


% AoD 저장 (phi, theta)
result_angle = zeros(4,tap_len);

if sum( sum( t_angle(1:2,:), 1 ), 2 ) == 0
    result_angle(1,:) = -pi/2 + rand(1,tap_len)*pi;
    result_angle(2,:) = rand(1,tap_len)*pi;
else
    % AoD, AoA 행렬이 채널 탭 수와 일치 하지 않을 경우 처리
    [~,angle_num] = size(t_angle);
    
    if angle_num < tap_len
        t_angle = [t_angle repmat( t_angle(end), 1, tap_len-angle_num )];
    elseif angle_num > tap_len
        t_angle = t_angle(:,tap_len);
    end

    result_angle(1:2,:) = t_angle(1:2,:);
end

% AoA 저장 (phi, theta)
if sum( sum( t_angle(3:4,:), 1 ), 2 ) == 0
    result_angle(3,:) = -pi/2 + rand(1,tap_len)*pi;
    result_angle(4,:) = rand(1,tap_len)*pi;
else
    result_angle(3:4,:) = t_angle(3:4,:);
end


% tap 생성
tap = tap.' * ones(1, sym_len);
tap = tap * sqrt(normalized_power) * sqrt(loss);

% 채널 계수 생성
%===============================================================

% 변환된 송수신 위치 행렬 (채널 탭 수 by 안테나 수)
if rad_pattern == 1, rad_pattern = @(x,y) ones(tap_len,1); end
temp_tx_v = exp( 1j*k* ant_mat1 * trans_f( result_angle(1,:), result_angle(2,:) ) ).';
tx_v = temp_tx_v .* repmat( rad_pattern(result_angle(1,:).', result_angle(2,:).'), 1, Ntx);

temp_rx_v = exp( 1j*k* ant_mat2 * trans_f( result_angle(3,:), result_angle(4,:) ) ).';
rx_v = temp_rx_v .* repmat( rad_pattern(result_angle(3,:).', result_angle(4,:).'), 1, Nrx);


% 랜덤 시퀀스 생성
if (ant_num(1,3) > 0.5) && (ant_num(2,3) > 0.5), corr = 0;
elseif ant_num(1,3) > 0.5, corr = Ntx;
elseif ant_num(2,3) > 0.5, corr = Nrx;
else corr = 1; 
end

temp_rnd = zeros(tap_len, sym_len, corr);
for i = 1:corr, temp_rnd(:,:,i) = Doppler_rnd(sym_len, tap_len, fs, dopp_f, mode); end


% 채널 계수 배열 생성
coeff = zeros(tap_len, sym_len, Nrx, Ntx);

% 채널 계수 저장
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