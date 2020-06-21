function W = beamformer(ant_num, rx_angle, Ns, mode, m_path)
% 송수신 각도를 통해 Beamforming matrix를 계산하는 함수
% ant_num: 안테나의 배열 [ 수평, 수직, 안테나 간격(lamda의 비율) ]
% rx_angle: 채널의 송수신 각도
% Ns: 데이터 스트림의 수
% mode: 빔포머 구조(1: All connected, 2: Sub connected)
% m_path: 각 스트림 전송마다 최대 전력으로 전송되는 경로

% 변수 초기화
if nargin < 5, m_path = ones(1,Ns); end

% 기본 변수 설정
%==========================================================
% 배열 변수
c = 3e8;                % 빛의 속도
f = 2.4e9;              % 캐리어 주파수
lamda = c/f;            % 신호의 파장
k = 2*pi / lamda;       % 파수
dy = ant_num(3)*lamda;          % 송신 안테나 간 거리 (가로)
dz = ant_num(3)*lamda;          % 송신 안테나 간 거리 (세로)
N = ant_num(1)*ant_num(2);

% 안테나 위치 행렬
temp1 = repmat(0:ant_num(1)-1, ant_num(2), 1);
temp2 = repmat(0:ant_num(2)-1, 1, ant_num(1));
ant_mat = [ zeros(N,1) temp1(:)*dy (temp2.')*dz];

% 좌표계 변환 함수
trans_f = @(t_phi, t_theta) [ sin(t_theta).*cos(t_phi); sin(t_theta).*sin(t_phi); cos(t_theta) ];


% 빔 matrix 계산 
if mode == 1
    % All connected array를 사용하는 경우
    
    W = zeros(N,Ns);
    for i = 1:Ns
        phi = rx_angle(1,m_path(i),i,i);    theta = rx_angle(2,m_path(i),i,i);
        W(:,i) = exp( -1j*k* ant_mat * trans_f(phi, theta) );
    end
    
elseif mode == 2
    % Sub connected array를 사용하는 경우
    
    W = zeros(N*Ns,Ns);
    for i = 1:Ns
        phi = rx_angle(1,m_path(i),i,i);    theta = rx_angle(2,m_path(i),i,i);
        W(1+(i-1)*N : i*N, i) = exp( -1j*k* ant_mat * trans_f(phi, theta) );
    end
end


end