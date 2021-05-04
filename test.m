% 3차원 채널 행렬 고유값 확인
clc, clear

model = SCM();
model.n_path = 1;
model.n_mray = 1;
model.ant(1,16);

% 채널 계수 생성
[temp, rx_angle] = model.FD_channel(64);
h(1,:) = temp(1,1,:,:);

% 각도 탐색
[~, idx] = max( abs( temp(:,1,1,1) ) );
angle = rx_angle(:,idx);

% 탐색 벡터
N = 16;
offset = [-N * pi/32: pi/32 : 0   pi/32 : pi/32 : N*pi/32];
ang = angle(1:2) + offset;
% ang(1,:) = ang(1,:) * 0 + pi/2;
W = steer_precoding(model.fc, model.tx_ant, ang);

% 출력
res = abs(h * W).^2;
[~, idx] = max(res)