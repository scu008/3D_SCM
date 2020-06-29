% 3차원 채널 행렬 고유값 확인
clc, clear

% 채널 모델 생성
model = SCM();
model.ant(4,4);


tmp = model.FD_channel(10);
h(:,:) = tmp(1,1,:,:);