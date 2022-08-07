% 다중 반송파 fading test
clc, clear
close all


% 변수 설정
bit_len = 1000;
mod_type = 4;
freq = [11e9 11e9];
model = SCM();
model.n_path = 3;
model.n_mray = 2;
ray_ang = cell(1,3);

% 각도 설정
ray_ang{1} = rand(4,3)*(pi/2);
ray_ang{2} = rand(4,2)*(pi/2);
ray_ang{3} = rand(4,7)*(pi/2);
% model.full_ang = ray_ang;

% 데이터 생성
bit = randi([0 1], 1, bit_len);
sym = base_mod(bit, mod_type);


% 채널 생성
h = model.MC_channel(freq, size(sym,2));
y = model.MC_fading(sym, h);

tmp1(1:1,1:252) = y(1,:,:);
tmp2(1:1,1:252) = y(2,:,:);
sum(sum(abs(tmp1 - tmp2).^2))
