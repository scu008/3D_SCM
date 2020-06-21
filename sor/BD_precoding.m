function [sym_hat, K, rx_H, W] = BD_precoding(sym, H, A)
% Block diagonalization을 위한 precoding matrix를 계산하는 함수
% K: tr(W'W) - 정규화 factor
% rx_H: precoding으로 인한 변환 채널 matrix
% W: 부반송파에 따른 precoding matrix
% sym: 송신 신호
% H: MIMO 채널 matrix
% A: 수신기의 안테나 구성 ex) 2개의 수신기가 각각 2개의 안테나를 가진 경우 - [2 2]
% rep: 아날로그 빔포밍을 위한 브렌치의 수

% 채널 행렬의 사이즈를 구한다.
[fft_len, Mr, Mt] = size(H);
% 수신기 구성을 초기화

if nargin < 3, A = ones(1,Mr); end

% 반복 행렬 선언
rx_H = zeros(fft_len,Mr,Mt);
W = zeros(fft_len,Mt,Mt);
K = zeros(1,fft_len);
sym_hat = zeros(Mt,fft_len);


% 모든 서브케리어에 대해 precoding 수행
for k = 1 : fft_len
    
    % 서브케리어 하나의 채널 행렬을 저장
    t_H(:,:) = H(k,:,:);
    
    % 변수 초기화
    N = 1;  t_W = [];
    
    % Precoding matrix 계산
    for i = 1:length(A)
        % SVD를 이용한 Precoding matrix 계산
        temp = t_H;
        temp(N : N + A(i) - 1,:) = [];
        [~,~,V] = svd(temp);
        temp = V(:,end-A(i)+1:end);
        
        % i번째 수신기에 대한 precoding matrix 계산 후 참조점 갱신
        t_W = [t_W temp];
        N = N + A(i);
    end    
    
    % Precoding 수행
    sym_hat(:,k) = t_W * sym(:,k);
    W(k,:,:) = t_W;
    K(k) = sqrt( trace(t_W * t_W') );
    
    % Precoding matrix로 인한 변환 채널 matrix 계산
    rx_H(k,:,:) = t_H * t_W;
    
end



