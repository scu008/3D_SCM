function [sym_hat, K, rx_H, W] = BD_precoding(sym, H, A)
% Block diagonalization�� ���� precoding matrix�� ����ϴ� �Լ�
% K: tr(W'W) - ����ȭ factor
% rx_H: precoding���� ���� ��ȯ ä�� matrix
% W: �ιݼ��Ŀ� ���� precoding matrix
% sym: �۽� ��ȣ
% H: MIMO ä�� matrix
% A: ���ű��� ���׳� ���� ex) 2���� ���űⰡ ���� 2���� ���׳��� ���� ��� - [2 2]
% rep: �Ƴ��α� �������� ���� �귻ġ�� ��

% ä�� ����� ����� ���Ѵ�.
[fft_len, Mr, Mt] = size(H);
% ���ű� ������ �ʱ�ȭ

if nargin < 3, A = ones(1,Mr); end

% �ݺ� ��� ����
rx_H = zeros(fft_len,Mr,Mt);
W = zeros(fft_len,Mt,Mt);
K = zeros(1,fft_len);
sym_hat = zeros(Mt,fft_len);


% ��� �����ɸ�� ���� precoding ����
for k = 1 : fft_len
    
    % �����ɸ��� �ϳ��� ä�� ����� ����
    t_H(:,:) = H(k,:,:);
    
    % ���� �ʱ�ȭ
    N = 1;  t_W = [];
    
    % Precoding matrix ���
    for i = 1:length(A)
        % SVD�� �̿��� Precoding matrix ���
        temp = t_H;
        temp(N : N + A(i) - 1,:) = [];
        [~,~,V] = svd(temp);
        temp = V(:,end-A(i)+1:end);
        
        % i��° ���ű⿡ ���� precoding matrix ��� �� ������ ����
        t_W = [t_W temp];
        N = N + A(i);
    end    
    
    % Precoding ����
    sym_hat(:,k) = t_W * sym(:,k);
    W(k,:,:) = t_W;
    K(k) = sqrt( trace(t_W * t_W') );
    
    % Precoding matrix�� ���� ��ȯ ä�� matrix ���
    rx_H(k,:,:) = t_H * t_W;
    
end



