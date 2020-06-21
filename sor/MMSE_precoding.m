function [ sym_hat, K, W ] = MMSE_precoding(sym, H, No)
% MMSE precoding�� �����ϴ� �Լ�
% H: ��� ä�� ���
% sym_hat: precoding�� ����� �۽� ��ȣ
% K: �� �ιݼ��ĸ��� ���� precoding ����ȭ ����
% No: ���� ����

% ä�� ����� ����� ���Ѵ�.
[fft_len, Mr, Mt] = size(H);

% �ݺ� ��� ����
K = zeros(1,fft_len);
W = zeros(fft_len,Mt,Mt);
sym_hat = zeros(Mt,fft_len);

% ��� �����ɸ�� ���� precoding ����
for k = 1 : fft_len
    
    % �����ɸ��� �ϳ��� ä�� ����� ����
    t_H(:,:) = H(k,:,:);

    % Precoding matrix ���
    G = t_H' * inv( t_H * t_H' + eye(Mr)*No );
    
    % Precoding ����
    sym_hat(:,k) = G * sym(:,k);
    W(k,:,:) = G;
    K(k) = sqrt( trace(G * G') );
end


