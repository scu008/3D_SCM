function [ s_hat, W ] = MMSE_detect( Y, H, No )
    % MMSE detection�� �����ϴ� �Լ�

    % ä�� ����� ����� ���Ѵ�.
    [fft_len, Mr, Mt] = size(H);
    
    % �ݺ� ��� ����
    s_hat = zeros(Mt,fft_len);
    t_Y = [Y; zeros(Mt,fft_len)];
    W = zeros(fft_len,Mt,Mr);
    
    % ��� �����ɸ�� ���� ����
    for k = 1 : fft_len
        
        % �����ɸ��� �ϳ��� ä�� ����� ����
        t_H(:,:) = H(k,:,:);
        e_H = [t_H; sqrt(No) * eye(Mt)];
        
        % MMSE ����� ���
        G = inv( e_H' * e_H ) * e_H';
        
        % detect�� �ɺ� ����
        s_hat(:,k) = G * t_Y(:,k);
        W(k,:,:) = G;
        
    end

end

