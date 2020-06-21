function [ s_hat, W ] = ZF_detect( Y, H )
    % ZF detection�� �����ϴ� �Լ�
    
    % ä�� ����� ����� ���Ѵ�.
    [fft_len, Mr, Mt] = size(H);
    
    % �ݺ� ��� ����
    s_hat = zeros(Mt,fft_len);
    W = zeros(fft_len,Mt,Mr);
    
    
    % ��� �����ɸ�� ���� ����
    for k = 1 : fft_len        
        
        % �����ɸ��� �ϳ��� ä�� ����� ����
        t_H(:,:) = H(k,:,:);        
        
        % �ǻ� ����� ��� 
        G = inv( t_H' * t_H ) * t_H';
                
        % detect�� �ɺ� ����
        s_hat(:,k) = G * Y(:,k);
        W(k,:,:) = G;
        
    end

end

