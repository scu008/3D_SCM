function [ s_hat ] = ML_detect( Y, H, mod_type )
    % ML detection�� �����ϴ� �Լ�
    
    % ä�� ����� ����� ���Ѵ�.
    [fft_len, Mr, Mt] = size(H);
    all = 2^(Mt * mod_type);
    
    % �ɺ� pool ����
    sym_pool = base_mod( de2bi( 0 : all - 1, mod_type * Mt), mod_type).';
    
    % �ݺ� ��� ����
    s_hat = zeros(Mt,fft_len);
    
    
    % ��� �����ɸ�� ���� ����
    for k = 1 : fft_len
        
        % �����ɸ��� �ϳ��� ä�� ����� ����
        t_H(:,:) = H(k,:,:);
        
        % �Ÿ� ���
        t_dis = repmat( Y(:,k), 1, all ) - (t_H * sym_pool);
        t_dis = sum( real(t_dis).^2 + imag(t_dis).^2 );
        [~, idx] = min(t_dis);
        
        % detect�� �ɺ� ����
        s_hat(:,k) = sym_pool(:,idx);
        
    end
end
