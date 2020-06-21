function [ s_hat ] = DFE_detect( Y, H, mod_type, mode, No )
    % DFE detection�� �����ϴ� �Լ�
    % mode - 1: ZF, 2: MMSE
    
    if nargin <= 4, mode = 1; end
    
    % ä�� ����� ����� ���Ѵ�.
    [fft_len, Mr, Mt] = size(H);
    
    % �ݺ� ��� ����
    s_hat = zeros(Mt,fft_len);
    
    % ���� ���� Ȯ�� ���� ����
    if mode == 2, t_Y = [Y; zeros(Mt,fft_len)];
    else t_Y = Y; end
    
    
    % ��� �����ɸ�� ���� ����
    for k = 1 : fft_len
        
        % �����ɸ��� �ϳ��� ä�� ����� ����
        t_H(:,:) = H(k,:,:);
        
        % ä�� ��� Ȯ�� ���� ����
        if mode == 2, e_H = [t_H; sqrt(No) * eye(Mt)];
        else e_H = t_H; end
        
        % pesudo invers ��� �� ����
        G = inv( e_H' * e_H ) * e_H';
        [~, idx] = sort( sum( abs(G).^2, 2), 'descend' );
        [Q, R] = qr( e_H(:,idx) );
        Z = Q' * t_Y(:,k);
        
        % Post SNR�� ���� ���� ��ȣ ���� ����
        for kk = 1:Mt
            
            % ���� ��ġ 
            i = Mt - kk + 1;

            % ���� ���� ��ȣ ���� �� ���� ����
            s_hat(i,k) = base_mod( base_demod( ( Z(i) - R(i,:) * s_hat(:,k) ) / R(i,i), mod_type ), mod_type );
            
        end
        
        s_hat(idx,k) = s_hat(:,k);
        
    end

end

