function [ s_hat ] = OSIC_detect( Y, H, mod_type, mode, No )
    % OSIC detection�� �����ϴ� �Լ�
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
        
        
        % �����ؾ��� �� ��ȣ�� ����
        idx = 1:Mt;
        
        % Post SNR�� ���� ���� ��ȣ ���� ����
        for kk = 1:Mt

            % pesudo invers ��� �� ���� ����
            G = inv( e_H' * e_H ) * e_H';
            [~, i] = min( sum( abs(G).^2, 2) );
            s_hat(idx(i),k) = base_mod( base_demod( G(i,:) * t_Y(:,k), mod_type ), mod_type );
            
            % ä�� �� ���� ��ȣ ����
            t_Y(:,k) = t_Y(:,k) - e_H(:,i) * s_hat(idx(i),k);
            e_H(:,i) = [];  idx(i) = [];
            
        end
    end
end

