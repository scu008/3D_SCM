function [ s_hat ] = QRM_detect( Y, H, mod_type, M, mode, No)
    % QRM detection�� �����ϴ� �Լ� ( mode: 1 - ZF, 2 - MMSE, 3 - sorted H, 4 - sorted MMSE)

    % ä�� ����� ����� ���Ѵ�.
    [fft_len, Mr, Mt] = size(H);
    point = 2^mod_type;
    
    % �ִ� �˻� Ƚ��
    all = M * point;
    
    % �ɺ� pool ����
    sym_pool = base_mod( de2bi( 0 : point - 1, mod_type), mod_type).';    
    
    % M�� �ݺ��� �ĺ� �ɺ� pool
    M_sym_pool = ones(M,1) * sym_pool;
    M_sym_pool = M_sym_pool(:).';
    
    % �ݺ� ��� ����
    t_Y = [Y; zeros(Mt, fft_len)];
    s_hat = zeros(Mt,fft_len);
    
    
    
    % ��� �����ɸ�� ���� ����
    for k = 1 : fft_len
        
        % �����ɸ��� �ϳ��� ä�� ����� ����
        t_H(:,:) = H(k,:,:);
        
        if mode == 1
            
            % QR���� ����
            [Q, R] = qr(t_H);
            % Z��� ���
            Z = Q' * Y(:,k);
            
        elseif mode == 2
            
            % noise ǥ�������� �߰��� ��� ���
            e_H = [t_H; sqrt(No) * eye(Mt)];
            
            % QR���� ����
            [Q, R] = qr(e_H);
            % Z��� ���
            Z = Q' * t_Y(:,k);
            
            
        elseif mode == 3
            
            % ä�� ��Ʈ������ ����
            G = inv(t_H' * t_H) * t_H';
            [~, G_idx] = sort( sum( ( G .* conj(G) ).' ), 2, 'descend' );
            t_H = t_H(:,G_idx);
            
            % QR���� ����
            [Q, R] = qr(t_H);
            % Z��� ���
            Z = Q' * Y(:,k);
            
        elseif mode == 4
            
            % noise ǥ�������� �߰��� ��� ���
            e_H = [t_H; sqrt(No) * eye(Mt)];
            
            % ä�� ��Ʈ������ ����
            G = inv(e_H' * e_H) * e_H';
            [~, G_idx] = sort( sum( ( G .* conj(G) ).' ), 2, 'descend' );
            e_H = e_H(:,G_idx);
            
            % QR���� ����
            [Q, R] = qr(e_H);
            % Z��� ���
            Z = Q' * t_Y(:,k);
            
        end
        
        
        % ��� �� �Ÿ� ���� ���
        x = zeros(Mt, all);    dis = zeros(1, all);
        
        % �˻� ��ȣ ���� �ʱ�ȭ
        d = 1;
        
        % ���� �˰���
        for t = Mt : -1 : 1
            
            % ���� �ĺ� �ɺ� ����
            if d == 1
                x(t, 1 : point) = sym_pool;
            elseif d ~= M       
                temp_sym_pool = repmat( sym_pool, d, 1);
                x(t, 1 : d*point) = temp_sym_pool(:).';
            else
                x(t,:) = M_sym_pool;
            end
            
            % �Ÿ� ���
            temp = ( Z(t) - R(t,:) * x(:,1:d*point) );
            dis(1:d*point) = dis(1:d*point) + ( real(temp).^2 + imag( temp ).^2 );            
            
            % ���� ��� �������� ����
            [dis(1 : d*point), idx] = sort(dis(1 : d*point));
            x(:, 1 : d*point) = x(:,idx);
            
            % �˻� ���� ���� �� ����
            d = d * point;
            if d > M, d = M; end
            
            % ���� �� �ݺ�
            x(:,1 : d*point ) = repmat(x(:,1:d), 1, point);
            dis(1 : d*point ) = repmat(dis(1:d), 1, point);
            
        end
        
        % detect�� �ɺ� ����
        if mode == 1 || mode == 2
            s_hat(:,k) = x(:,1);
        elseif mode == 3 || mode == 4
            s_hat(G_idx,k) = x(:,1);
        end
        
    end

end