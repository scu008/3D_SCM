function [ s_hat ] = QRM_detect( Y, H, mod_type, M, mode, No)
    % QRM detection을 수행하는 함수 ( mode: 1 - ZF, 2 - MMSE, 3 - sorted H, 4 - sorted MMSE)

    % 채널 행렬의 사이즈를 구한다.
    [fft_len, Mr, Mt] = size(H);
    point = 2^mod_type;
    
    % 최대 검사 횟수
    all = M * point;
    
    % 심볼 pool 생성
    sym_pool = base_mod( de2bi( 0 : point - 1, mod_type), mod_type).';    
    
    % M번 반복한 후보 심볼 pool
    M_sym_pool = ones(M,1) * sym_pool;
    M_sym_pool = M_sym_pool(:).';
    
    % 반복 행렬 선언
    t_Y = [Y; zeros(Mt, fft_len)];
    s_hat = zeros(Mt,fft_len);
    
    
    
    % 모든 서브케리어에 대해 수행
    for k = 1 : fft_len
        
        % 서브케리어 하나의 채널 행렬을 저장
        t_H(:,:) = H(k,:,:);
        
        if mode == 1
            
            % QR분해 수행
            [Q, R] = qr(t_H);
            % Z행렬 계산
            Z = Q' * Y(:,k);
            
        elseif mode == 2
            
            % noise 표준편차를 추가한 행렬 계산
            e_H = [t_H; sqrt(No) * eye(Mt)];
            
            % QR분해 수행
            [Q, R] = qr(e_H);
            % Z행렬 계산
            Z = Q' * t_Y(:,k);
            
            
        elseif mode == 3
            
            % 채널 메트릭스를 정렬
            G = inv(t_H' * t_H) * t_H';
            [~, G_idx] = sort( sum( ( G .* conj(G) ).' ), 2, 'descend' );
            t_H = t_H(:,G_idx);
            
            % QR분해 수행
            [Q, R] = qr(t_H);
            % Z행렬 계산
            Z = Q' * Y(:,k);
            
        elseif mode == 4
            
            % noise 표준편차를 추가한 행렬 계산
            e_H = [t_H; sqrt(No) * eye(Mt)];
            
            % 채널 메트릭스를 정렬
            G = inv(e_H' * e_H) * e_H';
            [~, G_idx] = sort( sum( ( G .* conj(G) ).' ), 2, 'descend' );
            e_H = e_H(:,G_idx);
            
            % QR분해 수행
            [Q, R] = qr(e_H);
            % Z행렬 계산
            Z = Q' * t_Y(:,k);
            
        end
        
        
        % 결과 및 거리 저장 행렬
        x = zeros(Mt, all);    dis = zeros(1, all);
        
        % 검사 신호 변수 초기화
        d = 1;
        
        % 검출 알고리즘
        for t = Mt : -1 : 1
            
            % 현재 후보 심볼 저장
            if d == 1
                x(t, 1 : point) = sym_pool;
            elseif d ~= M       
                temp_sym_pool = repmat( sym_pool, d, 1);
                x(t, 1 : d*point) = temp_sym_pool(:).';
            else
                x(t,:) = M_sym_pool;
            end
            
            % 거리 계산
            temp = ( Z(t) - R(t,:) * x(:,1:d*point) );
            dis(1:d*point) = dis(1:d*point) + ( real(temp).^2 + imag( temp ).^2 );            
            
            % 현재 결과 오름차순 정렬
            [dis(1 : d*point), idx] = sort(dis(1 : d*point));
            x(:, 1 : d*point) = x(:,idx);
            
            % 검사 변수 증가 및 제한
            d = d * point;
            if d > M, d = M; end
            
            % 선택 후 반복
            x(:,1 : d*point ) = repmat(x(:,1:d), 1, point);
            dis(1 : d*point ) = repmat(dis(1:d), 1, point);
            
        end
        
        % detect한 심볼 저장
        if mode == 1 || mode == 2
            s_hat(:,k) = x(:,1);
        elseif mode == 3 || mode == 4
            s_hat(G_idx,k) = x(:,1);
        end
        
    end

end