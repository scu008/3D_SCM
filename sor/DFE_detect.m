function [ s_hat ] = DFE_detect( Y, H, mod_type, mode, No )
    % DFE detection을 수행하는 함수
    % mode - 1: ZF, 2: MMSE
    
    if nargin <= 4, mode = 1; end
    
    % 채널 행렬의 사이즈를 구한다.
    [fft_len, Mr, Mt] = size(H);
    
    % 반복 행렬 선언
    s_hat = zeros(Mt,fft_len);
    
    % 수신 벡터 확장 여부 결정
    if mode == 2, t_Y = [Y; zeros(Mt,fft_len)];
    else t_Y = Y; end
    
    
    % 모든 서브케리어에 대해 수행
    for k = 1 : fft_len
        
        % 서브케리어 하나의 채널 행렬을 저장
        t_H(:,:) = H(k,:,:);
        
        % 채널 행렬 확장 여부 결정
        if mode == 2, e_H = [t_H; sqrt(No) * eye(Mt)];
        else e_H = t_H; end
        
        % pesudo invers 계산 및 정렬
        G = inv( e_H' * e_H ) * e_H';
        [~, idx] = sort( sum( abs(G).^2, 2), 'descend' );
        [Q, R] = qr( e_H(:,idx) );
        Z = Q' * t_Y(:,k);
        
        % Post SNR이 가장 좋은 신호 부터 검출
        for kk = 1:Mt
            
            % 검출 위치 
            i = Mt - kk + 1;

            % 이전 검출 신호 제거 후 검출 수행
            s_hat(i,k) = base_mod( base_demod( ( Z(i) - R(i,:) * s_hat(:,k) ) / R(i,i), mod_type ), mod_type );
            
        end
        
        s_hat(idx,k) = s_hat(:,k);
        
    end

end

