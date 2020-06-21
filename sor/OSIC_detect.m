function [ s_hat ] = OSIC_detect( Y, H, mod_type, mode, No )
    % OSIC detection을 수행하는 함수
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
        
        
        % 제거해야할 열 번호를 저장
        idx = 1:Mt;
        
        % Post SNR이 가장 좋은 신호 부터 검출
        for kk = 1:Mt

            % pesudo invers 계산 및 선택 검출
            G = inv( e_H' * e_H ) * e_H';
            [~, i] = min( sum( abs(G).^2, 2) );
            s_hat(idx(i),k) = base_mod( base_demod( G(i,:) * t_Y(:,k), mod_type ), mod_type );
            
            % 채널 및 수신 신호 갱신
            t_Y(:,k) = t_Y(:,k) - e_H(:,i) * s_hat(idx(i),k);
            e_H(:,i) = [];  idx(i) = [];
            
        end
    end
end

