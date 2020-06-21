function W = beamformer(ant_num, rx_angle, Ns, mode, m_path)
% �ۼ��� ������ ���� Beamforming matrix�� ����ϴ� �Լ�
% ant_num: ���׳��� �迭 [ ����, ����, ���׳� ����(lamda�� ����) ]
% rx_angle: ä���� �ۼ��� ����
% Ns: ������ ��Ʈ���� ��
% mode: ������ ����(1: All connected, 2: Sub connected)
% m_path: �� ��Ʈ�� ���۸��� �ִ� �������� ���۵Ǵ� ���

% ���� �ʱ�ȭ
if nargin < 5, m_path = ones(1,Ns); end

% �⺻ ���� ����
%==========================================================
% �迭 ����
c = 3e8;                % ���� �ӵ�
f = 2.4e9;              % ĳ���� ���ļ�
lamda = c/f;            % ��ȣ�� ����
k = 2*pi / lamda;       % �ļ�
dy = ant_num(3)*lamda;          % �۽� ���׳� �� �Ÿ� (����)
dz = ant_num(3)*lamda;          % �۽� ���׳� �� �Ÿ� (����)
N = ant_num(1)*ant_num(2);

% ���׳� ��ġ ���
temp1 = repmat(0:ant_num(1)-1, ant_num(2), 1);
temp2 = repmat(0:ant_num(2)-1, 1, ant_num(1));
ant_mat = [ zeros(N,1) temp1(:)*dy (temp2.')*dz];

% ��ǥ�� ��ȯ �Լ�
trans_f = @(t_phi, t_theta) [ sin(t_theta).*cos(t_phi); sin(t_theta).*sin(t_phi); cos(t_theta) ];


% �� matrix ��� 
if mode == 1
    % All connected array�� ����ϴ� ���
    
    W = zeros(N,Ns);
    for i = 1:Ns
        phi = rx_angle(1,m_path(i),i,i);    theta = rx_angle(2,m_path(i),i,i);
        W(:,i) = exp( -1j*k* ant_mat * trans_f(phi, theta) );
    end
    
elseif mode == 2
    % Sub connected array�� ����ϴ� ���
    
    W = zeros(N*Ns,Ns);
    for i = 1:Ns
        phi = rx_angle(1,m_path(i),i,i);    theta = rx_angle(2,m_path(i),i,i);
        W(1+(i-1)*N : i*N, i) = exp( -1j*k* ant_mat * trans_f(phi, theta) );
    end
end


end