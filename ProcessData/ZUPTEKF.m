function [x_h, q_hat]=ZUPTEKF(u,zupt)

N=length(u);
alpha=100;
R_acc =alpha*(0.009).^2*eye(3)';
R_gyr =alpha*(0.5*pi/180).^2*eye(3)'; % [rad/s]
Q=[R_acc, zeros(3); zeros(3), R_gyr];
P=zeros(9);
P(1:6,1:6)=alpha*(1e-10)*eye(6);
P(7:9,7:9)=alpha*(0.1*pi/180).^2*eye(3);
R=alpha*(0.01).^2*eye(4)';

% Allocate vecors
x_h=zeros(9,N);
quats=zeros(4,N);
q_hat=zeros(4,N);

q_hat_1=init_orient(u);
q_hat(:,1)=q_hat_1;

%Set gravity
g=[0; 0; 9.814];
Ts=1./200;

for k=2:N
    
    %Dynamic update
    u_h=u(:,k);
    q_hat(:, k)=exp_q_R(Ts*u_h(4:6),q_hat(:,k-1));
    Rb2t=quat2Rot(q_hat(:, k));
    f_t=Rb2t*u_h(1:3);
    acc_t=f_t+g;
    A=eye(6); A(1,4)=Ts; A(2,5)=Ts; A(3,6)=Ts;
    B=[(Ts^2)/2*eye(3);Ts*eye(3)];
    x_h(1:6,k)=A*x_h(1:6,k-1)+B*acc_t;
    St=[0 -f_t(3) f_t(2); f_t(3) 0 -f_t(1); -f_t(2) f_t(1) 0];
    O=zeros(3);
    I=eye(3);
    F=[I Ts*I (Ts^2)/2*St;
        O I Ts*St;
        O O I];
    G=[1/2*Ts^2*eye(3) O; Ts*eye(3) O; O -Ts*eye(3)];
    P=F*P*F'+G*Q*G';
    P=(P+P')/2;
    
    % Check if a zero velocity update should be done. If so, do the
    % following
    if zupt(k)==true
        
        K=(P(3:6,:)')/(P(3:6,3:6)'+R);
        z=-x_h(3:6,k);
        dx=K*z;
        x_h(1:6,k)=x_h(1:6,k)+dx(1:6);
        q_hat(:,k)=exp_q_L(-dx(7:9),q_hat(:,k));  
        quat=[q_hat(2:4,k); q_hat(1,k)];
        quats(:, k) = quat;
        P=P-K*P(3:6,:);
        P=(P+P')/2;
        
    end
    
end

end
