function [x_h, q_hat]=ZUPTEKFUWB(u,zupt,UWB,R_UWB,beacons)
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
q_hat=zeros(4,N);

q_hat_1=init_orient(u);
q_hat(:,1)=q_hat_1;

%Set gravity
g=[0; 0; 9.814];
Ts=1./200;

%Number of beacons
N_b=size(beacons,2);

for t=2:N
    
    u_h=u(:,t);
    
    % Update the navigation equations.
    q_hat(:, t)=exp_q_R(Ts*u_h(4:6),q_hat(:,t-1));
    Rb2t=quat2Rot(q_hat(:, t));
    f_t=Rb2t*u_h(1:3);
    acc_t=f_t+g;
    A=eye(6); A(1,4)=Ts; A(2,5)=Ts; A(3,6)=Ts;
    B=[(Ts^2)/2*eye(3);Ts*eye(3)];
    x_h(1:6,t)=A*x_h(1:6,t-1)+B*acc_t;
    St=[0 -f_t(3) f_t(2); f_t(3) 0 -f_t(1); -f_t(2) f_t(1) 0];
    O=zeros(3);
    I=eye(3);
    F=[I Ts*I O;
        O I Ts*St;
        O O I];
    G=[O O; Ts*quat2Rot(q_hat(:, t)) O; O -Ts*quat2Rot(q_hat(:, t))];
    P=F*P*F'+G*Q*G';
    P=(P+P')/2;
    
    % Check if a zero velocity update should be done. If so, do a
    % measurement update
    if zupt(t)==true
        
        K=(P(3:6,:)')/(P(3:6,3:6)'+R);
        z=-x_h(3:6,t);
        dx=K*z;
        x_h(1:6,t)=x_h(1:6,t)+dx(1:6);
        q_hat(:,t)=exp_q_L(-dx(7:9),q_hat(:,t));
        P=P-K*P(3:6,:);
        P=(P+P')/2;
        
    end
    
    %Check if a beacon-measurement is available. If so, do a measurement
    %update
    for i=1:N_b
        if (~isnan(UWB(i,t)))
            
            diffs=(x_h(1:3,t)-beacons(:,i));
            dists=sqrt(sum(diffs.^2));
            H=(diffs./dists)';
            K=((H*P(1:3,:))')/(H*P(1:3,1:3)*H'+R_UWB(i,i));
            z=UWB(i,t)-dists';
            dx=K*z;
            x_h(1:6,t)=x_h(1:6,t)+dx(1:6);
            q_hat(:,t)=exp_q_L(-dx(7:9),q_hat(:,t));
            P=P-K*H*P(1:3,:);
            P=(P+P')/2;
        end
    end
    
end

end
