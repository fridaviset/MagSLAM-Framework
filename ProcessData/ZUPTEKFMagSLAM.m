function [x_h, q_hat, m, P_m]=ZUPTEKFMagSLAM(u,zupt,sigma_y,Lambda,Indices,sigma_lin,N_m,xl,xu,yl,yu,zl,zu)

N=length(u);
alpha=100;
R_acc =alpha*(0.009).^2*eye(3)';
R_gyr =alpha*(0.5*pi/180).^2*eye(3)'; % [rad/s]
Q=[R_acc, zeros(3); zeros(3), R_gyr];
P=zeros(N_m+12,N_m+12);
P(1:6,1:6)=alpha*(1e-10)*eye(6);
P(7:9,7:9)=alpha*(0.1*pi/180).^2*eye(3);
R=alpha*(0.01).^2*eye(4)';
y_mag=u(7:9,:);

% Allocate vecors
x_h=zeros(9,N);
quats=zeros(4,N);
q_hat=zeros(4,N);
q_hat_1=init_orient(u); 
q_hat(:,1)=q_hat_1;

%Set gravity
g=[0; 0; 9.814];
Ts=1./200;

%Initialise the magnetic field map
P(10:end,10:end)=[sigma_lin^2*eye(3), zeros(3,N_m);
    zeros(N_m,3), Lambda];
m=zeros(N_m+3,1);

for k=2:N
    
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
    F=[I Ts*I O;
        O I Ts*St;
        O O I];
    G=[O O; Ts*quat2Rot(q_hat(:, k)) O; O -Ts*quat2Rot(q_hat(:, k))];
    P(1:9,1:9)=F*P(1:9,1:9)*F'+G*Q*G';
    P(1:9,1:9)=(P(1:9,1:9)+P(1:9,1:9)')/2;
    
    % Check if a zero velocity update should be done. If so, do the
    % following
    if zupt(k)==true
        
        K=(P(3:6,1:9)')/(P(3:6,3:6)'+R);
        z=-x_h(3:6,k);
        dx=K*z;
        x_h(1:6,k)=x_h(1:6,k)+dx(1:6);
        q_hat(:,k)=exp_q_L(-dx(7:9),q_hat(:,k));
        quat=[q_hat(2:4,k); q_hat(1,k)];
        quats(:, k) = quat;
        P(1:9,1:9)=P(1:9,1:9)-K*P(3:6,1:9);
        P(1:9,1:9)=(P(1:9,1:9)+P(1:9,1:9)')/2;
        P(1:9,1:9)=(P(1:9,1:9)+P(1:9,1:9)')/2;
        
    end
    
    if (~isnan(y_mag(1,k)))
        if mod(k,101)==0
            %KF meas update
            NablaPhi=[eye(3),Nabla_Phi3D(x_h(1:3,k),N_m,xl,xu,yl,yu,zl,zu,Indices)];
            f=NablaPhi*m;
            J=reshape(JacobianPhi3D(x_h(1:3,k),N_m,xl,xu,yl,yu,zl,zu,Indices),9,N_m)*m(4:end);
            J = reshape(J,3,3);
            
            %Prepare scew-symmetric magnetic field vector
            fx=[0 -f(3) f(2);
                f(3) 0 -f(1);
                -f(2) f(1) 0];
            
            %Meas state update full EKF
            H_t=[J, zeros(3), fx, NablaPhi];
            S_t=H_t*P*H_t'+sigma_y^2*eye(3);
            K_t=P*H_t'*inv(S_t);
            S_t=0.5*(S_t+S_t');
            eps_t=quat2Rot(q_hat(:,k))*y_mag(:,k)-(f);
            eta_t=K_t*eps_t;
            x_h(1:6,k)=x_h(1:6,k)+eta_t(1:6);
            q_hat(:,k)=exp_q_L(eta_t(7:9),q_hat(:,k));
            m=m+eta_t(10:end);
            
            P=P-K_t*S_t*K_t';
            P=1/2*(P+P');
            
        end
        
    end
    
end

P_m=P(10:end,10:end);

end
