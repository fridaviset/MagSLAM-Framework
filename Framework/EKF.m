function [p_hat,v_hat,q_hat,m,P]=EKF(N,T,u,y_mag,q_0,p_0,v_0,zupt,sigma_y,Lambda,Indices,sigma_lin,N_m,xl,xu,yl,yu,zl,zu)
%Copyright (C) 2022 by Frida Viset

%Pre-allocate position trajectory for all particles
q_hat=zeros(4,N);
v_hat=zeros(3,N);
p_hat=zeros(3,N);

%Pre-allocate space for the full state covariance
P=zeros(N_m+12,N_m+12);
P(1:9,1:9)=0.001*eye(9);

%Initialise orientation, velocity and position estimates
q_hat(:,1)=q_0;
v_hat(:,1)=v_0;
p_hat(:,1)=p_0;

%Initialise the magnetic field map
P(10:end,10:end)=[sigma_lin^2*eye(3), zeros(3,N_m);
    zeros(N_m,3), Lambda];
m=zeros(N_m+3,1);

%Set gravity
g=[0; 0; 9.7941];

% IMU process noise
R_acc =(0.009).^2*eye(3)';
R_gyr =(0.5*pi/180).^2*eye(3)'; % [rad/s]
R_vel = (0.01).^2*eye(3)';

%iterate through all timesteps
for t=2:N
    
    % Transform measured force to force in
    % the navigation coordinate system.
    f_t=quat2Rot(q_hat(:,t-1))*u(1:3,t-1);

    % Create a ske symmetric matrix of the specific fore vector
    St=[0 -f_t(3) f_t(2); f_t(3) 0 -f_t(1); -f_t(2) f_t(1) 0];
    
    %KF Dyn update
    O=zeros(3); I=eye(3);
    F=[O I O;
        O O St;
        O O O];
    
    % Noise gain matrix
    G=[O O; quat2Rot(q_hat(:,t-1)) O; O -quat2Rot(q_hat(:,t-1))];
    
    delta_acc=f_t-g;
    v_hat(:,t)=v_hat(:,t-1)+T*delta_acc;
    p_hat(:,t)=p_hat(:,t-1)+T*v_hat(:,t-1);
    delta_q=T*u(4:6,t-1);
    q_hat(:,t)=exp_q_L(delta_q,q_hat(:,t-1));
    
    P(1:9,1:9)=F*P(1:9,1:9)*F'+G*[R_acc, zeros(3); zeros(3), R_gyr]*G';
    
    %Perform ZUPT-meas-update
    if zupt(t)
    
        %Meas state update full EKF
        S_t=P(4:6,4:6)+R_vel;
        K_t=P(:,4:6)*inv(S_t);
        S_t=0.5*(S_t+S_t');
        eps_t=-v_hat(:,t);
        eta_t=K_t*eps_t;
        p_hat(:,t)=p_hat(:,t)+eta_t(1:3);
        v_hat(:,t)=v_hat(:,t)+eta_t(4:6);
        q_hat(:,t)=exp_q_L(eta_t(7:9),q_hat(:,t));
        m=m+eta_t(10:end);
        
        P=P-K_t*S_t*K_t';
        P=1/2*(P+P');
    
    if (isnan(y_mag(1,t)))
        %Skip meas update
    elseif 1==0
        
        %KF meas update
        NablaPhi=[eye(3),Nabla_Phi3D(p_hat(:,t),N_m,xl,xu,yl,yu,zl,zu,Indices)];
        f=NablaPhi*m;
        J=reshape(JacobianPhi3D(p_hat(:,t),N_m,xl,xu,yl,yu,zl,zu,Indices),9,N_m)*m(4:end);
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
        eps_t=quat2Rot(q_hat(:,t))*y_mag(:,t)-(f);
        eta_t=K_t*eps_t;
        p_hat(:,t)=p_hat(:,t)+eta_t(1:3);
        v_hat(:,t)=p_hat(:,t)+eta_t(4:6);
        q_hat(:,t)=exp_q_L(eta_t(7:9),q_hat(:,t));
        m=m+eta_t(10:end);
        
        P=P-K_t*S_t*K_t';
        P=1/2*(P+P');
        
    end
    
    
    
end





end