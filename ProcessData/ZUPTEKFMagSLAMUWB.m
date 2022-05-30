function [x_h, q_hat, m, P_m]=ZUPTEKFMagSLAMUWB(u,zupt,sigma_y,Lambda,Indices,sigma_lin,N_m,xl,xu,yl,yu,zl,zu,UWB,R_UWB,beacons)
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

% Initialize the navigation state vector x_h, and the quaternion vector
% quat.
q_hat_1=init_orient(u); % Subfunction located further down in the file.
q_hat(:,1)=q_hat_1;

%Set gravity
g=[0; 0; 9.814];
Ts=1./200;

%Initialise the magnetic field map
P(10:end,10:end)=[sigma_lin^2*eye(3), zeros(3,N_m);
    zeros(N_m,3), Lambda];
m=zeros(N_m+3,1);

%Number of beacons
N_b=size(beacons,2);

for t=2:N
    
    u_h=u(:,t);
    
    % Update the navigation equations.
    q_hat(:, t)=exp_q_R(Ts*u_h(4:6),q_hat(:,t-1));
    
    % Convert quaternion to a rotation matrix
    Rb2t=quat2Rot(q_hat(:, t));
    
    % Transform measured force to force in
    % the navigation coordinate system.
    f_t=Rb2t*u_h(1:3);
    
    % Subtract (add) the gravity, to obtain accelerations in navigation
    % coordinat system.
    acc_t=f_t+g;
    
    % State space model matrices
    A=eye(6); A(1,4)=Ts; A(2,5)=Ts; A(3,6)=Ts;
    B=[(Ts^2)/2*eye(3);Ts*eye(3)];
    
    % Update the position and velocity estimates.
    x_h(1:6,t)=A*x_h(1:6,t-1)+B*acc_t;
    
    % Create a skew symmetric matrix of the specific force vector
    St=[0 -f_t(3) f_t(2); f_t(3) 0 -f_t(1); -f_t(2) f_t(1) 0];
    
    % Zero matrix
    O=zeros(3);
    
    % Identity matrix
    I=eye(3);
    
    F=[I Ts*I O;
        O I Ts*St;
        O O I];
    
    % Noise gain matrix
    G=[O O; Ts*quat2Rot(q_hat(:, t)) O; O -Ts*quat2Rot(q_hat(:, t))];
    
    % Update the filter state covariance matrix P.
    P(1:9,1:9)=F*P(1:9,1:9)*F'+G*Q*G';
    
    % Make sure the filter state covariance matrix is symmetric.
    P(1:9,1:9)=(P(1:9,1:9)+P(1:9,1:9)')/2;
    
    % Check if a zero velocity update should be done. If so, do the
    % following
    if zupt(t)==true
        
        % Calculate the Kalman filter gain
        K=(P(3:6,1:9)')/(P(3:6,3:6)'+R);
        
        % Calculate the prediction error. Since the detector hypothesis
        % is that the platform has zero velocity, the prediction error is
        % equal to zero minus the estimated velocity.
        z=-x_h(3:6,t);
        
        % Estimation of the perturbations in the estimated navigation
        % states
        dx=K*z;
        
        % Correct the navigation state using the estimated perturbations.
        x_h(1:6,t)=x_h(1:6,t)+dx(1:6);
        q_hat(:,t)=exp_q_L(-dx(7:9),q_hat(:,t));
        
        quat=[q_hat(2:4,t); q_hat(1,t)];
        quats(:, t) = quat;
        
        % Update the filter state covariance matrix P.
        P(1:9,1:9)=P(1:9,1:9)-K*P(3:6,1:9);
        
        % Make sure the filter state covariance matrix is symmetric.
        P(1:9,1:9)=(P(1:9,1:9)+P(1:9,1:9)')/2;
        
    end
    
    for i=1:N_b
        if (~isnan(UWB(i,t)))
            
            %Find the matrix H
            diffs=(x_h(1:3,t)-beacons(:,i));
            dists=sqrt(sum(diffs.^2));
            H=(diffs./dists)';
            
            % Calculate the Kalman filter gain
            K=((H*P(1:3,:))')/(H*P(1:3,1:3)*H'+R_UWB(i,i));
            
            % Calculate the prediction error.
            z=UWB(i,t)-dists';
            
            % Estimation of the perturbations in the estimated navigation
            % states
            dx=K*z;
            
            % Correct the navigation state using the estimated perturbations.
            x_h(1:6,t)=x_h(1:6,t)+dx(1:6);
            q_hat(:,t)=exp_q_L(-dx(7:9),q_hat(:,t));
            
            quat=[q_hat(2:4,t); q_hat(1,t)];
            quats(:, t) = quat;
            
            P=P-K*H*P(1:3,:);
            P=1/2*(P+P');
        end
    end
    
    if (~isnan(y_mag(1,t)))
        if mod(t,101)==0
            %KF meas update
            NablaPhi=[eye(3),Nabla_Phi3D(x_h(1:3,t),N_m,xl,xu,yl,yu,zl,zu,Indices)];
            f=NablaPhi*m;
            J=reshape(JacobianPhi3D(x_h(1:3,t),N_m,xl,xu,yl,yu,zl,zu,Indices),9,N_m)*m(4:end);
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
            x_h(1:6,t)=x_h(1:6,t)+eta_t(1:6);
            q_hat(:,t)=exp_q_L(eta_t(7:9),q_hat(:,t));
            m=m+eta_t(10:end);
            
            P=P-K_t*S_t*K_t';
            P=1/2*(P+P');
            
        end
        
    end
    
end

P_m=P(10:end,10:end);

end





