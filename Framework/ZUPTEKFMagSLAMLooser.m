function [x_h, q_hat, m, P_m]=ZUPTEKFMagSLAMLooser(u,zupt,sigma_y,Lambda,Indices,sigma_lin,N_m,xl,xu,yl,yu,zl,zu)

N=length(u);
alpha=100;
R_acc =alpha*(0.009).^2*eye(3)';
R_gyr =alpha*(0.5*pi/180).^2*eye(3)'; % [rad/s]
Q=[R_acc, zeros(3); zeros(3), R_gyr];
P=zeros(N_m+12,N_m+12);
P(1:6,1:6)=alpha*(1e-10)*eye(6);
P(7:9,7:9)=alpha*(0.1*pi/180).^2*eye(3);
R=alpha*(0.01).^2*eye(4)';
%R(1,1)=100000;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        Run the filter algorithm          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=2:N
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %           Time  Update         %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    u_h=u(:,k);
    
    % Update the navigation equations.
    q_hat(:, k)=exp_q_R(Ts*u_h(4:6),q_hat(:,k-1));
    
    %% Gravity vector
    
    % Convert quaternion to a rotation matrix
    Rb2t=quat2Rot(q_hat(:, k));
    
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
    x_h(1:6,k)=A*x_h(1:6,k-1)+B*acc_t;
    
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
    G=[O O; Ts*quat2Rot(q_hat(:, k)) O; O -Ts*quat2Rot(q_hat(:, k))];
    
    % Update the filter state covariance matrix P.
    P(1:9,1:9)=F*P(1:9,1:9)*F'+G*Q*G';
    
    % Make sure the filter state covariance matrix is symmetric.
    P(1:9,1:9)=(P(1:9,1:9)+P(1:9,1:9)')/2;
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %      Zero-velocity update      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Check if a zero velocity update should be done. If so, do the
    % following
    if zupt(k)==true
        
        % Calculate the Kalman filter gain
        K=(P(3:6,1:9)')/(P(3:6,3:6)'+R);
        
        % Calculate the prediction error. Since the detector hypothesis
        % is that the platform has zero velocity, the prediction error is
        % equal to zero minus the estimated velocity.
        z=-x_h(3:6,k);
        
        % Estimation of the perturbations in the estimated navigation
        % states
        dx=K*z;
        
        % Correct the navigation state using the estimated perturbations.
        x_h(1:6,k)=x_h(1:6,k)+dx(1:6);
        q_hat(:,k)=exp_q_L(-dx(7:9),q_hat(:,k));
        
        quat=[q_hat(2:4,k); q_hat(1,k)];
        quats(:, k) = quat;
        
        % Update the filter state covariance matrix P.
        P(1:9,1:9)=P(1:9,1:9)-K*P(3:6,1:9);
        
        % Make sure the filter state covariance matrix is symmetric.
        P(1:9,1:9)=(P(1:9,1:9)+P(1:9,1:9)')/2;
        
        %Update the cross-covariance terms of the covariance matrix P
        %P(10:end,1:9)=P(10:end,1:9)*F';
        %P(1:9,10:end)=F*P(1:9,10:end);
    
        % Make sure the filter state covariance matrix is symmetric.
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

%% SUBFUNCTIONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  funtion [x quat]=orient(u)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [quat]=init_orient(u)

% Under the assumption that the system is stationary during the first 20
% samples, the initial roll and pitch is calculate from the 20 first
% accelerometer readings.
f_u=mean(u(1,1:20));
f_v=mean(u(2,1:20));
f_w=mean(u(3,1:20));

roll=atan2(-f_v,-f_w);
pitch=atan2(f_u,sqrt(f_v^2+f_w^2));

% Set the attitude vector
attitude=[roll pitch 0]';

% Calculate quaternion corresponing to the initial attitude
Rb2t=Rt2b(attitude)';
quat=dcm2q(Rb2t);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  function q=dcm2q(R)
%
%>
%> @brief Function that converts a directional cosine matrix (rotation
%> matrix) in to a quaternion vector.
%>
%> @param[out]    q      Quaternion vector.
%> @param[in]   R      Rotation matrix.
%>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q=dcm2q(R)

T = 1 + R(1,1) + R(2,2) + R(3,3);

if T > 10^-8
    
    S = 0.5 / sqrt(T);
    qw = 0.25 / S;
    qx = ( R(3,2) - R(2,3) ) * S;
    qy = ( R(1,3) - R(3,1) ) * S;
    qz = ( R(2,1) - R(1,2) ) * S;
    
else
    
    if (R(1,1) > R(2,2)) && (R(1,1) > R(3,3))
        
        S = sqrt( 1 + R(1,1) - R(2,2) - R(3,3)) * 2; % S=4*qx
        qw = (R(3,2) - R(2,3)) / S;
        qx = 0.25 * S;
        qy = (R(1,2) + R(2,1)) / S;
        qz = (R(1,3) + R(3,1)) / S;
        
    elseif (R(2,2) > R(3,3))
        
        S = sqrt( 1 + R(2,2) - R(1,1) - R(3,3) ) * 2; %S=4*qy
        qw = (R(1,3) - R(3,1)) / S;
        qx = (R(1,2) + R(2,1)) / S;
        qy = 0.25 * S;
        qz = (R(2,3) + R(3,2)) / S;
        
    else
        
        S = sqrt( 1 + R(3,3) - R(1,1) - R(2,2) ) * 2; % S=4*qz
        qw = (R(2,1) - R(1,2)) / S;
        qx = (R(1,3) + R(3,1)) / S;
        qy = (R(2,3) + R(3,2)) / S;
        qz = 0.25 * S;
        
    end
    
end

%Store in vector
q = [qw qx qy qz]';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  function R=Rt2b(ang)
%
%>
%> @brief Function that calculates the rotation matrix for rotating a
%> vector from coordinate frame t to the coordinate frame b, given a
%> vector of Euler angles.
%>
%> @param[out]  R      Rotation matrix.
%> @param[in]   ang    Euler angles [roll,pitch,heading]
%>
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R=Rt2b(ang)


cr=cos(ang(1));
sr=sin(ang(1));

cp=cos(ang(2));
sp=sin(ang(2));

cy=cos(ang(3));
sy=sin(ang(3));

R=[cy*cp sy*cp -sp;
    -sy*cr+cy*sp*sr cy*cr+sy*sp*sr cp*sr;
    sy*sr+cy*sp*cr -cy*sr+sy*sp*cr cp*cr];




end

