function error=Check_error_for_aligned_GT(x_h,p_opti)

%Perform alignment of gt based on ZUPT-EKF estimate
A=zeros(4);
N=size(x_h,2);
for t=1:N
    p_h=[0; x_h(1:3,t)];
    p_o=[0; p_opti(:,t)];
    A=A+q_L(p_h)*q_R(p_o);
end
[V,~] = eig(A);
quat=V(:,1);
R=quat2Rot(quat);
p_opti_c=p_opti;

for t=1:N
   p_opti_c(1:3,t)=R*p_opti(1:3,t);
end

error=sum((p_opti_c-x_h(1:3,:)).^2);

end