function [zupt]=ZUPT(u,g,sigma_a,sigma_g,W,gamma)

% Allocate memmory
zupt=zeros(1,length(u));

N=length(u);
T=zeros(1,N-W+1);

for k=1:N-W+1   
    ya_m=mean(u(1:3,k:k+W-1),2);   
    for l=k:k+W-1
        tmp=u(1:3,l)-g*ya_m/norm(ya_m);
        T(k)=T(k)+u(4:6,l)'*u(4:6,l)/sigma_g+tmp'*tmp/sigma_a;    
    end    
end

T=T./W;

% Check if the test statistics T are below the detector threshold. If so, 
% chose the hypothesis that the system has zero velocity 
for k=1:length(T)
    if T(k)<gamma
       zupt(k:k+W-1)=ones(1,W); 
    end    
end

end



