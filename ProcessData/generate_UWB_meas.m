function [UWB, R_UWB]=generate_UWB_meas(p_opti_full,beacons,radius)
UWB=NaN(size(beacons,2),size(p_opti_full,2));
N_b=size(beacons,2);
N=size(p_opti_full,2);
R_UWB=0.1*eye(N_b);

for t=1:20:40000
    for i=1:size(beacons,2)
        dist=norm(p_opti_full(:,t)-beacons(:,i));
        if dist<radius     
            UWB(i,t)=dist;
        end
    end
    UWB(:,t)=UWB(:,t)+mvnrnd(zeros(N_b,1),R_UWB)';
end


end