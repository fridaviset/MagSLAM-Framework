%Loads the algorithm settings and the IMU data
close all; clear;

load('2022_01_20_1.mat');

seed=42;
rng(seed);

%Repeat the measurements to simulate repeated laps
u_full=[u_full, u_full, u_full, u_full];
p_opti_full=[p_opti_full, p_opti_full, p_opti_full, p_opti_full];

%Run the zero-velocity detector
g=9.8173;
sigma_a=(1)^2; 
sigma_g=(0.048*pi/180)^2;     
W=5;
gamma=4e5; 
[zupt_temp]=ZUPT(u_full,g,sigma_a,sigma_g,W,gamma);

%Run the Kalman filter
[x_h_ZUPT, quats]=ZUPTEKF(u_full,zupt_temp);

%Perform alignment of gt based on ZUPT-EKF estimate
A=zeros(4);
N=size(x_h_ZUPT,2);

%Assume the acceleration vector is only caused by gravity
for t=1:N
    p_h=[0; x_h_ZUPT(1:3,t)];
    p_o=[0; p_opti_full(:,t)];
    A=A+q_L(p_h)*q_R(p_o);
end
[V,~] = eig(A);
quat=V(:,1);
R=quat2Rot(quat);
p_opti_c=p_opti_full;

for t=1:N
   p_opti_c(1:3,t)=R*p_opti_full(1:3,t);
end

beacons=[1.5 1 0]';
N_b=size(beacons,2);
radius=1.5;
[UWB, R_UWB]=generate_UWB_meas(p_opti_c,beacons,radius);
[x_h_ZUPT_UWB, q_hat_UWB]=ZUPTEKFUWB(u_full,zupt_temp,UWB,R_UWB,beacons);

margin=2;

%Find the domain
%Estimate the dimensions of the problem based on the first round
traj=x_h_ZUPT(1:3,:);
xl=min(traj(1,:))-margin;
xu=max(traj(1,:))+margin;
yl=min(traj(2,:))-margin;
yu=max(traj(2,:))+margin;
zl=min(traj(3,:))-margin;
zu=max(traj(3,:))+margin;

%Find the magnetic field measurements
y_m=u_full(7:9,:);
y_m_norm=sum(y_m.^2);

%Number of basis functions used in Reduced-Rank approximation
N_m=850;

%Magnetic field params

%New params
sigma_SE=0.84;
l_SE=0.81;
sigma_lin=0.72;
sigma_y=0.1;

%Calculate Lambda and the order of indices used in the
%analytic basis functions of the Reduced-Rank Approximation
[Indices, Lambda]=Lambda3D(N_m,xl,xu,yl,yu,zl,zu,sigma_SE,l_SE);

y_mag=u_full(8:9,:);
q_0=quats(:,1);
v_0=zeros(3,1);
p_0=zeros(3,1);
Ts=1./200;
N=size(u_full,2);

[x_h, q_hat, m, P_m]=ZUPTEKFMagSLAM(u_full,zupt_temp,sigma_y,Lambda,Indices,sigma_lin,N_m,xl,xu,yl,yu,zl,zu);

[x_h_all, q_hat, m_all, P_m_all]=ZUPTEKFMagSLAMUWB(u_full,zupt_temp,sigma_y,Lambda,Indices,sigma_lin,N_m,xl,xu,yl,yu,zl,zu,UWB,R_UWB,beacons);

%Plot results

res=0.05; z=0;
[X,Y,Z,pointsVec,PhiVec,NablaPhi3D]=prepare_magnetic_field_plots(xl,xu,yl,yu,zl,zu,N_m,Indices,res,z);
fontsize=13;
f=figure; clf;
linewidth=0.6;
plot3(x_h_ZUPT(1,:),x_h_ZUPT(2,:),x_h_ZUPT(3,:)+10,'b','LineWidth',linewidth);
hold on;
plot3(p_opti_c(1,:),p_opti_c(2,:),p_opti_c(3,:)+10,'r','LineWidth',linewidth)
plot3(x_h_all(1,:),x_h_all(2,:),x_h_all(3,:)+10,'k','LineWidth',linewidth);
view(2);
for i=1:N_b
d = radius*2;
px = beacons(1,i)-radius;
py = beacons(2,i)-radius;
h = rectangle('Position',[px py d d],'Curvature',[1,1],'linewidth',2,'EdgeColor','r','Linestyle','--');

%Smaller circle for center of beacon
scatter3(beacons(1,i),beacons(2,i),beacons(3,i)+10,'r','Linewidth',2);

end
plot_projection_norm_magnetic_field(m_all,P_m_all,NablaPhi3D,X,Y,Z,fontsize);
caxis([0.6 1.5]);
grid off;
axis equal;
xlim([-1.5 6]);
ylim([-2 7.5]);

%Plot pure magnetic field map and save
fontsize=13;
[MagNorm, Variance]=plot_projection_norm_magnetic_field(m_all,P_m_all,NablaPhi3D,X,Y,Z,fontsize);
caxis([0.6 1.5]);
view(2);

A=parula;
%Normalise Y data
MagNorm=MagNorm-min(min(MagNorm));
MagNorm=MagNorm./max(max(MagNorm));
indices=ceil((MagNorm).*254+1);

%Normalise Opacity
Variance=Variance-min(min(Variance));
Variance=Variance./max(max(Variance));

%Quick calculation
ColorPicture=reshape(A(indices,:),size(MagNorm,1),size(MagNorm,2),3);

%Make a squeezed image
MagNormSqueezed=zeros(size(MagNorm).*10);
VarianceSqueezed=zeros(size(MagNorm).*10);

%Flip the map
MagNorm=flip(MagNorm,1);
Variance=flip(Variance,1);

%Iterate over all the new coordinates
imagewidth=size(MagNormSqueezed);
offset=[0; 0];
domain_size=size(MagNorm);

for i=1:imagewidth(1)
    for j=1:imagewidth(2)
        pos=[i; (j-imagewidth(2)./2)./((i)./imagewidth(1))+imagewidth(2)./2];
        rel_pos=pos-offset;
        unrounded=domain_size'./imagewidth'.*rel_pos;
        shift=round(unrounded);
        dimension_1_ok=(shift(1)<=size(MagNorm,1)) && (shift(1)>=1);
        dimension_2_ok=(shift(2)<=size(MagNorm,2)) && (shift(2)>=1);
        if dimension_1_ok && dimension_2_ok
            MagNormSqueezed(i,j)=MagNorm(shift(1),shift(2));
            VarianceSqueezed(i,j)=Variance(shift(1),shift(2));
        end
    end
end

imagename='MagneticFieldMap';
%Normalise Y data
indices_squeezed=ceil((MagNormSqueezed).*254+1);
ColorPictureSqueezed=reshape(A(indices_squeezed,:),size(MagNormSqueezed,1),size(MagNormSqueezed,2),3);
imwrite(ColorPictureSqueezed,[imagename,'.png'],'Alpha',VarianceSqueezed);

error_SLAM=Check_error_for_aligned_GT(x_h,p_opti_c);
error_ZUPT=Check_error_for_aligned_GT(x_h_ZUPT,p_opti_c);
error_UWB=Check_error_for_aligned_GT(x_h_ZUPT_UWB,p_opti_c);
error_ALL=Check_error_for_aligned_GT(x_h_all,p_opti_c);

f=figure; clf;
color_SLAM=[0, 184, 129]./255;
color_UWB=[255, 61, 90]./255;
linewidth=1.5;
plot(error_SLAM,'Color',color_SLAM,'LineWidth',linewidth);
hold on;
plot(error_ZUPT,'b','LineWidth',linewidth);
plot(error_UWB,'Color',color_UWB,'LineWidth',linewidth);
plot(error_ALL,'k','LineWidth',linewidth);
legend({'SLAM','ZUPT','ZUPT-UWB','SLAM-UWB'},'Interpreter','Latex','Fontsize',fontsize,'Location','northwest');
ylabel('Position RMSE (m)','Interpreter','Latex','Fontsize',fontsize);
xlabel('Timestep','Interpreter','Latex','Fontsize',fontsize);

fileID=fopen('results.txt','w');
fprintf(fileID,['IMU+ZUPT RMSE: ',num2str(sqrt(mean(error_ZUPT.^2))),'\n']);
fprintf(fileID,['IMU+ZUPT+Magnetic field RMSE: ',num2str(sqrt(mean(error_SLAM.^2))),'\n']);
fprintf(fileID,['IMU+ZUPT+UWB RMSE: ',num2str(sqrt(mean(error_UWB.^2))),'\n']);
fprintf(fileID,['IMU+ZUPT+Magnetic field+UWB RMSE: ',num2str(sqrt(mean(error_ALL.^2))),'\n']);
fclose(fileID);
