clear; close all;
MagNorm=-peaks(200);
Variance=linspace(0,8,size(MagNorm,1))'*ones(1,size(MagNorm,2));


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
imagename='testImage';
imwrite(ColorPicture,[imagename,'.png'],'Alpha',Variance);

%Try to make a squeezed image
MagNormSqueezed=zeros(size(MagNorm));
VarianceSqueezed=zeros(size(MagNorm));

%Iterate over all the old coordinates
imagewidth=size(MagNorm);
offset=[0; 0];
domain_size=size(MagNorm);

for i=1:imagewidth(1)
    for j=1:imagewidth(2)
        pos=[i; (j-100)./(i./200)+100];
        rel_pos=pos-offset;
        unrounded=imagewidth'./domain_size'.*rel_pos;
        shift=floor(unrounded);
        dimension_1_ok=(shift(1)<=size(MagNorm,1)) && (shift(1)>=1);
        dimension_2_ok=(shift(2)<=size(MagNorm,2)) && (shift(2)>=1);
        if dimension_1_ok && dimension_2_ok
            MagNormSqueezed(i,j)=MagNorm(shift(1),shift(2));
            VarianceSqueezed(i,j)=Variance(shift(1),shift(2));
        end
    end
end

imagename='testImage2';
%Normalise Y data
indices_squeezed=ceil((MagNormSqueezed).*254+1);
ColorPictureSqueezed=reshape(A(indices_squeezed,:),size(MagNorm,1),size(MagNorm,2),3);
imwrite(ColorPictureSqueezed,[imagename,'.png'],'Alpha',VarianceSqueezed);

%surf(X,0.5*(yl+yu)+(Y-0.5*(yl+yu)).*0.1.*(X-xl),Z,Norm,'AlphaData',1./Var,'FaceAlpha','flat','EdgeColor','none');