function write_to_image(imagename,MagNorm,Variance)

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

imwrite(ColorPicture,[imagename,'.png'],'Alpha',Variance);

end