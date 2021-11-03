function Img_pad = padding(Img,KK)

hw=size(Img);
a=size(KK);
a=floor(a(1)/2);

Img_pad = zeros(hw(1)+2*a,hw(2)+2*a);
Img_pad(1+hw(1):end-hw(1),1+hw(2):end-hw(2)) = Img;