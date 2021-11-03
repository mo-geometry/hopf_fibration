function pts_new = rotatePts9ac(pts)


[nC,~,nF] = size(pts);
pts_new = zeros(size(pts));


rotAngle = linspace(0,6*pi,nF);
for ii=1:nF
    for jj=1:nC
        rotor = rodrigues2([1,0,0],rotAngle(ii));
        pts_new(jj,:,ii) = pts(jj,:,ii)*rotor';
    end    
end

