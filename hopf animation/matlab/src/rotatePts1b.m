function pts_new = rotatePts1b(pts)


[nC,~,nF] = size(pts);
pts_new = zeros(size(pts));


rotAngle = linspace(pi/4,pi/2,nF);
for ii=1:nF
    for jj=1:nC
        rotor = rodrigues2([0,1,0],rotAngle(ii));
        pts_new(jj,:,ii) = pts(jj,:,ii)*rotor';
    end    
end

