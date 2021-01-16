function pts_new = rotatePts7(pts)


[nC,~,nF] = size(pts);
pts_new = zeros(size(pts));



seq1 = floor(nF/8);
rotAngle1 = linspace(0,pi/2,seq1);

for ii=1:seq1
    for jj=1:nC
        rotor = rodrigues2([0,0,1],rotAngle1(ii));
        pts_new(jj,:,ii) = pts(jj,:,ii)*rotor';
    end    
end
rotor2 = rodrigues2([0,0,1],pi/2);
rotAngle2 = linspace(0,2*pi,4*seq1);
for ii=seq1+1:5*seq1
    for jj=1:nC
        rotor = rodrigues2([1,0,0],rotAngle2(ii-seq1));
        pts_new(jj,:,ii) = pts(jj,:,ii)*rotor2'*rotor';
    end    
end
rotor3 = rodrigues2([1,0,0],2*pi);
rotAngle3 = linspace(pi/2,2*pi,3*seq1);
for ii=5*seq1+1:nF
    for jj=1:nC
        rotor = rodrigues2([0,0,1],rotAngle3(ii-5*seq1));
        pts_new(jj,:,ii) = pts(jj,:,ii)*rotor3'*rotor';
    end    
end