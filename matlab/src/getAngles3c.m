function angles = getAngles3c(nC,nF)


tmpAng = pi/2*ones(nF,1);

phi0 = linspace(0,pi,nC)';
inc  = 0.025*pi*ones(nC,1);
inc0  = 2*pi*ones(nC,1);

%% return angles
angles.theta(:,1,1) = tmpAng(1)*ones(nC,1);
angles.phi(:,1,1) = phi0;

lim1 = floor(nF/8);
for ii=2:lim1

    %% return angles
    angles.theta(:,1,ii) = linspace(pi/2-pi/2*(ii)/lim1,pi/2+pi/2*(ii)/lim1,nC)';
    angles.phi(:,1,ii) = phi0;

end

for ii=lim1+1:4*lim1

    %% return angles
    angles.theta(:,1,ii) = linspace(0,pi,nC)';
    angles.phi(:,1,ii) =  mod(linspace(0,pi+3*pi*(ii-lim1)/(3*lim1),nC)'+0.5*(ii-lim1)/(3*lim1)*inc0,2*pi);
%
end

for ii=4*lim1+1:5*lim1

    %% return angles
    angles.theta(:,1,ii) = linspace(0,pi,nC)';
    angles.phi(:,1,ii) =  mod(angles.phi(:,1,ii-1)+2*inc,2*pi);

end

for ii=5*lim1+1:nF

    %% return angles
    angles.theta(:,1,ii) = linspace(pi/2-pi/2*(nF-ii)/(3*lim1),pi/2+pi/2*(nF-ii)/(3*lim1),nC)';
    angles.phi(:,1,ii) =  mod(angles.phi(:,1,ii-1)+2*inc,2*pi);

end




