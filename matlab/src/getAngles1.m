function angles = getAngles1(nC,nF)

tmp = linspace(0,2*pi,nF);
tmpAng = pi*(0.5+0.5*(5/6)*sin(tmp));

phi0 = linspace(0,2*pi,nC)';
inc  = 0.025*pi*ones(nC,1);
%% return angles
angles.theta(:,1,1) = tmpAng(1)*ones(nC,1);
angles.phi(:,1,1) = phi0;
for ii=2:nF

    %% return angles
    angles.theta(:,1,ii) = tmpAng(ii)*ones(nC,1);
    angles.phi(:,1,ii) = mod(angles.phi(:,1,ii-1)+inc,2*pi);

end