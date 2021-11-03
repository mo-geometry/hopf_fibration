function angles = getAngles2(nC,nF)

tmpAng = pi/2*ones(nF,1);

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