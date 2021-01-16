function angles = getAngles5c(nC,nF)

tmpAng = (pi/4)*ones(nF,1);

% theta = pi/6*ones(nC,1);
phi0 = linspace(0,2*pi,nC)';
inc  = 0.0125*pi*ones(nC,1);
%% return angles
angles.theta(:,1,1) = tmpAng(1)*ones(nC,1);
angles.phi(:,1,1) = phi0;
for ii=2:nF

    %% return angles
    angles.theta(:,1,ii) = tmpAng(ii)*ones(nC,1);
    angles.phi(:,1,ii) = mod(angles.phi(:,1,ii-1)-inc,2*pi);

end