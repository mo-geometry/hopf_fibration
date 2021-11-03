function angles = getAngles3b(nC,nF)


tmpAng = pi/2*ones(nF,1);

phi0 = linspace(0,pi,nC)';
%% return angles
angles.theta(:,1,1) = tmpAng(1)*ones(nC,1);
angles.phi(:,1,1) = phi0;
for ii=2:nF

    %% return angles
    angles.theta(:,1,ii) = tmpAng(ii)*ones(nC,1);
    angles.phi(:,1,ii) = phi0;

end

