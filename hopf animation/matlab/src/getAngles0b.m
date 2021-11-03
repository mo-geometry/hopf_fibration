function angles = getAngles0b(nC,nF)

tmpAng = pi/2*ones(nF,1);
%% return angles
angles.theta(:,1,1) = tmpAng(1)*ones(nC,1);
angles.phi(:,1,1) = zeros(nC,1);
for ii=2:nF
    %% return angles
    angles.theta(:,1,ii) = tmpAng(ii)*ones(nC,1);
    angles.phi(:,1,ii) = linspace(0,2*pi*ii/nF,nC)';
end