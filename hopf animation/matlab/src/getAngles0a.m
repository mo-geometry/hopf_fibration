function angles = getAngles0a(nF)

phi0=linspace(0,2*pi,nF)';
phi1=linspace(0,4*pi,nF)';
phi2=linspace(0,8*pi,nF)';
theta0=linspace(0,2*pi,4*nF)';
angles.theta(1,1,:)=[pi/2*ones(nF,1);pi/2*sin(theta0)+pi/2;pi/2*ones(nF,1)];
angles.phi(1,1,:)=[phi0;phi1;phi2;phi2;phi1;phi0];