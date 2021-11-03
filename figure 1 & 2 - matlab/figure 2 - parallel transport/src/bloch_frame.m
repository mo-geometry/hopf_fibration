function [fig, eV, eB, TanVec, Theta, Phi] = bloch_frame(UofT, tvec, theta0, phi0)
theta0 = theta0 * pi / 180;
phi0 = phi0 * pi / 180;
%% Extracting the quaternion
Dt = tvec(2)-tvec(1); % time step
a=UofT(:,1);b=UofT(:,2);c=UofT(:,3);d=UofT(:,4); % quaternion elements

%% Derivative of the quaternion
da = gradient(a,Dt);
db = gradient(b,Dt);
dc = gradient(c,Dt);
dd = gradient(d,Dt);

d2a = gradient(da,Dt);
d2b = gradient(db,Dt);
d2c = gradient(dc,Dt);
d2d = gradient(dd,Dt);


%% Evolving the bloch vector rN
r0 = [sin(theta0)*cos(phi0);sin(theta0)*sin(phi0);cos(theta0)];
fig.blochVector = zeros(length(tvec),3);

for kk=1:length(tvec)    
Rotor=[ a(kk)^2+b(kk)^2-c(kk)^2-d(kk)^2,    2*(b(kk)*c(kk)-a(kk)*d(kk)),        2*(b(kk)*d(kk)+a(kk)*c(kk));
        2*(b(kk)*c(kk)+a(kk)*d(kk)),        a(kk)^2-b(kk)^2+c(kk)^2-d(kk)^2,    2*(c(kk)*d(kk)-a(kk)*b(kk));
        2*(b(kk)*d(kk)-a(kk)*c(kk)),        2*(c(kk)*d(kk)+a(kk)*b(kk)),        a(kk)^2-b(kk)^2-c(kk)^2+d(kk)^2];
fig.blochVector(kk,:)=Rotor*r0;
end

%% Hamiltonian Elements
Hi = 2*(a.*db - da.*b - dc.*d + c.*dd);
Hj = 2*(a.*dc - da.*c + db.*d - b.*dd);
Hk = 2*(a.*dd - da.*d - db.*c + b.*dc);
H = [Hi, Hj, Hk];

Hi_dt = 2*(a.*d2b - d2a.*b - d2c.*d + c.*d2d);
Hj_dt = 2*(a.*d2c - d2a.*c + d2b.*d - b.*d2d);
Hk_dt = 2*(a.*d2d - d2a.*d - d2b.*c + b.*d2c);

%% The bloch vector
ri = fig.blochVector(:,1);
rj = fig.blochVector(:,2);
rk = fig.blochVector(:,3);
R = [ri, rj, rk];

ri_dt = gradient(ri,Dt);
rj_dt = gradient(rj,Dt);
rk_dt = gradient(rk,Dt);

%% Phase
Ddyn = Hi.*ri+Hj.*rj+Hk.*rk;

H_dot_H =  Hi.*Hi + Hj.*Hj + Hk.*Hk;
dR_dot_dR = ri_dt .* ri_dt + rj_dt .* rj_dt + rk_dt .* rk_dt;
dH_dot_dR = Hi_dt .* ri_dt + Hj_dt .* rj_dt + Hk_dt .* rk_dt;
denom = (H_dot_H - Ddyn.^2);
Dgeo = (dH_dot_dR - dR_dot_dR .* Ddyn) ./ denom;

%% Bloch angles
dTheta = (Hj .* ri - Hi .* rj) ./ sqrt(ri.^2 + rj.^2);
dPhi = Hk - rk .* (Hi .* ri + Hj .* rj) ./ (ri.^2 + rj.^2);

Theta = cumsum(dTheta)*Dt - dTheta(1)*Dt + theta0;
Phi = cumsum(dPhi)*Dt - dPhi(1)*Dt + phi0;


%% Initialisation
fig.geometricPhase = zeros(size(Dgeo));
fig.dynamicPhase = zeros(size(Ddyn)); 

%% Integrating the angles
fig.geometricPhase=cumsum(Dgeo)*Dt;
fig.geometricPhase=fig.geometricPhase-fig.geometricPhase(1);

fig.dynamicPhase=cumsum(Ddyn)*Dt;
fig.dynamicPhase=fig.dynamicPhase-fig.dynamicPhase(1);

%% moving frame
eV = cross(H, R) ./ sqrt(denom);
eB = (H - R .* Ddyn) ./ sqrt(denom);

%% Numerical integration of the tangent vector
TanVec=zeros(2,length(tvec));
Rotor=zeros(2,2,length(tvec));
TanVec(:,1)=[1/sqrt(2);1/sqrt(2)];

for kk=1:length(tvec)
    Rotor(1:2,1:2,kk)=[cos(fig.geometricPhase(kk)),sin(fig.geometricPhase(kk));-sin(fig.geometricPhase(kk)),cos(fig.geometricPhase(kk))];
end

for kk=1:length(tvec)-1
%     TanVec(:,kk+1)=expm(DSop(:,:,kk)*Dt)*TanVec(:,kk);
    TanVec(:,kk+1)=Rotor(1:2,1:2,kk)'*TanVec(:,1);
end

TanVec=TanVec';