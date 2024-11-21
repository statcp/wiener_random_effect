function sam=WienerNIG_L_sam(para)
% WienerNIG_L_sam is used to generate samples from the Wiener process model
% X(t)=v*t+w^(1/2)*B(t), where v~N(mu,kappa2), w~IG(eta,zeta^(-1))
mu0=para.mu;
kappa2_0=para.kappa2;
kappa0=sqrt(kappa2_0);
eta0=para.eta;
zeta0=para.zeta;
T=para.T; %this is the observation scheme of the original data
          %it is a cell vector of size n*1, where n is the number of
          % observations
n=length(T);
sam=cell(n,1);
v_V=randn(n,1)*kappa0+mu0;
w_V=random('inversegaussian',eta0,1/zeta0,n,1);
for i=1:n
    dT=diff([0;T{i}]);
    dX=v_V(i)*dT+sqrt(w_V(i)*dT).*randn(size(dT));
    X=cumsum(dX);
    sam{i}={T{i},X};
end