function Rt=WienerNIG_Rt(para,t,Df,funcLam)
% the reliability function at t for degradation model
% X(t)=vLam(t)+omega^(1/2)B(Lam(t)) with respect to the threshold D
%
if nargin<4
    tmpLam=@(t) t;
else
    tmpLam=@(t) funcLam(t,para.alpha);
end
% the pdf of Lambda(T) is f(x;v,w)=D/sqrt(2*pi*w*x^3)*exp(-(v*x-D)^2/2/w/x)
% after integrating w out,
% f(x;v)=D*sqrt(zeta)/(2*pi*sqrt(x^3))*exp(zeta/eta)*2*Kp(sqrt(ab))*(b/a)^(p/2)
% where a=zeta/eta, b=(v*x-D)^2/x+zeta
% p=-1
a=para.zeta/para.eta^2;
b=@(x,v) (v.*x-Df).^2./x+para.zeta; 
tmpft=@(x,v) Df*sqrt(para.zeta)./(pi*sqrt(x.^3)).*exp(para.zeta/para.eta).*besselk(-1,sqrt(a*b(x,v))).*(a./b(x,v)).^(1/2);
t_v=t(:);
Rt=zeros(size(t_v));
for i=1:length(t_v)
    tmpF=@(x,v) tmpft(x,v).*pdf('normal',v,para.mu,para.kappa2.^0.5);
    Rt(i)=integral2(tmpF,0,tmpLam(t_v(i)),para.mu-para.kappa2^0.5*6,para.mu+para.kappa2^0.5*6);
end
Rt=1-reshape(Rt,size(t));