function loglike=WienerNIG_lnL(Obs,estP)
% Normal and IG random effects
%obs={{T1,X1},{T2,X2},...{Tn,Xn}}.
maxVBite=100;
VBto1=1e-6;
%
n=length(Obs);
n_V=zeros(1,n); %number of observations in each group
for i=1:n
    [~,n_V(i)]=size(Obs{i}{2});
end
N=sum(n_V);
m_V=zeros(1,N); % number of observations for each unit
Tend_V=zeros(1,N); % the end of observation for each unit
Xend_V=zeros(1,N); % the degration level at the end of observation for each unit
sumdX2overdT_V=zeros(1,N); % the sum of delta(X) for each unit
sumlndT_V=zeros(1,N); % the sum of log(delta(T)) for each unit
v_V=zeros(1,N);
w_V=zeros(1,N);
for i=1:n
    tmpT=Obs{i}{1};
    tmpX=Obs{i}{2};
    if size(tmpT,2)~=1
        errormessage('The input format is incorrect!');
        return;
    end
    range_i=(sum(n_V(1:i-1))+1):sum(n_V(1:i));% the indices of the units in the i-th group
    m_V(range_i)=length(tmpT);
    Tend_V(range_i)=tmpT(end,:);
    Xend_V(range_i)=tmpX(end,:);
    tmpdT=tmpT-[0;tmpT(1:end-1)];
    tmpdX=tmpX-[zeros(1,n_V(i));tmpX(1:end-1,:)];
    sumdX2overdT_V(range_i)=sum((tmpdX.^2)./(tmpdT*ones(1,n_V(i))),1);
    sumlndT_V(range_i)=sum(log(tmpdT));
    v_V(range_i)=Xend_V(range_i)./Tend_V(range_i);
    w_V(range_i)=var(tmpdX-tmpdT*v_V(range_i),0,1);
end
%% Parameter assign
mu_0=estP.mu;
kappa2_0=estP.kappa2;
eta_0=estP.eta;
zeta_0=estP.zeta;
%% one VB approximation
Ev_V0=mu_0*ones(1,N);
Ev2_V0=kappa2_0*ones(1,N)+Ev_V0.^2;
a_V=1/zeta_0/eta_0^2;
b_V0=sumdX2overdT_V-2*Ev_V0.*Xend_V+Ev2_V0.*Tend_V+1/zeta_0;
p_V=-(m_V+1)/2;
tmpR=besselk(p_V+1,sqrt(a_V.*b_V0),1)./besselk(p_V,sqrt(a_V.*b_V0),1);
if any(isnan(tmpR)|isinf(tmpR))
    tmpID=find(isnan(tmpR)|isinf(tmpR));
    for kk=1:length(tmpID)
        tmpindID=tmpID(kk);
        tmpR(tmpindID)=Rbesselk(p_V(tmpindID)+1,sqrt(a_V*b_V0(tmpindID)));
    end
end
%Ew_V0=sqrt(b_V0./a_V).*tmpR;
Ew1_V0=sqrt(a_V./b_V0).*tmpR-2*p_V./b_V0;
for vbite=1:maxVBite
    Ev_V1=(Xend_V.*Ew1_V0+mu_0/kappa2_0)./(Tend_V.*Ew1_V0+1/kappa2_0);
    Ev2_V1=1./(Tend_V.*Ew1_V0+1/kappa2_0)+Ev_V1.^2;
    b_V1=sumdX2overdT_V-2*Ev_V1.*Xend_V+Ev2_V1.*Tend_V+1/zeta_0;
    tmpR=besselk(p_V+1,sqrt(a_V.*b_V1),1)./besselk(p_V,sqrt(a_V.*b_V1),1);
    if any(isnan(tmpR)|isinf(tmpR))
        tmpID=find(isnan(tmpR)|isinf(tmpR));
        for kk=1:length(tmpID)
            tmpindID=tmpID(kk);
            tmpR(tmpindID)=Rbesselk(p_V(tmpindID)+1,sqrt(a_V*b_V1(tmpindID)));
        end
    end
    %         Ew_V1=sqrt(b_V1./a_V).*tmpR;
    Ew1_V1=sqrt(a_V./b_V1).*tmpR-2*p_V./b_V1;
    if sum(abs(Ew1_V1-Ew1_V0))...
            +sum(abs(Ev_V1-Ev_V0))+sum(abs(Ev2_V1-Ev2_V0))<VBto1
        break;
    else
        Ev_V0=Ev_V1;
        Ev2_V0=Ev2_V1;
        %             Ew_V0=Ew_V1;
        Ew1_V0=Ew1_V1;
    end
end
Ew_V1=sqrt(b_V1./a_V).*tmpR;
S2v_V1=Ev2_V1-Ev_V1.^2;
%
loglike=-m_V/2*log(2*pi)-1/2*sumlndT_V-1/2*Ew1_V1.*(sumdX2overdT_V-2*Ev_V1.*Xend_V+Ev2_V1.*Tend_V)...
    -1/2*log(kappa2_0)-1/2/kappa2_0*(Ev2_V1-2*mu_0*Ev_V1+mu_0^2)...
    -1/2*log(zeta_0)-1/2*log(2*pi)-1/zeta_0/2/eta_0^2*(Ew_V1-2*eta_0+eta_0^2*Ew1_V1)...
    +1/2*(log(S2v_V1)+1)...
    -p_V/2.*log(a_V./b_V1)+log(2)+lnbesselk(p_V,sqrt(a_V.*b_V1))+(a_V.*Ew_V1+b_V1.*Ew1_V1)/2;
loglike=sum(loglike);
    function y=lnbesselk(p,x)
        y=log(besselk(p,x,1))-x;
    end
end

