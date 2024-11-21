function estP=WienerNIG_L_MCMC(Obs)
% Random effects-Wiener process model 
% X(t)=v*t+w^(1/2)B(t)
% v~N(mu,kappa2), w~IG(eta, zeta^(-1))
% obs={{T1,X1},{T2,X2},...{Tn,Xn}}.
% EM setting
maxEMite=100;
maxVBite=100;
EMtol=1e-6;
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
%% Initialization
mu_0=mean(v_V);
kappa2_0=var(v_V);
eta_0=mean(w_V);
zeta_0=(mean(1./w_V)-1/eta_0);
%% EM iteration
for ite=1:maxEMite
    % E-step using VB approximation
    %{
    if ite==1
        Ev_V0=mu_0*ones(1,N);
        Ev2_V0=(kappa2_0+mu_0^2)*ones(1,N);
    else
        Ev_V0=Ev_V1;
        Ev2_V0=Ev2_V1;
    end
    a_V=1/zeta_0/eta_0^2;
    b_V0=sumdX2overdT_V-2*Ev_V0.*Xend_V+Ev2_V0.*Tend_V+1/zeta_0;
    p_V=-(m_V+1)/2;
    tmpR=besselk(p_V+1,sqrt(a_V.*b_V0),1)./besselk(p_V,sqrt(a_V.*b_V0),1);
    if any(isnan(tmpR)|isinf(tmpR))
        tmpID=find(isnan(tmpR)|isinf(tmpR));
        for kk=1:length(tmpID)
            tmpindID=tmpID(kk);
            tmpR(tmpindID)=Rbesselk(p_V(tmpindID)+1,sqrt(a_V*b_V0(tmpindID)));
    %             tmpR2(tmpindID)=tmpR(tmpindID)*Rbesselk(P_flatten+2,sqrt(A_flatten*B_flatten(tmpindID)));
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
    %             tmpR2(tmpindID)=tmpR(tmpindID)*Rbesselk(P_flatten+2,sqrt(A_flatten*B_flatten(tmpindID)));
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
    %}
    % E-step using MCMC sampling
    Ev_V1_MCMC=zeros(1,N);    
    Ev2_V1_MCMC=zeros(1,N);
    Ew_V1_MCMC=zeros(1,N);
    Ew1_V1_MCMC=zeros(1,N);
    for i=1:N
        logpdf = @(vw) -m_V(i)/2*log(vw(2))-(sumdX2overdT_V(i)-2*vw(1)*Xend_V(i)+vw(1)^2*Tend_V(i))/2/vw(2)...
            -(vw(1)-mu_0).^2/2/kappa2_0-3/2*log(vw(2))-(vw(2)-eta_0)^2/2/zeta_0/eta_0^2/vw(2); % Target distribution
        logproppdf = @(vw1,vw0)  -(vw1(1)-mu_0).^2/2/kappa2_0+(1/eta_0/zeta_0-1)*log(vw1(2))-vw1(2)/(eta_0^2*zeta_0);
%         proprnd = @(vw)  [randn(1)*sqrt(kappa2_0)+mu_0,random('inversegaussian',eta_0,1/zeta_0)];
        proprnd = @(vw0)  [randn(1)*sqrt(kappa2_0)+mu_0,gamrnd(1/eta_0/zeta_0,eta_0^2*zeta_0)];
        [vw_smpl,vm_r] = mhsample([mu_0,eta_0],1000,'logpdf',logpdf,'proprnd',proprnd,'logproppdf',logproppdf,'burnin',1000,'thin',10);
%         [vw_smpl,vm_r] = mhsample([mu_0,eta_0],1000,'logpdf',logpdf,'proprnd',proprnd,'symmetric',1,'burnin',10,'thin',10);
%         disp(vm_r);
        Ev_V1_MCMC(i)=mean(vw_smpl(:,1));
        Ev2_V1_MCMC(i)=mean(vw_smpl(:,1).^2);
        Ew_V1_MCMC(i)=mean(vw_smpl(:,2));
        Ew1_V1_MCMC(i)=mean(1./vw_smpl(:,2));
    end
    % M-step
    mu_1=mean(Ev_V1_MCMC);
    kappa2_1=mean(Ev2_V1_MCMC-mu_1^2);
    eta_1=mean(Ew_V1_MCMC);
    zeta_1=mean(Ew1_V1_MCMC-1/eta_1);
    % check convergence
    EMdif=abs(mu_1-mu_0)+abs(kappa2_1-kappa2_0)+abs(eta_1-eta_0)+abs(zeta_1-zeta_0);
    if EMdif<EMtol
        break;
    else
        mu_0=mu_1;
        kappa2_0=kappa2_1;
        eta_0=eta_1;
        zeta_0=zeta_1;
    end
end
%
estP=struct('mu',mu_1,'kappa2',kappa2_1,'eta',eta_1,'zeta',zeta_1,'EMite',ite,'EMdif',EMdif);
loglike=WienerNIG_lnL(Obs,estP);
estP.loglike=loglike;
