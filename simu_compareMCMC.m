%test WienerNIG
mu_0=2.6;
kappa_0=1.0;
eta_0=5.4;
zeta_0=0.8;
kappa2_0=kappa_0^2;
para_0=struct('mu',mu_0,'kappa2',kappa2_0,'eta',eta_0,'zeta',zeta_0);
p_cover=0.9;
%
Nsim=1e3;
%
n=[25,50,100];
m=[10,20,40];
%
Simu_mu=zeros(Nsim,9);
Simu_kappa2=zeros(Nsim,9);
Simu_eta=zeros(Nsim,9);
Simu_zeta=zeros(Nsim,9);
Simu_EMite=zeros(Nsim,9);
Simu_dt=zeros(Nsim,9);
%
Simu_mu_MCMC=zeros(Nsim,9);
Simu_kappa2_MCMC=zeros(Nsim,9);
Simu_eta_MCMC=zeros(Nsim,9);
Simu_zeta_MCMC=zeros(Nsim,9);
Simu_EMite_MCMC=zeros(Nsim,9);
Simu_dt_MCMC=zeros(Nsim,9);
%
count_nm=0;
for n1=1:length(n)
    for m1=1:length(m)
        count_nm=count_nm+1;
        fprintf([num2str(count_nm),repmat('.',1,Nsim-1)]);
        fprintf('\n\r');
        m_V=m(m1)*ones(n(n1),1);
        T=(1:m(m1))';
        T=T.^0.8;
        T_cell=cell(n(n1),1);
        for ii=1:n(n1)
            T_cell{ii}=T;
        end
        para_0.T=T_cell;
        parfor rep=1:Nsim
            % generate samples
            obs_test=WienerNIG_L_sam(para_0);
            % point estimates
            s=tic;
            estP_test=WienerNIG_L(obs_test);
            Simu_dt(rep,count_nm)=toc(s);
            Simu_mu(rep,count_nm)=estP_test.mu;
            Simu_kappa2(rep,count_nm)=estP_test.kappa2;
            Simu_eta(rep,count_nm)=estP_test.eta;
            Simu_zeta(rep,count_nm)=estP_test.zeta;
            Simu_EMite(rep,count_nm)=estP_test.EMite;
            % MCMC
            s=tic;
            estP_MCMC=WienerNIG_L_MCMC(obs_test);
            Simu_dt_MCMC(rep,count_nm)=toc(s);
            Simu_mu_MCMC(rep,count_nm)=estP_MCMC.mu;
            Simu_kappa2_MCMC(rep,count_nm)=estP_MCMC.kappa2;
            Simu_eta_MCMC(rep,count_nm)=estP_MCMC.eta;
            Simu_zeta_MCMC(rep,count_nm)=estP_MCMC.zeta;
            Simu_EMite_MCMC(rep,count_nm)=estP_MCMC.EMite;
            fprintf('\b|\n');
        end
        fprintf('\r');
    end
end
Simu_kappa=Simu_kappa2.^0.5;
bias_MSE=[mean(Simu_mu,1)-mu_0;mean((Simu_mu-mu_0).^2,1);...
    mean(Simu_kappa)-kappa_0;mean((Simu_kappa-kappa_0).^2);...
    mean(Simu_eta)-eta_0;mean((Simu_eta-eta_0).^2);...
    mean(Simu_zeta)-zeta_0;mean((Simu_zeta-zeta_0).^2)]';
Simu_kappa_MCMC=Simu_kappa2_MCMC.^0.5;
bias_MSE_MCMC=[mean(Simu_mu_MCMC,1)-mu_0;mean((Simu_mu_MCMC-mu_0).^2,1);...
    mean(Simu_kappa_MCMC)-kappa_0;mean((Simu_kappa_MCMC-kappa_0).^2);...
    mean(Simu_eta_MCMC)-eta_0;mean((Simu_eta_MCMC-eta_0).^2);...
    mean(Simu_zeta_MCMC)-zeta_0;mean((Simu_zeta_MCMC-zeta_0).^2)]';

% disp([mu_0,mean(Simu_mu),std(Simu_mu)]);
% disp([kappa2_0,mean(Simu_kappa2),std(Simu_kappa2)]);
% disp([eta_0,mean(Simu_eta),std(Simu_eta)]);
% disp([zeta_0,mean(Simu_zeta),std(Simu_zeta)]);