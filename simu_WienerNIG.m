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
Simu_dt=zeros(Nsim,9);
%
nboot=5e3;
Simu_boot_mu=cell(Nsim,9);
Simu_boot_kappa2=cell(Nsim,9);
Simu_boot_eta=cell(Nsim,9);
Simu_boot_zeta=cell(Nsim,9);
% percentile coverage
Simu_boot_mu_CP1=zeros(Nsim,9);
Simu_boot_kappa2_CP1=zeros(Nsim,9);
Simu_boot_eta_CP1=zeros(Nsim,9);
Simu_boot_zeta_CP1=zeros(Nsim,9);
%
count_nm=0;
for n1=1:3
    for m1=1:3
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
            %
            Simu_mu(rep,count_nm)=estP_test.mu;
            Simu_kappa2(rep,count_nm)=estP_test.kappa2;
            Simu_eta(rep,count_nm)=estP_test.eta;
            Simu_zeta(rep,count_nm)=estP_test.zeta;
            estP_test.T=T_cell;
            %bootstrapped samples
            bootstat_test=bootstrapDM(nboot,@WienerNIG_L_sam,estP_test,@WienerNIG_L);
            tmp_boot_mu=zeros(nboot,1);
            tmp_boot_kappa2=zeros(nboot,1);
            tmp_boot_eta=zeros(nboot,1);
            tmp_boot_zeta=zeros(nboot,1);
            for kk=1:nboot
                tmp_boot_mu(kk)=bootstat_test{kk}.mu;
                tmp_boot_kappa2(kk)=bootstat_test{kk}.kappa2;
                tmp_boot_eta(kk)=bootstat_test{kk}.eta;
                tmp_boot_zeta(kk)=bootstat_test{kk}.zeta;
            end
            Simu_boot_mu{rep,count_nm}=tmp_boot_mu;
            Simu_boot_kappa2{rep,count_nm}=tmp_boot_kappa2;
            Simu_boot_eta{rep,count_nm}=tmp_boot_eta;
            Simu_boot_zeta{rep,count_nm}=tmp_boot_zeta;
            % cp
            %interval estimation by percentile
            boot_mu_sort=sort(tmp_boot_mu);
            if boot_mu_sort(round(nboot*(1-p_cover)/2))<=mu_0&&...
                    boot_mu_sort(round(nboot*(1+p_cover)/2))>=mu_0
                Simu_boot_mu_CP1(rep,count_nm)=1;
            end
            boot_kappa2_sort=sort(tmp_boot_kappa2);
            if boot_kappa2_sort(round(nboot*(1-p_cover)/2))<=kappa2_0&&...
                    boot_kappa2_sort(round(nboot*(1+p_cover)/2))>=kappa2_0
                Simu_boot_kappa2_CP1(rep,count_nm)=1;
            end
            boot_eta_sort=sort(tmp_boot_eta);
            if boot_eta_sort(round(nboot*(1-p_cover)/2))<=eta_0&&...
                    boot_eta_sort(round(nboot*(1+p_cover)/2))>=eta_0
                Simu_boot_eta_CP1(rep,count_nm)=1;
            end
            boot_zeta_sort=sort(tmp_boot_zeta);
            if boot_zeta_sort(round(nboot*(1-p_cover)/2))<=zeta_0&&...
                    boot_zeta_sort(round(nboot*(1+p_cover)/2))>=zeta_0
                Simu_boot_zeta_CP1(rep,count_nm)=1;
            end
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
CP_perc=[sum(Simu_boot_mu_CP1);sum(Simu_boot_kappa2_CP1);sum(Simu_boot_eta_CP1);sum(Simu_boot_zeta_CP1)]/Nsim;