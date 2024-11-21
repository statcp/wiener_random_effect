function estP=WienerNIG(Obs)
%obs={{T1,X1},{T2,X2},...{Tn,Xn}}.
myoptimset=optimset('display','off');
funcLam=@(t,alpha) t.^alpha;
alpha_0=1;
% n=length(Obs); % number of groups: all the units, even in the same group, are i.i.d
% % the units in the same group has the same observation times
tmpfunc=@(a) myfunc(Obs,funcLam,a,1);
alpha1=fmincon(tmpfunc,alpha_0,[],[],[],[],0,10,[],myoptimset);
estP=myfunc(Obs,funcLam,alpha1);
estP.alpha=alpha1;

    function y=myfunc(obs,funclam,a,opt)
        n=length(obs);
        new_obs=cell(1,n);
        for i=1:n
            new_obs{i}={funclam(obs{i}{1},a),obs{i}{2}};
        end
        tmp_res=WienerNIG_L(new_obs);
        if nargin>=4
            if opt==1
                y=-tmp_res.loglike;
            else
                y=tmp_res;
            end
        else
            y=tmp_res;
        end
    end
end