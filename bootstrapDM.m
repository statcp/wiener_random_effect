function [bootstat,bootsam]=bootstrapDM(nboot,bootfun,para0,estfun,paraflag,npara)
% This is the parametric bootstrap method to obtain resamples for the
% estimates of a degradation model
% NBOOT: is the number of resamples;
% BOOTFUN: is the function to generate resamples for the degradation model
% PARA0: is the model parameters (point estimates) for the degradation
% model that will be used to generate resamples
% ESTFUN: is the function used to obtain resamples for the point estimates
% by estimating the model parameters using ESTFUN with the resamples for
% degradation data as input
% PARAFLAG: parallel bootstrap is used if PARAFLAG is TRUE. The default is
% FALSE.
% NPARA: is the number of parallel cores in parallel bootstrap
if nargin<5
    paraflag=false;
end
if nargout>1
    sam_flag=true;
else
    sam_flag=false;
end
bootstat=cell(nboot,1);
bootsam=cell(nboot,1);
if paraflag % parallel bootstrap
    try
         parpool(npara);
    end
    fprintf(repmat('.',1,nboot));
    fprintf('\n\r');
    parfor rep=1:nboot
        obs=bootfun(para0);
        bootest=estfun(obs);
        bootstat{rep}=bootest;
        if sam_flag
            bootsam{rep}=obs;
        end
        fprintf('\b|\n');
    end
else
    for rep=1:nboot
        obs=bootfun(para0);
        bootest=estfun(obs);
        bootstat{rep}=bootest;
        if sam_flag
            bootsam{rep}=obs;
        end
    end
end