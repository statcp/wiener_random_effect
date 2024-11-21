%% Ratio of Modified Bessel funtion of the Second Kind
% The ratio between the modified Bessel functions of the second kind:
% besselk(p,x)/besselk(p-1,k)
%%
function tout=Rbesselk(p,x)
% tout=ones(size(x));
% tmpI=(x~=Inf);
% tmpx=x(tmpI);
% tout(tmpI)=besselmx(double('K'),p,tmpx,1)./besselmx(double('K'),p-1,tmpx,1);
    %the last argument is the scale function, 
    %if scale=0, then besselmx returns Kp(x)
    %if scale=1, then besselmx returns Kp(x)/e^(-x)
if p<0
    p=abs(p)+1;
    tout=1./Rbesselk(p,x);
else if p>2
% else if isnan(besselk(p,x,1)./besselk(p-1,x,1))
        tout=2*(p-1)./x+1./Rbesselk((p-1),x);
    else
        if besselk(p,x,1)<=0||x==Inf
            tout=ones(size(x));
        else
            tout=besselk(p,x,1)./besselk(p-1,x,1);
            %the last argument is the scale function, 
            %if scale=0, then besselmx returns Kp(x)
            %if scale=1, then besselmx returns Kp(x)/e^(-x)
        end
    end
end