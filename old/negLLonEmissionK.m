function [Q,dQ,d2Q]=negLLonEmissionK(k,w,x,y,dt,model)

    %
    % function [Q,dQ,d2Q]=negLLonEmissionK(k,w,x,y,dt,model)
    %
    %  INPUTS:
    %  k           The receptive field to be optimized.
    %  w           The weights for being in the relevant state.
    %  x           The stimulus as an NxT matrix, rows are pixels, columns
    %              are time points.  To allow for a bias term in the k12
    %              and k21 vectors, the final row of x should be all 1's
    %  y           The spike-train.
    %  dt          The time-step in seconds.
    %  model       If model is 'Bernoulli', then a Bernoulli spiking model
    %              is used.  Otherwise, the default Poisson spiking model
    %              is used.
    %
    %  OUTPUTS:
    %  Q           The negative log-likelihood.
    %  dQ          The gradient of the negative log-likelihood.
    %  d2Q         The Hessian of the negative log-likelihood.
    %

    if nargin==6 && ischar(model) && strcmpi(model,'bernoulli')
        ldt=log(dt);
        d=length(k);
        x=x';
        kx=k'*x;
        silent=~y;
        pos=kx>0;
        pos0=pos&silent;
        pos1=pos&y;
        neg=~pos;
        neg0=neg&silent;
        neg1=neg&y;
        
        wfdt_n0=w(neg0).*exp(kx(neg0)+ldt);
        
        wdt_p0=w(pos0)*dt;
        oneplus_p0=1+kx(pos0);
        
        fdt_n1=exp(kx(neg1)+ldt);
        expfdt_n1=exp(fdt_n1);
        oneminus_n1=1-1./expfdt_n1;
        
        oneplus_p1=1+kx(pos1);
        expfdt_p1=exp((oneplus_p1+.5*kx(pos1).^2)*dt);
        oneminus_p1=1-1./expfdt_p1;
        
        Q=sum(wfdt_n0)+sum(wdt_p0.*(oneplus_p0+.5*kx(pos0).^2))-sum(w(neg1).*log(oneminus_n1))-sum(w(pos1).*log(oneminus_p1));

        if nargout > 1
            tmp_n1=w(neg1).*fdt_n1./(expfdt_n1-1);
            tmp_p1=w(pos1)*dt./(expfdt_p1-1);
            dQ=x(:,neg0)*wfdt_n0'+x(:,pos0)*(wdt_p0.*oneplus_p0)'-x(:,neg1)*tmp_n1'-x(:,pos1)*(tmp_p1.*oneplus_p1)';
            if nargout > 2
                d2Q=(repmat(wfdt_n0,d,1).*x(:,neg0))*x(:,neg0)';
                d2Q=d2Q+(repmat(wdt_p0,d,1).*x(:,pos0))*x(:,pos0)';
                d2Q=d2Q-(repmat(tmp_n1.*(1-fdt_n1./oneminus_n1),d,1).*x(:,neg1))*x(:,neg1)';
                d2Q=d2Q-(repmat(tmp_p1.*(1-oneplus_p1.^2*dt./oneminus_p1),d,1).*x(:,pos1))*x(:,pos1)';
            end
        end       
    else
        d=length(k);
        x=x';
        kx=k'*x;
        index1=kx<=0;
        x1=x(:,index1);
        kx1=kx(index1);
        w1=w(index1);
        temppy1=-w1.*y(index1);
        temppexp=w1.*exp(kx1)*dt;

        index2=~index1;
        x2=x(:,index2);
        kx2=kx(index2);
        w2=w(index2);
        temppy2=w2.*y(index2);
        temppdt=w2*dt;
        onepluskx=1+kx2;
        lambda=onepluskx+.5*kx2.^2;

        Q=sum(temppexp+temppy1.*kx1) + sum(lambda.*temppdt-temppy2.*log(lambda));

        if nargout > 1
            yoverlambda=temppy2./lambda;
            bigterm=temppdt-yoverlambda;
            dQ=x1*(temppexp+temppy1)' + x2*(bigterm.*onepluskx)';
            if nargout > 2
                d2Q=(repmat(temppexp,d,1).*x1)*x1';
                d2Q=d2Q + (repmat(bigterm+yoverlambda.*onepluskx.^2./lambda,d,1).*x2)*x2';
            end
        end
    end
end