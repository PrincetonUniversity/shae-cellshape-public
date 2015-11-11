function y = qmean(x,quant)

%qmean  takes a matrix or vector x and returns
%the mean of some upper quantile of values 
q=quantile(x,quant);
for icol=1:size(x,2);
    xcol=x(:,icol);
    if all(xcol==q(icol))
        q(icol)=-Inf;
    end
    y(icol)=nanmean(xcol(xcol>q(icol)));
end

    