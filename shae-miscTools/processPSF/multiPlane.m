function [out,bSave] = multiPlane(a,x,zList,zPrime,imgSize);

% x(:,1) and x(:,2) are the paired xy coordinates to evaluate
% zList is a list of the zPositions in a single stack
% zPrime is a single value indicating where the sample is relative to the
% surface (in um of stage motion? BPB 20120910 1055 am)
if nargout>1
    bSave =[];
end

c = a;
% linear in zPrime
for ii = 1:8
    a(ii) = c(ii) + c(ii+8).*zPrime;
end

for iZ = 1:size(zList,2);
% (b1,b2,b3,b4,b5,b6,b7,b10,b11)

    % b(1) = amplitude;
    % gaussian with offset plus linear
    b(1) = a(1) + a(2).*exp(-((zList(iZ)).^2/(2*a(3)^2)))+a(4).*zList(iZ);
    % b(2) = width of primary gaussian
    b(2) = polyval(a([5,6,7,8]),zList(iZ));
    
    b(isnan(b)) = 0;
    
    % store coeffs
    if nargout>1
        bSave = cat(1,bSave,b);
    end
    out(:,:,iZ) = reshape(singlePlane(b,x(:,1:2)),imgSize);
end
