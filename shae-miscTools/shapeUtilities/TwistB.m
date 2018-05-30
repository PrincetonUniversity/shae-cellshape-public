function T=TwistB(V,p)
% this function take the centerline, V of a cell and a smoothing length P
% and calculated a 2D twist

 if length(V)<10
     T=[];
 else
N=length(V);
if p>N/4
    p=round(N/4);
end
V2=zeros(N+4*p,3);
for i=1:3
    V2(:,i)=interp1(1:N,V(:,i),1-2*p:N+2*p,'linear','extrap')';
end
V=V2;
if size(V,2)==1
    
end


[~,dS]=gradient(V,p);
dS=dS;
dS=medfilt2(dS,[5,1]);
[~,ddS]=gradient(dS,p);
ddS=ddS;

[~,dddS]=gradient(ddS,p);
dddS=dddS;





C=cross(dS,ddS);


T=dot(C',dddS')./dot(C',C');
T=T(2*p+1:end-p*2);
 end

