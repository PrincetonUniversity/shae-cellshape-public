function K=Curvature(V,p)
% this function take the centerline, V of a cell and a smoothing length P
% and calculated a 2D curvature

 if length(V)<10
     K=[];
 else
N=length(V);
if p>N/4
    p=round(N/4);
end
if size(V,2)==1
    
end

dS=zeros(size(V));
ddS=dS;

    for i=p+1:N-p
    dS(i,:)=(V(i+p,:)-V(i-p,:))/(2*p);
    ddS(i,:)=(V(i+p,:)+V(i-p,:)-2*V(i,:))/p^2;
    end
    for i=5:p

        
    dS(i,:)=(V(2*i-1,:)-V(1,:))/(2*(i-1));
    ddS(i,:)=(V(2*i-1,:)+V(1,:)-2*V(i,:))/(i-1)^2;
    end

        
    for i=N-p:N-5
    dS(i,:)=(V(N,:)-V(2*i-N,:))/(2*(N-i));
    ddS(i,:)=(V(N,:)+V(2*i-N,:)-2*V(i,:))/(N-i)^2;
    end
    for i=1:5
        dS(i,:)=dS(5,:);
        ddS(i,:)=ddS(5,:);
        dS(N-i+1,:)=dS(N-5,:);
        ddS(N-i+1,:)=ddS(N-5,:);
    end


K=zeros(length(V),1);
for i=1:length(V)
    K(i)=norm(cross(dS(i,:),ddS(i,:)))/norm(dS(i,:))^3;
    if isnan(K(i))
        K(i)=0;
    end
    
end
 end