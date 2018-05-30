function K=Curvature(V,p)
% this function take the centerline, V of a cell and a smoothing length P
% and calculated a 2D curvature

% check to make sure the length of the curve has sufficient number of
% points
if length(V)<10
    K=[];
else
    N=length(V);
    % rounding length
    if p>N/4
        p=round(N/4);
    end
   
    % preallocate derivative vectors
    dS=zeros(size(V));
    ddS=dS;
    
    % estimate discrete derivatives taken at lengthscale p
    for i=p+1:N-p
        dS(i,:)=(V(i+p,:)-V(i-p,:))/(2*p);
        ddS(i,:)=(V(i+p,:)+V(i-p,:)-2*V(i,:))/p^2;
    end
    
    
    %  estimate discrete derivatives near the start of the curve
    for i=5:p
        dS(i,:)=(V(2*i-1,:)-V(1,:))/(2*(i-1));
        ddS(i,:)=(V(2*i-1,:)+V(1,:)-2*V(i,:))/(i-1)^2;
    end
    
    %  estimate discrete derivatives near the end of the curve
    for i=N-p:N-5
        dS(i,:)=(V(N,:)-V(2*i-N,:))/(2*(N-i));
        ddS(i,:)=(V(N,:)+V(2*i-N,:)-2*V(i,:))/(N-i)^2;
    end
    
    % estimate discrete derivatives at the start and end of the curve
    for i=1:5
        dS(i,:)=dS(5,:);
        ddS(i,:)=ddS(5,:);
        dS(N-i+1,:)=dS(N-5,:);
        ddS(N-i+1,:)=ddS(N-5,:);
    end
    
    % pre-allocate K to hold calculated curvatures
    K = zeros(length(V),1);
    for i=1:length(V)
        % at each point along the curve calculate the norm of the
        % ||s' X s''||/||s'||^3
        K(i)=norm(cross(dS(i,:),ddS(i,:)))/norm(dS(i,:))^3;
        if isnan(K(i))
            K(i)=0;
        end
        
    end
end