function V=surfVolume(x,y,z)

%Takes a closed contour of the form S=(x(u,v),y(u,v),z(u,v)) and calculates
%the Volume. 
cl=[mean(x)',mean(y)',mean(z)'];
    at=zeros(size(x,2),1);
    for i=1:size(x,2)
        rt=zeros(size(x,1),1);
        for j=1:size(x,1)
            rt(j)=mean(sqrt((x(j,i)-mean(x(:,i))).^2+(y(j,i)-mean(y(:,i))).^2+(z(j,i)-mean(z(:,i))).^2));
        end
        at(i)=sum(rt.*circshift(rt,1))*pi/ cline_para.pts;
    end
    
    ds=sqrt(sum(diff(cl).^2,2));
    at=(at+circshift(at,1))/2;
    at(end)=[];
    V=sum(ds.*at);