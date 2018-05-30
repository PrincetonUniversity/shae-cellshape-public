function [FL,FU]=Finternal3D(x,y,z,cline_para,weight)
if ~exist('weight','var')
    weight='uniform';
end

K=numel(x);
index=1:K;
positionMap=zeros(size(x,1),size(x,2));
positionMap(ind2sub(size(x),index))=index;
positionMap2=cat(3, circshift(positionMap,[1,0]),circshift(positionMap,[-1,0]),...
    circshift(positionMap,[0,1]),circshift(positionMap,[0,-1]));
circscale=.5;

for i=1:4
    positionMap2(:,1,1)=circshift(positionMap(:,2),round(size(positionMap,1)/4));
    positionMap2(:,1,2)=circshift(positionMap(:,2),round(3*size(positionMap,1)/4));
    positionMap2(:,1,3)=circshift(positionMap(:,2),round(size(positionMap,1)/2));
    
    positionMap2(:,end,1)=circshift(positionMap(:,end-1),round(size(positionMap,1)/4));
    positionMap2(:,end,2)=circshift(positionMap(:,end-1),round(3*size(positionMap,1)/4));
    positionMap2(:,end,4)=circshift(positionMap(:,end-1),round(size(positionMap,1)/2));
    
end

    
positionMap2=reshape(positionMap2,K,4);
%%
nMatrix=spalloc(K,K,K*5);
for i=1:K
    
    if strcmp(weight,'distance')
%         if i<size(x,1)
%             
%             
%             nMatrix(i,size(x,1)+1:2*size(x,1))=1/size(x,1);
%         elseif i>=K-size(x,1)+1
%             
%           nMatrix(i,end-2*size(x,1)+1:end-size(x,1))=1/size(x,1);
%      else
if i<size(x,1)
    ind=positionMap2(1:size(x,1),:);
ind=unique(ind(:))';
    neighpos=[x(ind)',y(ind)',z(ind)'];
    
elseif i>=K-size(x,1)+1
    ind=positionMap2(K-size(x,1)+1:K,:);
ind=unique(ind(:))';    
    neighpos=[x(ind)',y(ind)',z(ind)'];
    
else
    ind=positionMap2(i,:);
    neighpos=[x(ind)',y(ind)',z(ind)'];
    
end

%neighpos=[x(positionMap2(i,:))',y(positionMap2(i,:))',z(positionMap2(i,:))'];
           pos=[x(i),y(i),z(i)];
           ds=bsxfun(@minus,neighpos,pos);
           ds=dot(normalizeRows(ds),ds,2);
           ds=ds.^-1*4/length(ind);
           ds(isnan(ds)|isinf(ds))=0;
           w=ds./sum(ds);
nMatrix(i,ind)=ds; % 1 2 are neighbours along the axis.
nMatrix(i,i)=-sum(ds);
  %3,4 are neighbours in a single slice
   %      end
         else
         if i<size(x,1) 
            nMatrix(i,size(x,1)+1:2*size(x,1))=1/size(x,1);
     elseif i>=K-size(x,1)+1
        
          nMatrix(i,end-2*size(x,1)+1:end-size(x,1))=1/size(x,1);
     else 
nMatrix(i,squeeze(positionMap2(i,[1,2])))=1/2*circscale; % 1 2 are neighbours along the axis.
nMatrix(i,squeeze(positionMap2(i,[3,4])))=1/2*(1-circscale);  %3,4 are neighbours in a single slice
         end
             nMatrix(i,i)=-1;

    end
end
nMatrix=nMatrix;
% Cvec=nMatrix*[x(:),y(:),z(:)];
% C=dot(normr(Cvec),Cvec,2);
% C=reshape(C,size(x));
% d=reshape(d,size(x));
% c=reshape(c,size(x));
%%
f=(cline_para.gamma*speye(K,K)-cline_para.zalpha3d*nMatrix+cline_para.zbeta3d*nMatrix*nMatrix);
[FL,FU]=lu(f);
%force=inv(cline_para.gamma*speye(K,K)-cline_para.alpha*nMatrix+cline_para.beta*nMatrix*nMatrix);
% 
% [xu,xv]=gradient(x);
% [yu,yv]=gradient(y);
% [zu,zv]=gradient(z);
% 
% [xuu,xuv]=gradient(xu);
% [xvu,xvv]=gradient(xv);
% [yuu,yuv]=gradient(yu);
% [yvu,yvv]=gradient(yv);
% [zuu,zuv]=gradient(zu);
% [zvu,zvv]=gradient(zv);
