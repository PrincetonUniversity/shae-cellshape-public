function [E,Elength]=EdgeLength(F)

%takes patch object F and returns indices of edge end points and length of edges

E=[F.faces(:,[1,3]);F.faces(:,2:3);F.faces(:,1:2)];
E=sort(E,2);
[E,~,~]=unique(E,'rows');


if nargout==2
E1=F.vertices(E(:,1),:);
E2=F.vertices(E(:,2),:);
Elength=E1-E2;
Elength=sqrt(sum(Elength.^2,2));
end