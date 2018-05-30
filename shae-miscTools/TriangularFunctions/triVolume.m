function V= triVolume(F)

%Calculates volume of simple triangular mesh given by patch object F.
%Object cannot fold in on itself yet. I need to figure out how to implement
%the sign convention for this. For now, all i do is add up the volume of
%all tetrahedrons given by a single face and the center of mass. 
F.vertices=bsxfun(@minus,F.vertices,mean(F.vertices));
V=dot(cross(F.vertices(F.faces(:,1),:),F.vertices(F.faces(:,2),:)),...
    F.vertices(F.faces(:,3),:));
V=sum(abs(V))/6;