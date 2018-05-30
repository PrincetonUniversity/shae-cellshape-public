function [t_vector,b_vector,n_vector]=tbnVector(xyzs)

%takes an N by 3 matrix representing a set of coordinates in xyz and
%calculates the tangent, normal, and binormal vectors at each point. not
%that T cross N = B;

t_vector=diff(xyzs);

    t_vector=normalizeRows(t_vector);
if size(t_vector,1)>1 && sum(sum(abs(diff(t_vector))))>.0001
b_vector=cross(t_vector,[[0,0,0]',diff(t_vector)']');
else 
    b_vector=cross(t_vector,repmat([0,0,1],size(t_vector,1),1),2);
    b_vector=b_vector/norm(b_vector);
end

if all(all(b_vector==0,2))
    b_vector(1,:)=[0,0,1];
end

for i=find(~all(b_vector==0,2),1,'last'):size(b_vector,1)
    if all(b_vector(i,:)==0)&& i~=1
        b_vector(i,:)=b_vector(i-1,:);
    end   
end
first_blank=find(all(b_vector==0,2),1,'last');
for i=first_blank:-1:1
    if all(b_vector(i,:)==0)
        b_vector(i,:)=b_vector(i+1,:);
    end
    
end

for i=1:size(b_vector,1)
    b_vector(i,:)=b_vector(i,:)/norm(b_vector(i,:));
end
n_vector= cross(t_vector,b_vector,2);


for i=3:size(t_vector,1)
    newN=cross(t_vector(i-1,:),b_vector(i,:));
    theta=subspace(newN',n_vector(i-1,:)');
    
    sinsign=sign((dot(cross(newN',n_vector(i-1,:)'),t_vector(i-1,:))));
    cossign=sign(dot(newN',n_vector(i-1,:)'));
    

    if sinsign>0 && cossign<0
            theta=pi-theta;
    end
    if sinsign<=0 && cossign<=0
            theta=pi+theta;
    end
    if sinsign<0 && cossign>0
            theta=2*pi-theta;
    end
    
	b_vector(i,:) =  rodrigues_rot(b_vector(i,:),t_vector(i,:),theta);
    n_vector(i,:) =  rodrigues_rot(n_vector(i,:),t_vector(i,:),theta);
    
end
b_vector=b_vector;
