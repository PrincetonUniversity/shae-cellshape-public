function [T,C]=Twist(V,tol)
% this function take the centerline, V of a cell and a smoothing length P
% and calculated a 2D twist
% 
% the input 'tol' specifies a tolerance for how far the B-spline (quintic) 
% can vary from the input data 

if nargin<2
    % use the rms deviation from one centerline point to the next as the
    % tolerance for the B-splines
tol = double(nanmean(sqrt(sum((diff(V).^2),2))))/3;
% disp(tol)
    
end

 if length(V)<10
     T=[];
 else
N=length(V);

% a simple linear extrapolation seems really miss what's happening at the ends
V2=zeros(N,3);
% establish the linear variable that defines contour length
s = 1:size(V,1);
% recast as a double
V = double(V);

% fit each dimension independently as a smoothing spline, then use this
% smoothing spline to evaluate the derivatives
[x1] = spaps(s,V(:,1),tol,[],3);
[x2] = spaps(s,V(:,2),tol,[],3);
[x3] = spaps(s,V(:,3),tol,[],3);
x1 = fnxtr(x1,3);
x2 = fnxtr(x2,3);
x3 = fnxtr(x3,3);

% evaluate the extrapolation 
% s2 = 1-2*p:N+2*p;
s2 = 1:N;
V2(:,1)=fnval(x1,s2)';
V2(:,2)=fnval(x2,s2)';
V2(:,3)=fnval(x3,s2)';

% figure(); scatter3(V(:,1),V(:,2),V(:,3),'rx'); hold on; scatter3(V2(:,1),V2(:,2),V2(:,3))


%%

% tau, the torsion is 
% |xdot xdoubledot xtripledot|/|xdot xcross xdoubledot|^2
% which is equal to rho^2 * |xdot xdoubledot xtripledot| where rho is
% radius of curvature

% first derivative
x1dot = fnval(fnder(x1,1),s2);
x2dot = fnval(fnder(x2,1),s2);
x3dot = fnval(fnder(x3,1),s2);
% concatenate together
bigXdot = [x1dot',x2dot',x3dot'];

% second derivative
x1doubledot = fnval(fnder(x1,2),s2);
x2doubledot = fnval(fnder(x2,2),s2);
x3doubledot = fnval(fnder(x3,2),s2);
% concatenate together
bigXdoubledot = [x1doubledot',x2doubledot',x3doubledot'];

% third derivative
x1tripledot = fnval(fnder(x1,3),s2);
x2tripledot = fnval(fnder(x2,3),s2);
x3tripledot = fnval(fnder(x3,3),s2);
% concatenate together
bigXtripledot = [x1tripledot',x2tripledot',x3tripledot'];

C1 = dot(cross(bigXdot',bigXdoubledot'),bigXtripledot')';
C2 = dot(cross(bigXdot',bigXdoubledot'),cross(bigXdot',bigXdoubledot'))';

C3 = cross(bigXdot',bigXdoubledot');
C3 = sqrt(dot(C3,C3));
C4 = sqrt(dot(bigXdot',bigXdot')).^3;
C = (C3./C4);
% C = medfilt1(C,p);
% T = medfilt1(C1./C2,p);
T = C1./C2;


% % % %%
% % % 
% % % [~,dS]=gradient(V,p);
% % % dS=dS;
% % % dS=medfilt2(dS,[5,1]);
% % % [~,ddS]=gradient(dS,p);
% % % ddS=ddS;
% % % 
% % % [~,dddS]=gradient(ddS,p);
% % % dddS=dddS;
% % 
% % 
% % 
% % 
% % 
% % C=cross(dS,ddS);
% % 
% % 
% % T=dot(C',dddS')./dot(C',C');
% the ends are poorly defined, set them to nan. For "average" torsion
% calculations, there should probably be even more removed
T = T(7:end-6)';
T = cat(1,nan(6,1),T',nan(6,1))';
C = cat(1,nan(6,1),C(7:end-6)',nan(6,1));
 end

