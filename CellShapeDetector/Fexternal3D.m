function [Fx,Fy,Fz]=Fexternal3D(V)



            x = -12:12;
            h = exp(-x.^2/2/.9^2); 
            h = h/sum(h);
            V = convn(V, h, 'same');
            V = convn(V, h', 'same');
            V=convn(V,reshape(h,[1,1,length(h)]),'same');
%             h = exp(-x.^2/2/1^2); 
%             h = h/sum(h);
%             V = convn(V, reshape(h, 1,1,length(h)), 'same');
[Fy,Fx,Fz]=gradient(V);

