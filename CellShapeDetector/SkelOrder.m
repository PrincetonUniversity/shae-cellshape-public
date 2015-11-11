function pt_order=  SkelOrder(bw)
if size(size(bw),2)==2
    for ii=1:2
        bw2=bw;
        pts=find(bw2==1);
        [pts(:,1),pts(:,2)]=ind2sub(size(bw2),pts); 
        pt_order=zeros(size(pts));
        %pt_order(1,:)=pts(1,:);

        for i=1:length(pts)

            if i==1
            if ii==1
                pt_order(1,:)=pts(1,:);
                bw2(pt_order(1,1),pt_order(1,2))=0;
                if size(pts,1)==1
      last_pt=pts(1,:);
                    break
                end
                
            else if ii==2
                    pt_order(1,:)=last_pt;
                    bw2(pt_order(1,1),pt_order(1,2))=0;
                end
            end
            end
            l_edge=pt_order(i,1)-1;
            r_edge=pt_order(i,1)+1;
            t_edge=pt_order(i,2)+1;
            b_edge=pt_order(i,2)-1;
            
            if all(~bw2(l_edge:r_edge,b_edge:t_edge))
       
                break
            end
            
            [next_pt(1),next_pt(2)]=find(bw2(l_edge:r_edge,b_edge:t_edge)==1,1,'first');
            
            next_pt=next_pt-[2,2];
            pt_order(i+1,:)=next_pt+pt_order(i,:);
            last_pt=next_pt+pt_order(i,:);
            bw2(pt_order(i+1,1),pt_order(i+1,2))=0;
            %imagesc(bw2);figure(gcf);

        end
    end
end

if size(size(bw),2)==3
    for ii=1:2
        bw2=bw;
        pts=find(bw2==1);
        [pts(:,1),pts(:,2),pts(:,3)]=ind2sub(size(bw2),pts);
        
        pt_order=zeros(size(pts));
        
        if ii==1
            pt_order(1,:)=pts(1,:);
            bw2(pt_order(1,1),pt_order(1,2),pt_order(1,3))=0;
            if size(pts,1)==1
                break
            end
            
        else if ii==2
                pt_order(1,:)=last_pt;
                bw2(pt_order(1,1),pt_order(1,2),pt_order(1,3))=0;
            end
        end
        
        
        for i=1:length(pts)
            
            l_edge=pt_order(i,1)-1;
            r_edge=pt_order(i,1)+1;
            t_edge=pt_order(i,2)+1;
            b_edge=pt_order(i,2)-1;
            i_edge=pt_order(i,3)-1;
            o_edge=pt_order(i,3)+1;
            if all(~bw2(l_edge:r_edge,b_edge:t_edge,i_edge:o_edge))
                break
            end
            [next_pt(1),next_pt(2),next_pt(3)]=ind2sub([3,3,3],find(bw2(l_edge:r_edge,b_edge:t_edge,i_edge:o_edge)==1,1,'first'));
            
            next_pt=next_pt-[2,2,2];
            pt_order(i+1,:)=next_pt+pt_order(i,:);
            last_pt=next_pt+pt_order(i,:);
            bw2(pt_order(i+1,1),pt_order(i+1,2),pt_order(i+1,3))=0;
        end
        
    end
end

