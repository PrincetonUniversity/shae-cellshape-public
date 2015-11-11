function varargout=  stackLoad(stack_name, stack_z_size,Z_scale)
% takes a tiff stack and loads them into 3D matrices. Will take the entire
% stack and load them in chuncks of stack_z_size, if the stack is larger.
% Use multiple outputs if you expect multiple stacks because of this.
% default stack_z_size is size of entire stack. If there are no outputs,
% will display matrix with sliceBrowser.


%metadata patch
try
    dashPos = find(double(stack_name)==45,1,'last');
    metaFileName = stack_name;
    metaFileName(dashPos:end) = [];
    metaFileName = [metaFileName,'-Metadata.mat'];
    load(metaFileName);
    isMetaData = 1;
catch ME
    isMetaData = 0;
end



[row, col] = size(imread(stack_name, 'tif', 1));
image_info=imfinfo(stack_name);
totalStack=length(image_info);
if nargin==1;
    if nargout
    stack_z_size =totalStack/nargout;
    else
        stack_z_size=totalStack;
    end
    Z_scale=1;         %.6   %scaling between distance in the slide to distance moved in the image
end



if nargin==2
    Z_scale=1;         %.6   %scaling between distance in the slide to distance moved in the image
    if isempty(stack_z_size)
        stack_z_size =totalStack;
    end
end

if isempty(stack_z_size)
    stack_z_size=totalStack/nargout;
end

for h=1:min(floor(totalStack/stack_z_size))
    hyper_stack=zeros(row,col,stack_z_size);
    for k = 1:stack_z_size
        hyper_stack(:,:,k) = double(imread(stack_name, 'tif',...
            (h-1)*stack_z_size+k,'info',image_info));
    end
    
    hyper_stack=image_resize(hyper_stack,size(hyper_stack,1),...
        size(hyper_stack,2),round(Z_scale*size(hyper_stack,3)));
    
    if isMetaData
        % background subtract, shiftImgFilter, background add
        backgndTmp = quantile(hyper_stack(:),0.005);
        hyper_stack = hyper_stack-backgndTmp;
        for k2 = 1:image_para.stack_z_size
            if flag.reverse==1 && flag.F1==1
                planeIndex = image_para.Fstack_t_size*image_para.Fstack_z_size+k2;
            else
                planeIndex = k2;
            end
            
            filterString = totalImageInfo{planeIndex}.filterID;
            filterIdx = str2double(filterString(end-1));
            if isnan(filterIdx)
                filterIdx = 3; % no shifting of image
            end
            
            hyper_stack(:,:,k2) = shiftImgFilter(hyper_stack(:,:,k2),filterIdx);
            
        end
        hyper_stack = hyper_stack+backgndTmp;
    end
    
    
    if nargout==0
        SliceBrowser(hyper_stack)
    end
    varargout{h}=hyper_stack;
    
end
