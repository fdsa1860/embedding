% weizmann load data
% manually delete 
clear
load ('weizzmann_classification_masks.mat','original_masks');
field_name = fieldnames(original_masks);


cnt_skip=1;
cnt_run=1;
for i =1: length(field_name)
    if strcmp(field_name{i}(end-4:end),'walk1')  || strcmp(field_name{i}(end-3:end),'walk')
        skip_mask{cnt_skip} = getfield(original_masks, field_name{i});
        cnt_skip = cnt_skip+1;
    elseif strcmp(field_name{i}(end-4:end),'walk2') % flip left and right of the frames.
        temp = getfield(original_masks, field_name{i});
        for k =1:size(temp,3)
            temp(:,:,k) = fliplr(temp(:,:,k));
        end
        skip_mask{cnt_skip} = temp;
        cnt_skip = cnt_skip+1;
    end
        
    if strcmp(field_name{i}(end-3:end),'run2') || strcmp(field_name{i}(end-2:end),'run')
        run_mask{cnt_run} = getfield(original_masks, field_name{i});
        cnt_run = cnt_run+1;
    elseif strcmp(field_name{i}(end-3:end),'run1') % flip left and right of the frames.
        temp = getfield(original_masks, field_name{i});
        for k =1:size(temp,3)
            temp(:,:,k) = fliplr(temp(:,:,k));
        end
        run_mask{cnt_skip} = temp;
        cnt_run = cnt_run+1;
    end
end
cnt_skip = cnt_skip - 1;
cnt_run = cnt_run - 1;
% count the number of total frames
num_frame =0;
for i =1 : cnt_skip
    num_frame =num_frame + size(skip_mask{i},3);
end
for i =1 : cnt_run
    num_frame =num_frame + size(run_mask{i},3);
end
I = zeros(5000, num_frame);
class = zeros( num_frame, 1);
% get the boundingbox scale the boundingbox
flag_skip_flip = [1 1 0 0 1 0 0 1 0 0];
flag_run_flip = [1 1 1 0 1 1 1 1 1 0];
num_frame =1;
for i =1 : cnt_skip
    for k =1:size(skip_mask{i},3)
        Itemp = im2bw(squeeze(skip_mask{i}(:,:,k)));
        stat =regionprops(Itemp,'BoundingBox');
        rect_temp = stat(1).BoundingBox;
        Itemp = Itemp(rect_temp(2):rect_temp(4)+rect_temp(2)-1, rect_temp(1):rect_temp(1)+rect_temp(3)-1 );
        Itemp = imresize(Itemp, [100, 50]);
        if ~flag_skip_flip(i)
            Itemp = fliplr(Itemp);
        end
        I(:,num_frame)= Itemp(:);
        class(num_frame, 1) = 1;
        num_frame = num_frame+1;
        imshow(Itemp,[]);
        drawnow;
    end
%     pause;
end
for i =1 : cnt_run
    for k =1:size(run_mask{i},3)
        Itemp = im2bw(squeeze(run_mask{i}(:,:,k)));
        stat =regionprops(Itemp,'BoundingBox');
        rect_temp = stat(1).BoundingBox;
        Itemp = Itemp(rect_temp(2):rect_temp(4)+rect_temp(2)-1, rect_temp(1):rect_temp(1)+rect_temp(3)-1 );
        Itemp = imresize(Itemp, [100, 50]);
        if ~flag_run_flip(i)
            Itemp = fliplr(Itemp);
        end
        I(:,num_frame)= Itemp(:);
        class(num_frame, 1) = 2;
        num_frame = num_frame+1;
        imshow(Itemp,[]);
        drawnow;
    end
%     pause;
end
%%
save('weizmann_run_walk_all','I','class');
I = I(:,1:3:end);
class = class(1:3:end,1);
save('weizmann_run_walk_medium','I','class');