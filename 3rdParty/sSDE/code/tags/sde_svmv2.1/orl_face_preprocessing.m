% Script to generate the shifted face images
% By Fei Xiong
%    ECE Dept.,
%    Northeastern University
%    2013-10-05
N =960;
% xx =random('uniform', 0,1, [N 1]);
% yy =random('uniform', 0,1, [N 1]);
[xx yy]= meshgrid(1:27,1:22);
xx= xx*4;
yy=yy*4;
N =length(xx(:));
for i =1:N
    % xxx = ceil(xx(i)*(200-92));
    % yyy = ceil(yy(i)*(200-112));
    xxx= xx(i);
    yyy= yy(i);
    Itemp = uint8(randn(200,200)*0+127);
    Itemp(yyy:yyy+size(x1,1)-1,xxx:xxx+size(x1,2)-1) = x1;
    I(:,i)=Itemp(:);
end
imwrite(imresize(reshape(I(:,1),200,200),[75,75]),'left_top.png')
imwrite(imresize(reshape(I(:,22),200,200),[75,75]),'left_bottom.png')
imwrite(imresize(reshape(I(:,594),200,200),[75,75]),'right_bottom.png')
imwrite(imresize(reshape(I(:,573),200,200),[75,75]),'right_top.png')
imwrite(imresize(reshape(I(:,297),200,200),[75,75]),'center.png')
norm_data = double(I);
for i =1:N
    if (xx(i)/112 -0.5)^2 + (yy(i)/88 -0.5)^2<0.15
        class (i)=1;
    else
        class (i)=2;
    end
end
figure
y = class*2-3;
idx = y ==1;
plot(xx(idx), yy(idx), 'ro', 'MarkerSize', 5, 'LineWidth',2)
hold on;
idx = y ==-1;
plot(xx(idx), yy(idx), 'bo', 'MarkerSize', 5, 'LineWidth',2)
grid on;

figure
[LD] = slle(I, class, 4, 6, 1);
Y = LD;
Y=Y';
idx = y ==1;
plot3(Y(1,idx), Y(2,idx), Y(3,idx),'ro', 'MarkerSize', 5, 'LineWidth',2)
hold on;
idx = y ==-1;
plot3(Y(1,idx), Y(2,idx), Y(3,idx),'bo', 'MarkerSize', 5, 'LineWidth',2)
grid on;
%%
k=0;
for i =1:27
    for j=1:22
        k =k+1;
        % text(Y(1, k), Y(2, k), Y(3, k), num2str(k))
        if j<22
            plot3([Y(1, k) Y(1, k+1)], [Y(2, k) Y(2, k+1)], [Y(3, k) Y(3, k+1)], 'k-'), hold on;
        end
        if i< 27
            plot3([Y(1, k) Y(1, k+22)], [Y(2, k) Y(2, k+22)], [Y(3, k) Y(3, k+22)], 'k-'), hold on;
        end
    end
end
%% ORL load
path = 'D:\DownloadedDataset\orl_faces\'
dirc =dir(path);
for i =3: length(dirc)
    if dirc(i).isdir
        dirc_temp =dir([path dirc(i).name]);
        k =1;
        for j =1:length(dirc_temp)
            
            if ~dirc_temp(j).isdir
                Itemp =imread([path dirc(i).name '\' dirc_temp(j).name]);
                I{i}(k,:)= Itemp(:);
                k=k+1;
            end
        end
    end
end
idx =ones(length(I),1);
for i =1:length(I)
    if isempty(I{i})
        idx(i) =0;
    end
end
I = I(find(idx));

%% UMIST face preprocessing
count = 0;
for i =1:length(facedat)
    count =count + size(facedat{i},3);
end

norm_data= zeros(size(facedat{i},1)*size(facedat{i},2), count);
k=1;
for i =1:length(facedat)
    for j =1: size(facedat{i},3)
        temp = squeeze(facedat{i}(:,:,j));
        temp = histeq(temp);
        norm_data(:,k) = temp(:);
%         class(k) = -1;
%         if j> size(facedat{i},3)/2
%             class(k) = 1;
%         end
        view_point(k) = j / size(facedat{i},3);
        k = k+1;
    end
end

norm_data = norm_data - repmat(mean(norm_data,2), 1, count);