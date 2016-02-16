% Script to read and preprocess the silhouette image from the Tabletop_3D
% dataset( published by U Mich)
% By Fei Xiong
%    ECE Dept.,
%    Northeastern University
%    2013-10-05

%% load mask
% path = 'D:\DownloadedDataset\table_top_object_detection_dataset_UMIG\masks\';
path = '/home/fei/Downloads/Dataset/TableTopObjectDataset/Table-Top-Pose/mouse/masks/';
% path = '/home/fei/Downloads/Dataset/TableTopObjectDataset/Table-Top-Pose/mug/masks/';
% path = '/home/fei/Downloads/Dataset/TableTopObjectDataset/Table-Top-Pose/stapler/masks/';
dirc =dir(path);

strct = struct('Angle',[],'Height',[],'Scale',[]);
view_point = repmat(strct, length(dirc)-2,1);
k =1;
for i =3: length(dirc)
    if ~dirc(i).isdir
        pos = strfind(dirc(i).name,'_'); %'%s_%d_A%d_H%d_S%d-%s');
        view_point(k).Angle = str2num(dirc(i).name(pos(2)+2));
        view_point(k).Height = str2num(dirc(i).name(pos(3)+2));
        view_point(k).Scale = str2num(dirc(i).name(pos(4)+2));
        Itemp =imread([path dirc(i).name ]);
        Itemp = im2bw(Itemp);
        stat =regionprops(Itemp,'BoundingBox');
        rect_temp = stat(1).BoundingBox;
        
        Itemp = Itemp(rect_temp(2):rect_temp(4)+rect_temp(2)-1, rect_temp(1):rect_temp(1)+rect_temp(3)-1 );
        Itemp = imresize(Itemp, [128, 96]); % scale the image size.
        norm_data(:,k)= Itemp(:);
        class(k,1)= view_point(k).Height;
        k =k+1;
    end
end
norm_data = double(norm_data(:,1:k-1));
view_point = view_point(1:k-1);
norm_data = norm_data -repmat(mean(norm_data, 2),1, size(norm_data,2));
%%
cmap = colormap('jet');
Y = lle(double(norm_data),8, 3);
Y=Y';
subplot(1,2,1)
for i =1: size(Y,2)
c_idx =round(view_point(i).Height*8);
plot3(Y(1,i), Y(2,i), Y(3,i),'o', 'MarkerSize', 5, 'LineWidth',2, 'Color', cmap(c_idx,:)); hold on;
text(Y(1,i), Y(2,i), Y(3,i),num2str(view_point(i).Height));
end
subplot(1,2,2)
for i =1: size(Y,2)
c_idx =round(view_point(i).Angle*8);
plot3(Y(1,i), Y(2,i), Y(3,i),'o', 'MarkerSize', 5, 'LineWidth',2, 'Color', cmap(c_idx,:)); hold on;
text(Y(1,i), Y(2,i), Y(3,i),num2str(view_point(i).Angle));
end