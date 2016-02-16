% Main Script for testing the performance of SDE, sSDE, sSDE-l, SVDM, sLLE.
% By Fei Xiong
%    ECE Dept.,
%    Northeastern University
%    2013-10-05
% matlabpool open 8
path_data='/home/fei/Dropbox/Fei_Mengran/Matlab Code/sSDE/dataset/';
name_data={ 'mug_silhouette','stapler_silhouette','mouse_silhouette', 'Ionosphere','PAL',...
    'sonar','heart', 'weizmann_lena_run_skip', 'weizmann_run_skip_medium', 'swissroll3D_type2','face_grid'};
dim_thresh=0.985;
for i=1:11
    [validation_fold,test_fold,SDESVMKMFStruct,error_rate_SDESVM1] ...
        = test_SDESVMKMF([path_data name_data{i}],dim_thresh, 'LLW','S');
    [validation_fold,test_fold,SDEStruct,error_rate_SDE1] = test_SDE([path_data name_data{i}],dim_thresh, 'LLW','S');
    [validation_fold,test_fold,SDESVMKMFStruct,error_rate_SDESVM1] ...
        = test_SDESVM([path_data name_data{i}],dim_thresh, 'LLW','S');
    [validation_fold,test_fold,SDEStruct,error_rate_SDE1] = test_SDE_test_only([path_data name_data{i}],dim_thresh, 'LLW','S');
    
    [validation_fold,test_fold,SDEStruct,error_rate_SDE1] = test_sLLEv1([path_data name_data{i}], 'LLW','S');
    [validation_fold,test_fold,SDMStruct,error_rate_SDM] = test_SDM([path_data name_data{i}]);
    [validation_fold,test_fold,SDMStruct,error_rate_SDM] = test_SDMv2([path_data name_data{i}]);
    [validation_fold,test_fold,SDEStruct,error_rate_SDE1] = test_sLLEv2([path_data name_data{i}], 'LLW','S');
end