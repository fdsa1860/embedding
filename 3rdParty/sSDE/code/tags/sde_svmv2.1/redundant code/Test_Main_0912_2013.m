warning off all
clear;
clc
matlabpool close;
matlabpool open ;
clear
name_data='swissroll3D_type2';
dim_thresh=0.985;
[validation_fold,test_fold,SDESVMStruct,error_rate_SDESVM1] = test_SDESVM(name_data,dim_thresh, 'MF','S');
[validation_fold,test_fold,SDESVMStruct,error_rate_SDESVM1] = test_SDESVMKMF_nl_effect(name_data,fname, dim_thresh, 'LLW','S');
[validation_fold,test_fold,SDESVMStruct,error_rate_SDESVM1] = test_SDESVM_nl_only(name_data,fname, dim_thresh, 'MF','S');
[validation_fold,test_fold,SDEStruct,error_rate_SDE1] = test_SDE(name_data,dim_thresh, 'MF','S');
[validation_fold,test_fold,SDMStruct,error_rate_SDM] = test_SDM(name_data);
[validation_fold,test_fold,SDESVMKMFStruct,error_rate_SDESVM1] ...
    = test_SDESVMKMF(name_data,dim_thresh, 'MF','S')
%%
name_data='swissroll3D_type2';
dim_thresh=0.985;
display('============SDM==============')
[validation_fold,test_fold,SDEStruct,error_rate_SDE1] = test_sLLEv2(name_data, 'LLW','S');
[validation_fold,test_fold,SDEStruct,error_rate_SDE1] = test_SDE_test_only(name_data,dim_thresh, 'LLW','S');
[validation_fold,test_fold,SDMStruct,error_rate_SDM] = test_SDMv2(name_data);

%%
name_data={ 'Ionosphere', 'swissroll3D_type2','PAL','heart','sonar'}
for i=1:length(name_data)
    [validation_fold,test_fold,SDEStruct,error_rate_SDE1] = test_sLLE(name_data{i}, 'LLW','S');
end
%%
name_data={ 'Ionosphere', 'swissroll3D_type2','PAL','heart','sonar'}
for i=1:length(name_data)
    [validation_fold,test_fold,SDEStruct,error_rate_SDE1] = test_sLLEv1(name_data{i}, 'LLW','S');
end
%%
name_data={ 'mug_silhouette', 'Ionosphere','PAL','heart','sonar', 'swissroll3D_type2'};
for i=1:length(name_data)
    [validation_fold,test_fold,SDEStruct,error_rate_SDE1] = test_sLLEv2(name_data{i}, 'LLW','S');
end
%%
name_data={ 'mug_silhouette','face_grid', 'Ionosphere','PAL',...
    'sonar','heart', 'swissroll3D_type2','stapler_silhouette','mouse_silhouette'};
for i=1:length(name_data)
    [validation_fold,test_fold,SDMStruct,error_rate_SDM] = test_SDMv2(name_data{i});
end
%%
name_data={'swissroll3D_type2', 'Ionosphere', 'PAL','sonar','heart'};
dim_thresh=0.985;
for i=3:length(name_data)
    [validation_fold,test_fold,SDESVMKMFStruct,error_rate_SDESVM1] ...
    = test_SDESVMKMF(name_data{i},dim_thresh, 'LLW','S');
end

%%
name_data={ 'mug_silhouette','stapler_silhouette','mouse_silhouette', 'Ionosphere','PAL',...
    'sonar','heart', 'weizmann_lena_run_skip', 'weizmann_run_skip_medium', 'swissroll3D_type2','face_grid'};
dim_thresh=0.985;
for i=1:length(name_data)-1
    [validation_fold,test_fold,SDMStruct,error_rate_SDM] = test_SDE_test_only(name_data{i},dim_thresh, 'LLW','S');
end
%%
name_data={'Ionosphere', 'PAL','heart','sonar','swissroll3D_type2'};
dim_thresh=0.985;
for i=1:length(name_data)
    [validation_fold,test_fold,SDESVMKMFStruct,error_rate_SDESVM1] ...
    = test_SDESVM(name_data{i},dim_thresh, 'LLW','S');
end
%%
name_data={'swissroll3D_type2', 'Ionosphere', 'PAL','heart','sonar'};
dim_thresh=0.985;
for i=3:length(name_data)
    [validation_fold,test_fold,SDEStruct,error_rate_SDE1] = test_SDE(name_data{i},dim_thresh, 'LLW','S');
end
%%
name_data={'swissroll3D_type2', 'Ionosphere', 'PAL','heart','sonar'};
dim_thresh=0.985;
for i=3:length(name_data)
    [validation_fold,test_fold,SDMStruct,error_rate_SDM] = test_SDM(name_data{i});
end
%%
name_data={'weizmann_run_skip_medium','weizmann_lena_run_skip', 'mug_silhouette','PAL','face_grid'};
dim_thresh=0.985;
for i=1:1
    [validation_fold,test_fold,SDESVMKMFStruct,error_rate_SDESVM1] ...
        = test_SDESVMKMF(name_data{i},dim_thresh, 'LLW','S');
    [validation_fold,test_fold,SDEStruct,error_rate_SDE1] = test_SDE(name_data{i},dim_thresh, 'LLW','S');
    [validation_fold,test_fold,SDESVMKMFStruct,error_rate_SDESVM1] ...
        = test_SDESVM(name_data{i},dim_thresh, 'LLW','S');
    [validation_fold,test_fold,SDEStruct,error_rate_SDE1] = test_SDE_test_only(name_data{i},dim_thresh, 'LLW','S');
    
    [validation_fold,test_fold,SDEStruct,error_rate_SDE1] = test_sLLEv1(name_data{i}, 'LLW','S');
    [validation_fold,test_fold,SDMStruct,error_rate_SDM] = test_SDM(name_data{i});
    [validation_fold,test_fold,SDMStruct,error_rate_SDM] = test_SDMv2(name_data{i});
    [validation_fold,test_fold,SDEStruct,error_rate_SDE1] = test_sLLEv2(name_data{i}, 'LLW','S');
end
%%
