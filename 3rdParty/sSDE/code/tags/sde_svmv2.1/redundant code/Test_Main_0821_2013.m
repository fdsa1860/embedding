warning off all
clear;
clc
matlabpool close;
matlabpool open ;
clear
name_data='swissroll3D_type2';
dim_thresh=0.985;

[validation_fold,test_fold,SDESVMStruct,error_rate_SDESVM1] = test_SDESVM(name_data,dim_thresh, 'MF','S');
[validation_fold,test_fold,SDESVMStruct,error_rate_SDESVM1] = test_SDESVMKMF_nl_effect(name_data,dim_thresh, 'MF','S');
[validation_fold,test_fold,SDESVMStruct,error_rate_SDESVM1] = test_SDESVMKMF(name_data,dim_thresh, 'MF','S');
[validation_fold,test_fold,SDMStruct,error_rate_SDM] = test_SDM(name_data);
[validation_fold,test_fold,SDEStruct,error_rate_SDE1] = test_SDE(name_data,dim_thresh, 'MF','S');


%%

[validation_fold,test_fold,SDESVMStruct,error_rate_SDESVM1] = test_SDESVMKMF_nl_effect_fullconstraint(name_data,dim_thresh, 'MF','S');


%%
name_data='swissroll3D_type2';
display('============SDM==============')
[validation_fold,test_fold,SDMStruct,error_rate_SDM] = test_SDM(name_data);

name_data='PAL';
display('============SDM==============')
[validation_fold,test_fold,SDMStruct,error_rate_SDM] = test_SDM(name_data);

name_data='Ionosphere';
display('============SDM==============')
[validation_fold,test_fold,SDMStruct,error_rate_SDM] = test_SDM(name_data);

name_data='heart';
display('============SDM==============')
[validation_fold,test_fold,SDMStruct,error_rate_SDM] = test_SDM(name_data);

name_data='sonar';
display('============SDM==============')
[validation_fold,test_fold,SDMStruct,error_rate_SDM] = test_SDM(name_data);