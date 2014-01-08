%% AR generating example, %TODO: write function wrapper to generate all the required time series
% model1 = arima('Constant',0,'AR',0.3,'ARLags',1,'Variance',0.01);
% [Y1,E]=simulate(model1,1000);
% figure(1);
% plot(Y1);
% 
% model2 = arima('Constant',0,'AR',0.6,'ARLags',1,'Variance',0.01);
% [Y2,E]=simulate(model2,1000);
% figure(2);
% plot(Y2);
clear;
%% Simulate Process
%generate test time series ( 4 as an initial sample)
model1 = arima('Constant',0,'AR',0.3,'ARLags',1,'Variance',0.01);
[Y1,~]=simulate(model1,1000);
[Y2,~]=simulate(model1,999);
model2 = arima('Constant',0,'AR',0.6,'ARLags',1,'Variance',0.01);
[Y3,~]=simulate(model2,1000);
[Y4,~]=simulate(model2,999);
%full time series data cell D
D=cell(4,1);
D{1,1}=Y1;
D{2,1}=Y2;
D{3,1}=Y3;
D{4,1}=Y4;
clusterObj=cluster(D);
clusterObj.initialize(1,0,2);