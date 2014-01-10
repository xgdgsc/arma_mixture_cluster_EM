%% AR generating example, 
% TODO: write function wrapper to generate all the required time series
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
[Y3,~]=simulate(model1,800);
[Y4,~]=simulate(model1,999);
[Y5,~]=simulate(model1,1200);
model2 = arima('Constant',0,'AR',0.6,'ARLags',1,'Variance',0.01);
[Y6,~]=simulate(model2,1000);
[Y7,~]=simulate(model2,999);
[Y8,~]=simulate(model2,700);
[Y9,~]=simulate(model2,900);
[Y10,~]=simulate(model2,800);

%full time series data cell D
D=cell(10,1);
D{1,1}=Y1;
D{2,1}=Y2;
D{3,1}=Y3;
D{4,1}=Y4;
D{5,1}=Y5;
D{6,1}=Y6;
D{7,1}=Y7;
D{8,1}=Y8;
D{9,1}=Y9;
D{10,1}=Y10;

clusterObj=cluster(D);
clusterObj.initialize(1,0,2);
clusterObj.EM();