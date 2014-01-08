%% AR generating example, TODO: write function wrapper to generate all the required time series
% model1 = arima('Constant',0,'AR',0.3,'ARLags',1,'Variance',0.01);
% [Y1,E]=simulate(model1,1000);
% figure(1);
% plot(Y1);
% 
% model2 = arima('Constant',0,'AR',0.6,'ARLags',1,'Variance',0.01);
% [Y2,E]=simulate(model2,1000);
% figure(2);
% plot(Y2);

%% Simulate Process
model1 = arima('Constant',0,'AR',0.3,'ARLags',1,'Variance',0.01);
[Y1,E]=simulate(model1,1000);
[Y2,E]=simulate(model1,1000);
model2 = arima('Constant',0,'AR',0.6,'ARLags',1,'Variance',0.01);
[Y3,E]=simulate(model2,1000);
[Y4,E]=simulate(model2,1000);
D=[Y1,Y2,Y3,Y4];
clusterObj=cluster(D);
