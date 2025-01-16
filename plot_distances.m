clear; figure()
%{
csvFile = 'cmake-build-debug\distances2_N141120_w210_70_2g.csv'; % Your CSV file
data1 = csvread(csvFile);
plot(data1(5:end,1),data1(5:end,2));
hold on;
csvFile1 = 'cmake-build-debug\distances2_N141120_w210_95_2g.csv'; % Your CSV file
data2 = csvread(csvFile1);
plot(data2(5:end,1),data2(5:end,2));
legend(".7",".9")
hold on;
%}
csvFile1 = 'cmake-build-debug\distances2_N141120_w210_95_1g.csv'; % Your CSV file
data2 = csvread(csvFile1);
plot(data2(5:end,1),data2(5:end,2),LineWidth=2);
hold on;
csvFile1 = 'cmake-build-debug\distances_p_N141120_w210_95_1g.csv'; % Your CSV file
data2 = csvread(csvFile1);
plot(data2(5:end,1),data2(5:end,2),LineStyle="--",LineWidth=2);

csvFile1 = 'cmake-build-debug\distances_pf_N141120_w210_95_1g.csv'; % Your CSV file
data2 = csvread(csvFile1);
plot(data2(5:end,1),data2(5:end,2),LineStyle="--",LineWidth=2);
%legend(".7,2",".95,2",'.95,1',".95,1,P",Location="best")
legend("Boxed","Periodic","Periodic,fixed",Location="best");

ylabel("count");
xlabel("interparticle distance");
title("interdot spacing distribution")
grid on;
xlim([1.5,4]);