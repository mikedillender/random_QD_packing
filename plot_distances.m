clear; figure()
csvFile = 'cmake-build-debug\distances_N3276_w32.csv'; % Your CSV file
data1 = csvread(csvFile);
plot(data1(5:end,1),data1(5:end,2));