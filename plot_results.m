csvFile1 = "distances_z1998_p0_w90_g2_c999.csv"; % Your CSV file
%csvFile1 = "distances_z2151_p797_w90_g2_c95.csv"; % Your CSV file
data2 = csvread(csvFile1);
figure()
subplot(1,2,1);
plot(data2(5:end,1),data2(5:end,2),LineWidth=2); hold on;


ylabel("count");
xlabel("interparticle distance");
title("interdot spacing distribution")
xlim([1.5,4]);
grid on;

subplot(1,2,2);
bar(data2(1:30,3),data2(1:30,4)); hold on;

title("dot z-distribution")
ylabel("dots");
xlabel("z position")
grid on