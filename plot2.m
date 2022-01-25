close all
clear

data = load("output/data2_new.txt");

trials_num = 1000;

low_task_num_v = data(:,9);


x = data(:,1);


y = data(:,4) ./ low_task_num_v;
z = data(:,5) ./ low_task_num_v;
k = data(:,6) ./ low_task_num_v;
m = data(:,7) ./ low_task_num_v;
n = data(:,8) ./ low_task_num_v;

% plot survivability
plot(x, y, 'kx:')

hold on

plot(x, z, 'ko--')

hold on

plot(x, k, 'k^:')

hold on

plot(x, m, 'k+--')

hold on

plot(x, n, 'kx-')


xlabel('Σ Util');
ylabel('Survivability');
legend(["\gamma = 0.1", "\gamma = 0.3", "\gamma = 0.5", "\gamma = 0.7", "optimal"])

grid on