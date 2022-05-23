trials_num = 10000;

data = load("output/data1_new.txt");

[X,Y] = meshgrid(0.1:0.1:0.9, 0.5:0.05:0.9);
Z_1d = data(:,4) ./ data(:,5) * 100;
Z = reshape(Z_1d, 9, 9);
Z = Z' ./ trials_num;

s = surf(X,Y,Z);
xlabel('\bf \gamma');
ylabel('Σ Util');
zlabel('\Delta S');

colorbar

view(90,90)