trials_num = 10000;


[X,Y] = meshgrid(0.1:0.1:0.9, 0.5:0.05:0.9);
Z_1d = data(:,3);
Z = reshape(Z_1d, 9, 9);
Z = Z' ./ trials_num;

s = surf(X,Y,Z);
xlabel('\gamma');
ylabel('Σ Util');
zlabel('\Delta S');

colorbar