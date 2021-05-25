trials_num = 10000;

x = data2(:,1);
y1 = data2(:,2) ./ trials_num;
y2 = data2(:,3) ./ trials_num;

plot(x, [y1,y2], 'o-')

xlabel('Σ Util');
ylabel('Schedulability');

legend(["without MID mode","with MID mode"])
