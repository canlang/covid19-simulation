clearvars;close all;
load('rep100_result.mat')
su_list = 0:.1:1;


boxplot(extinct_date')
xticklabels(arrayfun(@(x) num2str(x, '%.1f'),su_list,'UniformOutput',false))