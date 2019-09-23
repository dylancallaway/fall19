close all
clear
clc

x = [1, 2, 3, 4];
for j = 1:length(x)
    j
    if j == 1
        x(1) = [];
    end
end