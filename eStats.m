function [RMSE, Residuals, Av_residual]= eStats(real, predict) 
% eStats(real,predict)
% Calculates the RMSE value and R squared value 

N = length(real); 
diff = zeros(N,1);
for i = 1:N
    diff(i) = abs(real(i) - predict(i)); 
end 
Residuals = diff; 
Av_residual = sum(diff)/N; 
%% RMSE
squared_diff = diff.^2; 
sum_squared_diff = sum(squared_diff);
MSE = sum_squared_diff / N; 

RMSE = sqrt(MSE); 

