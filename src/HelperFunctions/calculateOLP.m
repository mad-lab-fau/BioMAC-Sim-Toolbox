%======================================================================
%> @file calculateOLP.m
%> @brief Function to caluclate the OLP Regression
%>
%> @author Marlies Nitschke
%> @date December, 2018
%======================================================================

%======================================================================
%> @brief Function to caluclate the OLP Regression
%>
%> @details
%> The OLP regression analysis exposes systematic differences of two 
%> methods instead of similarities. The coefficient considers that both methods can 
%> be affected by random errors.
%>
%>
%>
%> \f$ \mathbf{x} = \left[x_1,x_2,\ldots,x_N \right]^T \f$
%>
%>
%> \f$ \mathbf{y} = \left[y_1,y_2,\ldots,y_N \right]^T  \f$ 
%>
%>
%> RMSE = \f$ \sqrt {\frac{1}{N} \displaystyle\sum_{i=1}^{N}\displaystyle (y_i - (\beta x_i + \alpha))^2 } \f$
%>
%>
%> \f$ \alpha = \bar{y} - \beta\bar{x} \f$
%>
%>
%> \f$ \beta = \displaystyle\sqrt {\frac {\beta_y}{\beta_x}} \f$
%> 
%>
%> \f$ \beta_y = \displaystyle\frac {\displaystyle\sum_{i=1}^{N}(x_i - \bar{x})(y_i - \bar{y}) }{\displaystyle\sum_{i=1}^{N}(x_i - \bar{x})^2} \f$
%>
%>
%> \f$ \beta_x = \displaystyle\frac {\displaystyle\sum_{i=1}^{N}(x_i - \bar{x})(y_i - \bar{y}) }{\displaystyle\sum_{i=1}^{N}(y_i - \bar{y})^2} \f$
%> 
%>
%> Relative RMSE = \f$ \displaystyle\frac{RMSE}{\displaystyle\frac{1}{2} \Bigg\{ [\max_{0<t<T}(x_i(t)) - \min_{0<t<T}(x_i(t)) ] + [\max_{0<t<T}(y_i(t)) - \min_{0<t<T}(y_i(t)) ] \Bigg\} } * 100\% \f$
%>
%>
%> @code
%> [ out_RMSE, out_beta, out_alpha, out_r,out_DRMSE] = calculateOLP(x, y)
%> @endcode
%>
%>
%>
%> @param  x               Double Vector: x signal 
%> @param  y               Double Vector: y signal
%> @retval out_RMSE        Double: RMSE : between 0 and 1
%> @retval out_beta        Double: slope: any real value
%> @retval out_alpha       Double: y-intercept: any real value
%> @retval out_r           Double: PPMCC (Pearson's r): between -1 and 1 
%> @retval out_DRMSE       Double: relative RMSE: between -100% and 100%
%======================================================================

function [ out_RMSE, out_beta, out_alpha, out_r,out_DRMSE] = calculateOLP(x, y)

%calculates the OLP regression analysis
% out_beta is slope
% out_alpha is y-intercept
% r_out is Pears product-moment correlation coefficient
% out_RMSE is root mean squared error


if (length(x) ~= length(y));
    msg = 'Vicon and Model trajectories have not the same length';
    error(msg);
end

x_i = x;
y_i = y;
n = length(x_i);

mean_x = mean(x_i);
mean_y = mean(y_i);

sum_x = 0;  % sum of (xi - mean_x)^
sum_y = 0;  %sum of (yi - mean_y)^
product_x_y = 0;    %sum of (xi - mean_x)*(yi - mean_y)

for i = 1 : n
    sum_x = sum_x + ((x_i(i)-mean_x)^2);
    sum_y = sum_y + ((y_i(i)-mean_y)^2);
    product_x_y = product_x_y + ((x_i(i)-mean_x)*(y_i(i)-mean_y));
end

r = product_x_y/( sqrt(sum_x)*sqrt(sum_y));

beta_y = product_x_y/sum_x;
beta_x = product_x_y/sum_y;

beta = sqrt(beta_y/beta_x);
alpha = mean_y - beta*mean_x;
% 
% beta = 1;
% alpha = 0;
rmse = 0;

for j = 1 : n
    rmse = rmse + (y_i(j) - (beta*x_i(j)+alpha))^2;
end


rmse = sqrt(rmse/n);
out_RMSE = rmse;
out_DRMSE = rmse/(1/2*(max(x)-min(x)+max(y)-min(y)))*100;
out_beta = beta;
out_alpha = alpha;
out_r = r;


end

