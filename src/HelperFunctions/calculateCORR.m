%======================================================================
%> @file calculateCORR.m
%> @brief Function to caluclate the Pearson correlation coefficient
%>
%> @author Marlies Nitschke
%> @date January, 2021
%======================================================================

%======================================================================
%> @brief Function to caluclate the Pearson correlation coefficient
%>
%> @details
%> The Pearson correlation coefficient does not detect systematic biases.
%> It quantifies the strength of linear combinaton and ranges from -1 indicating 
%> perfect negative correlaton to 1 indicating perfect positive correlation.  
%>
%> \f$ r= \displaystyle\frac{\displaystyle\sum_{i=1}^{N}(x_i - \bar{x})(y_i - \bar{y}) }{\displaystyle\sqrt{\sum_{i=1}^{N}(x_i - \bar{x})^2} \sqrt{\sum_{i=1}^{N}(y_i - \bar{y})^2}} \f$
%>
%> @code
%> out_r = calculateCORR(x, y)
%> @endcode
%>
%>
%>
%> @param  x       Double Vector: x signal 
%> @param  y       Double Vector: y signal
%> @retval out_r   Double: Pearson correlation coefficient which is between -1 and 1 
%======================================================================

function out_r = calculateCORR(x, y)

if (length(x) ~= length(y))
    msg = 'The two input variables do not have the same length.';
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

out_r = r;

end

