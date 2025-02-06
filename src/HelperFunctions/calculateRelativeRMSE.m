%======================================================================
%> @file calculateRelativeRMSE.m
%> @brief Function to calculate Relative RMSE
%>
%> @author Marlies Nitschke
%> @date December, 2018
%======================================================================

%======================================================================
%> @brief Function to calculate Relative RMSE
%>
%> @details
%> This function will calculate the relative Root-Mean-Squared Error (RRMSE) or
%> normalized RMSE. It indicates the relative fit of the model, that is, how close the observed 
%> data points are to the model's predicted values. It is normalized with respect to the average peak-to-
%> peak amplitude. If the input is not a vector but a scalar, the average is used for normalization.
%> As with the RMSE, lower relative RMSE values indicate better fit. 
%>
%>
%> \f$ \mathbf{x} = \left[x_1,x_2,\ldots,x_N \right]^T \f$ 
%>
%>
%> \f$ \mathbf{y} = \left[y_1,y_2,\ldots,y_N \right]^T \f$  
%>
%>
%> RMSE = \f$ \sqrt {\frac{1}{N} \displaystyle\sum_{i=1}^{N}\displaystyle (y_i - x_i )^2 } \f$
%>
%>
%> Relative RMSE = \f$ \displaystyle\frac{RMSE}{\displaystyle\frac{1}{2} \sum_{i=1}^{2}[\max_{0<t<T}(u_i(t)) - \min_{0<t<T}(u_i(t)]  } * 100\% \f$, if x and y are vectors
%>
%> Relative RMSE = \f$ \displaystyle\frac{RMSE}{\displaystyle\frac{1}{2} \sum_{i=1}^{2}u_i(t)  } * 100\% \f$, if x and y are scalars
%>
%>
%> @code
%> [RMSE] = calculateRelativeRMSE(x, y)
%> @endcode
%>
%>
%>
%> @param  x        Double Vector x
%> @param  y        Double Vector y
%> @retval RMSE     Double: relative RMSE, between -100% and 100%)
%======================================================================

function [ RMSE] = calculateRelativeRMSE(x, y)
%calculates relative RMSE value of two inputs
x(:,1) = x;
y(:,1) = y;

x(isnan(y)) = [];
y(isnan(y)) = [];
if (length(x) ~= length(y))
    msg = 'Invalid input: inputs do not have the same length';
    error(msg)
end

F = length(x); %number of frames

value = 0;
for k = 1:F
   value = value + (x(k) -y(k))^2;
end
max(x);

%RMSE = real(sqrt(value/F));
if numel(x) > 1
    RMSE = sqrt(mean((x-y).^2))/(1/2*(max(x)-min(x)+max(y)-min(y))) * 100;
else
    RMSE = sqrt(mean((x-y).^2))/(1/2*(x+y)) * 100;
end

end

