%======================================================================
%> @file calculateRMSE.m
%> @brief Function to calculate RMSE
%>
%> @author Marlies Nitschke
%> @date December, 2018
%======================================================================

%======================================================================
%> @brief Function to calculate RMSE
%>
%> @details
%> This function will calculate the Root-Mean-Squared Error (RMSE). 
%> It indicates the absolute fit of two signals.
%> Lower values of RMSE indicate better fit.
%>
%>
%> \f$ \mathbf{x} = \left[x_1,x_2,\ldots,x_N \right]^T \f$
%>
%>
%> \f$ \mathbf{y} = \left[y_1,y_2,\ldots,y_N \right]^T  \f$  
%>
%>
%> RMSE = \f$ \sqrt {\frac{1}{N} \displaystyle\sum_{i=1}^{N}\displaystyle (y_i - x_i )^2 } \f$
%>
%>
%> @code
%> RMSE = calculateRMSE(x, y)
%> @endcode
%>
%>
%>
%> @param  x         Double vector x
%> @param  y         Double vector y
%> @retval RMSE      Double: Root mean squared error
%======================================================================
function RMSE = calculateRMSE(x, y)

x(:,1) = x;
y(:,1) = y;

x(isnan(y)) = [];
y(isnan(y)) = [];
if (length(x) ~= length(y))
    msg = 'Invalid input: Signals do not have the same length';
    error(msg)
end

RMSE = sqrt(mean((x-y).^2));

end

