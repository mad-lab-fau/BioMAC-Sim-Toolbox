%======================================================================
%> @file calculateCMC.m
%> @brief Function to caluclate the Coefficient of Multiple Correlation 
%>
%> @author Marlies Nitschke
%> @date December, 2018
%======================================================================

%======================================================================
%>  @brief Function to caluclate the Coefficient of Multiple Correlation 
%>
%> @details
%> This function will caluclate the Coefficient of Multiple Correlation
%> or CMC for statistical analysis, in order to compare 
%> two signals \f$ Y_{1f} \f$ and \f$ Y_{2f} \f$. It is used to assess
%> the similarity of waveforms.
%> A higher value indicates higher similarity.
%>
%>
%> \f$
%> CMC = \sqrt{1 - \frac{\displaystyle\sum_{m=1}^{M}\sum_{f=1}^{F} (Y_{mf}-\bar{Y}_f)^2 / (F (M-1) )}{\displaystyle\sum_{m=1}^{M}\sum_{f=1}^{F} (Y_{mf}-\bar{Y})^2 / (MF-1) )}  }
%> \f$
%> 
%>
%> \f$ \bar{Y}_f = \frac {1}{M} \displaystyle\sum_{m=1}^{M}Y_{mf} \f$
%> 
%>
%> \f$ \bar{Y} = \frac {1}{MF} \displaystyle\sum_{m=1}^{M}\sum_{f=1}^{F}Y_{mf} \f$
%>
%>
%> In this function:
%> \f$ M=2 \Rightarrow \f$ 
%>
%>
%> \f$ Y_{1f} = \left[y_{11},y_{12},\ldots,y_{1F} \right]^T \f$ 
%>
%>
%> \f$ Y_{2f} = \left[y_{21},y_{22},\ldots,y_{2F} \right]^T \f$
%>
%>
%> @code
%> CMC_value = calculateCMC(y1f, y2f)
%> @endcode
%>
%>
%> Warning:
%> - If the IMC and OMC signals are not of the same size, an
%>   error will be thrown.
%>
%> @param  y1f             Double vector: Signal 1
%> @param  y2f             Double vector: Signal 2
%> @retval CMC_value       Double: Coefficient of mulitiple correlation which is between 0 and 1
%======================================================================




function CMC_value = calculateCMC(y1f, y2f)
%calculates the cmc value of two signals
y1f(:,1) = y1f;
y2f(:,1) = y2f;

if (length(y1f) ~= length(y2f))
    msg = 'Invalid input: Signals do not have the same length';
    error(msg)
end

F = length(y1f); %number of frames
numera = 0;
denom = 0;

grand_mean = 0;

for o = 1:F
    grand_mean = y1f(o,1) + y2f(o,1);
end
grand_mean = grand_mean * (1/(2*F));

for k = 1:F
    mean = 0.5 * (y1f(k,1) + y2f(k,1));
    numera = numera + (y1f(k,1)- mean)*(y1f(k,1)- mean) + (y2f(k,1)-mean)*(y2f(k,1)-mean);
    denom = denom + (y1f(k,1)-grand_mean)*(y1f(k,1)-grand_mean) + (y2f(k,1)-grand_mean)*(y2f(k,1)-grand_mean);
end
numera = numera/F;
denom = denom/ (2*F-1);

value = 1- (numera/denom);
CMC_value = real(sqrt(value));


end

