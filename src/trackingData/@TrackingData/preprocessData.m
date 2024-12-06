%======================================================================
%> @file @TrackingData/preprocessData.m
%> @brief TrackingData function to preprocess tracking data
%> @details
%> Details: TrackingData::preprocessData()
%>
%> @author Marlies Nitschke
%> @date May, 2018
%======================================================================

% ======================================================================
%> @brief Function to preprocess tracking data
%>
%> @details
%> This function calls:
%> - TrackingData::resampleData()
%> - TrackingData::correctGyroSum()
%> - TrackingData::correctVariance()
%>
%> @param   obj       TrackingData class object which should be processed
%> @param   nNodes    Double: Number of nodes
% ======================================================================
function preprocessData(obj,nNodes)

% resample data
obj.resampleData(nNodes);

% check sum of gyro signal
obj.correctGyroSum();

% check if the variance is large enough
obj.correctVariance();
end