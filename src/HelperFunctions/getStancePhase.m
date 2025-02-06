%======================================================================
%> @file getStancePhase.m
%> @brief Function to get the stance phase from vertical GRF
%> @details
%> Details: getStancePhase()
%>
%> @author Marlies Nitschke
%> @date September, 2018
%======================================================================

%======================================================================
%> @brief Function to get the stance phase from vertical GRF
%>
%> @details
%> Uses C. Maiwald, T. Sterzing, T. a. Mayer, and T. L. Milani,
%> Detecting foot-to-ground contact from kinematic data in running,
%> Footwear Sci., vol. 1, no. 2, pp. 111-118, Jun. 2009.
%> => "True touch down and take off events were initially determined from force
%>    plate data by applying a threshold of 10N for touch down and 5 N for take
%>    off (O'Connor et al. 2007) to the vertical GRF data."
%>
%> Warning: The algorithm is not perfect up to now. Sometime there are
%> errors in the result.
%>
%> @param   verticalGRF  Double array: Vertical ground reaction force in N
%> @retval  standing     Double array: 0: no standing; 1: standing (length of verticalGRF)
%> @retval  idxHSs       Double array: Indices of heel strikes (length is number of heel strikes)
%> @retval  idxTOs       Double array: Indices of toe offs (length is number of toe offs)
%======================================================================
function [standing, idxHSs, idxTOs] = getStancePhase(verticalGRF)

%% get HS at 10 N
idx10N = (verticalGRF > 10);
idxHSs = find(diff(idx10N) > 0.5)+1; % rising edge+1: first sample above 10 N

%% get TO at 5 N
idx05N = (verticalGRF > 5);
idxTOs = find(diff(idx05N) < -0.5)-1; % falling edge-1: last sample above 5 N

%% get standing phases
standing = zeros(size(verticalGRF));

% sort the indices in ascendig order
indices = [idxHSs', idxTOs'];
flags = [zeros(size(idxHSs')), ones(size(idxTOs'))]; % 0: HS; 1: TO

[indices_sort, idx] = sort(indices);
flags_sort = flags(idx);

% remove if there are consecutive values 
if isempty(flags_sort)
   return; 
end
indices_sort = indices_sort([1,diff(flags_sort)]~=0);
flags_sort= flags_sort([1,diff(flags_sort)]~=0);

% get standing
i = 1;
while  i <= length(indices_sort)
    
    if i == 1 && flags_sort(1) == 1 % TO as first event
        standing(1:indices_sort(1)) = 1;
        i = i +1;
    elseif i == length(indices_sort) && flags_sort(end) == 0 % HS as last event
        standing(indices_sort(i):end) = 1;
        i = i +1;
    elseif flags_sort(i) == 0 && flags_sort(i+1) == 1 % HS to TO
        standing(indices_sort(i):indices_sort(i+1)) = 1;
        i = i +2; % skip the TO
    else
        warning('There was a wrong detection');
    end
    
end


end