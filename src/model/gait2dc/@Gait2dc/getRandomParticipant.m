%======================================================================
%> @file @Gait2dc/getRandomParticipant.m
%> @brief Gait2dc function to generate a participant with random muscle parameters
%> @details
%> Details: Gait2dc::getRandomParticipant()
%>
%> @author Marlies Nitschke, Eva Dorschky
%> @date September 23, 2018
%======================================================================

%======================================================================
%> @brief Function to generate a participant with random muscle parameters
%>
%> @details
%> This function generates a new model object where muscle parameters 
%> were independenlty normally distributed around the muscle parameters of 
%> obj using the standard deviation SDmusParam. In Dorschky et al.,
%> "Optimal control simulation predicts effects of midsole materials on energy cost of running"
%> (2019), we change fmax, lceopt, width, L0, kPEE, umax, vmax, tact, tdeact, gmax, 
%> arel, and the moment arms at each joint.
%>
%> We use the same muscle parameters for right and left muscles (symmetric model). Otherwise
%> it does not make sense to use the symmetry option in the simulation.
%>
%> To generate this new model, we use the excelfile, bodyheight, and bodymass 
%> of obj. Other parameters which were changed in obj will not be copied since 
%> we have to create a new model object. If we would not create a new model object, 
%> but make a deep copy by inheriting from matlab.mixin.Copyable, this would cause 
%> problems with the listener callbacks (They would not be called after changes of 
%> the muscle parameters.).
%> => If you find a better solution, please adapt this function.
%>
%> This function does not shuffle or reset the random number generator 
%> (see https://de.mathworks.com/help/matlab/ref/rng.html).
%>
%>
%> @param  obj         Gait2dc model object
%> @param  SDmusParam  Double: Standard deviation (SD) for distribution
%>                     around default value in % of the default value.
%>
%> @retval randModel   Gait2dc model object: New model generated using excelfile,
%>                     bodyheight and bodymass of obj. Afterwards muscle parameters 
%>                     were changed randomly following a normal distribution.
%======================================================================
function  randModel = getRandomParticipant(obj, SDmusParam)

% Set muscle parameters which are changed 
musParam = {'fmax', 'lceopt', 'width', 'L0', ...
    'dRhip', 'dRknee', 'dRankle', ... % list only the right leg and match it to the parameter of the lft side in musParamSym
    'kPEE', 'umax', 'vmax', 'tact', 'tdeact', 'gmax', 'arel'};
musParamSym = {'dRhip', 'dRknee', 'dRankle'; 'dLhip', 'dLknee', 'dLank'}';

% Convert SD from xx% to 0.xx
SDmusParam = SDmusParam / 100;

% Make a new object using the excelfile, bodyheight and bodymass from obj
objClass = class(obj); % Get the class such that we can use this function also for Gait2dc_Exo
if contains(objClass,'Exo')
    randModel = feval(objClass, obj.norm_value, obj.type, obj.increase, obj.excelfile, obj.bodyheight, obj.bodymass); 
else
    randModel = feval(objClass, obj.excelfile, obj.bodyheight, obj.bodymass); % Call the constructor of the specific class
end

% Get the muscle table and make first all changes. Otherwise the mex
% function will be updated after each single change.
randMus = randModel.muscles;

% Get symmetry of muscles
i_muscles = [];
for imuscles_r = 1:obj.nMus
    % find right muscles
    if ismember(obj.muscles.Properties.RowNames{imuscles_r}(end-1:end),'_r')
        % find corresponding left muscle
        imuscles_l = find(strcmp(obj.muscles.Properties.RowNames, [obj.muscles.Properties.RowNames{imuscles_r}(1:end-2),'_l']));
        % add right and left
        i_muscles = [i_muscles; imuscles_r imuscles_l];
    end
end

% Go trough all muscles and all muscle parameters to change the values
% -> Using a normal distribution with the given value as mean and the given SD.
% -> Due to the multiplication of SD and default value, the moments which
% are 0 will stay zero.
for iParam = 1 : length(musParam)
    curParam = musParam{iParam};
    for iMus = 1 : size(i_muscles, 1)
      curMusRight = i_muscles(iMus, 1);
      curMusLeft = i_muscles(iMus, 2);
      % Get mean value for distribution
      theMean = randMus.(curParam)(curMusRight);
      % Set value for right muscle
      randMus.(curParam)(curMusRight) = normrnd(theMean, SDmusParam*abs(theMean));
      % Set value for left muscle
      if any(ismember(musParamSym(:, 1), curParam)) % parameter of right side is listes
          % Assign the value of the right leg to the corresponding parameter of the left leg
          curParamLeft = musParamSym{ismember(musParamSym(:, 1), curParam), 2};
          randMus.(curParamLeft)(curMusLeft) = randMus.(curParam)(curMusRight);
      else
          % Simply copy the value to the same parameter
          randMus.(curParam)(curMusLeft) = randMus.(curParam)(curMusRight);
      end
   end
end

% Change seeSlack to keep the difference between seeSlack and L0 constant
L0Diff = randMus.L0 - obj.muscles.L0;
randMus.seeSlack = randMus.seeSlack + L0Diff;

% Assign the muscle table to the random participant object
randModel.muscles = randMus;

end
