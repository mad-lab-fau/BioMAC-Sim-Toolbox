%======================================================================
%> @file @Gait3d/createRandomPerson.m
%> @brief Gait3d function to create random musculoskeletal parameters
%> @details
%> Details: Gait3d::createRandomPerson()
%>
%> @author Anne Koelewijn
%> @date December, 2020
%======================================================================

%======================================================================
%> @brief Function to make MEX functions
%> @static
%> @public
%>
%> @details
%> Creates a musculoskeletal model with random parameters
%>
%> @param  symmusrand             double matrix: matrix with random numbers
%>                                that will be used to randomze muscle parameters
%> @param  massratio (optional)   ratio to be used for randomizing the mass of the body
%> @param  lengthratio (optional) ratio to be used for randomizing the length of the body
%======================================================================
function randModel = createRandomPerson(obj, symmusrand, massratio,lengthratio)

% Make a new object using the opensim, bodyheight and bodymass from obj
objClass = class(obj); % Get the class such that we can use this function also for Gait2dc_Exo
randModel = feval(objClass, obj.osim.file); % Call the constructor of the specific class

if nargin == 4
    % also randomize the height and weight
    %         %> @todo test randomPerson generator
    randModel.bodymass = randModel.bodymass*massratio;
    randModel.bodyheight = randModel.bodyheight*lengthratio;

    for i = 1:randModel.nSegments
        randModel.segments.mass(i) = randModel.segments.mass(i)*massratio;
    end
    
    for i = 1:randModel.nSegments
        randModel.segments.mass_center(i) = randModel.segments.mass_center(i)*lengthratio;
    end
    
    %loop through joints, except first and do the correct thing at
    %knee
    for i = 1:randModel.nJoints
        if or(i==3,i==8) %todo: get through name
            randModel.joints.t1_coefs(i,end) = randModel.joints.t1_coefs(i,end)*lengthratio;
            randModel.joints.t2_coefs(i,end) = randModel.joints.t2_coefs(i,end)*lengthratio;
        end
        randModel.joints.location(i,:) = randModel.joints.location(i,:)*lengthratio;
    end

    % Vary the inertia
    for i = 1:randModel.nSegments
        randModel.segments.inertia{i} = randModel.segments.inertia{i}*massratio*lengthratio^2; 
    end
end

randparams = {'fmax', 'lceopt', 'vmax'};%, 'kPEE', 'umax', 'arel'};

% Vary muscle parameters, symmetric between the legs
%first adjust fmax to mass
randModel.muscles.fmax = randModel.muscles.fmax*massratio;
rndmus = ones(randModel.nMus, 1)+0.1*(symmusrand);
muslen = zeros(randModel.nMus, 1);
for i = 1:randModel.nMus
    muslen(i) = randModel.muscles.lceopt(i) + randModel.muscles.lslack(i);
end
for j = 1:length(randparams)
    musparam = randModel.muscles.(randparams{j});
    for i = 1:randModel.nMus
        musparam(i) = musparam(i)*rndmus(i,j);
    end
    randModel.muscles.(randparams{j}) = musparam;   
end
for i = 1:randModel.nMus
    randModel.muscles.lslack(i) = muslen(i) - randModel.muscles.lceopt(i);
end
    
%     %we want to vary the width
%     width    = sqrt(randModel.muscles.kactive(i))*rndmus(i,3);  % Width of hill curve
%     randModel.muscles.kactive(i) = width^2;
end
