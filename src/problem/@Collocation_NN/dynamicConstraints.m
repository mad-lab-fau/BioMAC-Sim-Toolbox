%======================================================================
%> @file dynamicConstraints.m
%> @brief Collocation function to compute dynamic constraint
%> @details
%> Details: Collocation::dynamicConstraints()
%>
%> @author Eva Dorschky
%> @date November, 2017
%======================================================================

%======================================================================
%> @brief Computes constaint violation demanding dynamic equilibrium
%>
%> @param obj           Collocation class object
%> @param option        String parsing the demanded output
%> @param X             Double array: State vector containing at least 'states' and 'controls' of
%>                      the model and speed and duration of the periodic movement
%======================================================================
function output = dynamicConstraints(obj,option,X)
    %% check input parameter
    if ~isfield(obj.idx,'states') || ~isfield(obj.idx,'controls') || ~isfield(obj.idx,'dur') % check whether controls are stored in X
        error('Model states and controls and duration need to be stored in state vector X.')
    end
    
    %% compute demanded output
    h = X(obj.idx.dur)/obj.nNodes;
    nconstraintspernode = obj.model.nConstraints;

    % displacement and velocity array
    x2_first = X(obj.idx.states(:,2));
    [~,~,rfoot,lfoot] = gait2dc('Stick',x2_first);
    px = [rfoot(2:end-1,1);lfoot(2:end-1,1)];
    Ncontactpoints = size(px,1);
    
    d_array = zeros(obj.nNodes, Ncontactpoints);
    dd_array = zeros(obj.nNodes, Ncontactpoints);

    for iNode = 1:obj.nNodes
        x1 = X(obj.idx.states(:,iNode));
        x2 = X(obj.idx.states(:,iNode+1));
        xd =(x2-x1)/h;

        if strcmp(obj.Euler,'BE')
            [d, ddot] = obj.model.getDisplacement_Velocity(x2, xd);
%                 for iNcon = 1:Ncontactpoints
%                     d_array(iNode, iNcon) = d(iNcon);
%                     dd_array(iNode, iNcon) = ddot(iNcon);
%                     yc_matrix(iNode, iNcon) = yc_array(iNcon);
%                 end
        elseif strcmp(obj.Euler,'ME')
            [d, ddot] = obj.model.getDisplacement_Velocity((x1+x2)/2,xd);
        end
        for iNcon = 1:Ncontactpoints
            d_array(iNode, iNcon) = d(iNcon);
            dd_array(iNode, iNcon) = ddot(iNcon);
        end
    end

    % TODO: somehow apply to dfdx and so on.
    % TODO: get prediction before dynmaics and then put it into the
    % getDynamics method. probably redo the function smh.

    % should look like [[d,d,d,d],[d,d,d,d], ...]
    input = d_array;
    input(end,end,2) = 0;
    input(:,:,2) = dd_array;
    % input should look like [[[disp, vel],[disp,vel]...]...]
    input_py = obj.np.array(input);

    % Fy should then also look like [[f,f,f,f],[f,f,f,f],....]
    % times a 1000 because of normalization used.
    Fy = input(:,:,1).^2+3*input(:,:,2)-6;%double(obj.net.predict(input_py).squeeze()) * 1000; 

    hh = 1e-6; % You should probably play a bit with this parameter. It is supposed to be quite small with respect to the value of d.
%     dnew = d_array;
    dnew = zeros(obj.nNodes, Ncontactpoints, Ncontactpoints);
    dvnew = zeros(size(dnew));
    for iNcon = 1:Ncontactpoints
        dnew(:,:,iNcon) = d_array;
        dnew(:, iNcon,iNcon) = d_array(:, iNcon) + hh;
        dvnew(:,:,iNcon) = dd_array;
        dvnew(:, iNcon,iNcon) = dd_array(:, iNcon);        
    end

    input2 = reshape(dnew, obj.nNodes, []);
    input2(end,end,2) = 0;
    input2(:,:,2) = reshape(dvnew, obj.nNodes, []);
    input_py = obj.np.array(input2);

    Fynew = input2(:,:,1).^2+3*input2(:,:,2)-6;%double(obj.net.predict(input_py)) * 1000;

    dfdd = zeros(size(d_array));
    [~,c] = size(Fynew);
    [~,c_Fy] = size(Fy);

    j=c_Fy;
    k=1;
    for i = 1:c_Fy:c
%         dfdd(:, i:j) = (Fynew(:,i:j) - Fy) / hh;
        test = (Fynew(:,i:j) - Fy) / hh;
        dfdd(:, k) = test(:,k);
        j = j + c_Fy;
        k=k+1;
    end
   
    dnew = zeros(obj.nNodes, Ncontactpoints, Ncontactpoints);
    dvnew = zeros(size(dnew));
    for iNcon = 1:Ncontactpoints
        dnew(:,:,iNcon) = d_array;
        dnew(:, iNcon,iNcon) = d_array(:, iNcon);
        dvnew(:,:,iNcon) = dd_array;
        dvnew(:, iNcon,iNcon) = dd_array(:, iNcon) + hh;        
    end

    input3 = reshape(dnew, obj.nNodes, []);
    input3(end,end,2) = 0;
    input3(:,:,2) = reshape(dvnew, obj.nNodes, []);

    % input should look like [[[disp, vel_new],[disp,vel_new]...]...]
    input2_py = obj.np.array(input3);
    Fynew1 = input3(:,:,1).^2+3*input3(:,:,2)-6;%double(obj.net.predict(input2_py).squeeze()) * 1000; 
    dfdvel = zeros(size(dd_array));
    [~,c] = size(Fynew);
    j=c_Fy;
    k=1;
    for i = 1:c_Fy:c
%         dfdvel(:, i:j) = (Fynew1(:,i:j) - Fy) / hh;
        test = (Fynew1(:,i:j) - Fy) / hh;
        dfdvel(:, k) = test(:, k);
        j = j + c_Fy;
        k=k+1;
    end

    if strcmp(option,'confun')
        output = zeros(nconstraintspernode*obj.nNodes,1);
        
        % dynamic equations must be zero
        for iNode=1:obj.nNodes
            ic = (1:nconstraintspernode) +  (iNode-1)*nconstraintspernode; %indices of constraints of iNode in c
            x1 = X(obj.idx.states(:,iNode));
            x2 = X(obj.idx.states(:,iNode+1));
            xd =(x2-x1)/h;
            u2 = X(obj.idx.controls(:,iNode+1));
            
            if strcmp(obj.Euler,'BE')
                output(ic) = obj.model.getDynamics(x2,xd,u2, d_array(iNode,:), dfdd(iNode,:), dfdvel(iNode,:), Fy(iNode,:));	% backward Euler discretization
            elseif strcmp(obj.Euler,'ME')
                % we're using u2 instead of (u1+u2)/2 because u is
                % open loop and it makes no difference except u2
                % will converge more easily
                output(ic) = obj.model.getDynamics((x1+x2)/2,xd,u2, d_array(iNode,:), dfdd(iNode,:), dfdvel(iNode,:), Fy(iNode,:));
            end
        end
    elseif strcmp(option,'jacobian')
        output = spalloc(nconstraintspernode*obj.nNodes,length(X),obj.Jnnz);
        
        for iNode = 1:obj.nNodes
            ic = (1:nconstraintspernode) +  (iNode-1)*nconstraintspernode; %indices of constraints of iNode in c
            
            x1 = X(obj.idx.states(:,iNode));
            ix = obj.idx.states(:,iNode);
            x2 = X(obj.idx.states(:,iNode+1));
            ix2 = obj.idx.states(:,iNode+1);
            xd =(x2-x1)/h;
            u2 = X(obj.idx.controls(:,iNode+1));
            iu2 = obj.idx.controls(:,iNode+1);
            
            if strcmp(obj.Euler,'BE')
                [~, dfdx, dfdxdot, dfdu] = obj.model.getDynamics(x2,xd,u2, d_array(iNode,:), dfdd(iNode,:), dfdvel(iNode,:), Fy(iNode,:));
                output(ic,ix) = -dfdxdot'/h;
                output(ic,ix2) = dfdx' + dfdxdot'/h;
            elseif strcmp(obj.Euler,'ME')
                [~, dfdx, dfdxdot, dfdu] = obj.model.getDynamics((x1+x2)/2,xd,u2, d_array(iNode,:), dfdd(iNode,:), dfdvel(iNode,:), Fy(iNode,:));
                output(ic,ix) = dfdx'/2 - dfdxdot'/h;
                output(ic,ix2) = dfdx'/2 + dfdxdot'/h;
            end
            output(ic,iu2) = dfdu';
            
            % derivative of constraints with respect to duration (because h is duration/N)
            output(ic,obj.idx.dur) = -dfdxdot' * (x2-x1) / h^2 / obj.nNodes;
            
        end
    else
        error('Unknown option.');
    end    
end