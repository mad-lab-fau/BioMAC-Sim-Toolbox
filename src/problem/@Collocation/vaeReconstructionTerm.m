function output = vaeReconstructionTerm(obj, option, X, vaeParams)

fctname = 'vaeReconstructionTerm';

if strcmp(option,'init')
    
    if ~isfield(obj.idx,'states')
        error('Model states are not stored in state vector X.')
    end
    obj.objectiveInit.(fctname).pyenv = pyenv(Version="/home/rzlin/ri94mihu/phd/BiomechPriorVAE/.venv/bin/python3.12", ExecutionMode="OutOfProcess");
    insert(py.sys.path, int64(0), what("BiomechPriorVAE").path);
    insert(py.sys.path, int64(0), strcat(what("BiomechPriorVAE").path,filesep,"src"));

    obj.objectiveInit.(fctname).idxJointsAllNodes = obj.idx.states(:, 1:obj.nNodes);
    obj.objectiveInit.(fctname).vaeParams = vaeParams;
    
    try
        %py.sys.path().append(vaeParams.pythonPath);
        obj.objectiveInit.(fctname).vaeModule = py.importlib.import_module('vaemodel');
        
        success = obj.objectiveInit.(fctname).vaeModule.initialize_vae(...
            vaeParams.modelPath, ...
            vaeParams.scalerPath, ...
            pyargs('num_dofs', int32(vaeParams.numDofs), ...
                   'latent_dim', int32(vaeParams.latentDim), ...
                   'hidden_dim', int32(vaeParams.hiddenDim), ...
                   'device', vaeParams.device));
        
        if ~success
            error('Failed to initialize VAE model');
        end
        
    catch ME
        error('Failed to initialize Python VAE interface: %s', ME.message);
    end
    
    obj.objectiveInit.(fctname).nJoints = length(obj.model.extractState('q'));
    
    output = NaN;
    return;
end


idxJointsAllNodes = obj.objectiveInit.(fctname).idxJointsAllNodes;
vaeModule = obj.objectiveInit.(fctname).vaeModule;
end_ = obj.model.nStates;
if obj.model.nCPs == 16
    currentIdx = idxJointsAllNodes([7:33 40:66 251:6:298 252:6:298 253:6:298 299:6:end_ 300:6:end_ 301:6:end_], 1:obj.nNodes);
elseif obj.model.nCPs == 12
    currentIdx = idxJointsAllNodes([7:33 40:66 251:6:286 252:6:286 253:6:286 287:6:end_ 288:6:end_ 289:6:end_], 1:obj.nNodes);
else
    error('Number of contact points is not supported.')
end
currentJoints = X(currentIdx);


x_ = X(obj.idx.states);
M = zeros(obj.nNodes, 33);

% Hardcoding mass and height here
mom_scale_fac = 1/9.81/1.83/75.16; % Hamner default values
currentJoints = [currentJoints; M'*mom_scale_fac]; % We don't use moments, but still fill in the zeros for Python


if strcmp(option,'objval')
    pyJoints = py.numpy.array(currentJoints);
    pyResult = vaeModule.reconstruct(pyJoints);
    nodeError = double(pyResult);
    output = nodeError;

    
elseif strcmp(option,'gradient')
    output = zeros(size(X));
    pyJoints = py.numpy.array(currentJoints);
    pyResult = vaeModule.reconstruct_withgrad(pyJoints);
    gradient = double(pyResult);
    if obj.model.nCPs == 16
        dVAEdx = gradient(:,1:102);
        dVEAdM = gradient(:,103:end);
    elseif obj.model.nCPs == 12
        dVAEdx = gradient(:,1:90);
        dVEAdM = gradient(:,91:end);
    end
    output(currentIdx) = dVAEdx';
end

end