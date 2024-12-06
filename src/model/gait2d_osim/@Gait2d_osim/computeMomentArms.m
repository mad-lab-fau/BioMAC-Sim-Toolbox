%======================================================================
%> @file @Gait2d_osim/computeMomentArms.m
%> @brief Gait2d_osim function to generate  muscle length polynomials from Opensim muscle paths
%> @details
%> Details: Gait2d_osim::computeMomentArms()
%> The code is copied from gait3d entirely, with adjustments for gait2d_osim
%>
%> @author Antonie J. (Ton) van den Bogert
%> @date October 26, 2015
%======================================================================

%======================================================================
%> @brief Function to generate  muscle length polynomials from Opensim muscle paths
%> @protected
%>
%> @details
%> For each muscle, the generalized coordinates are moved through their entire range, in
%> a specified number of steps (currently 6).  So for a muscle that spans 4 dofs, we will examine
%> 6^4 = 1296 poses.  For each pose, the Opensim API is used to obtain the moment arms of
%> the muscle.  Also the API is used to obtain the muscle-tendon length when all generalized
%> coordinates are zero.
%>
%> Stepwise regression is done on the moment arm data to produce the best set of polynomial terms.
%> The addition of terms stops when one of the following is true:
%> 1. RMS difference between moment arms from Opensim and polynomial, relative to the largest
%> moment arm, is less than "relerror" (currently 0.05)
%> 2. Adding another term would decrease the RMS difference by less than "minerrorchange" (currently
%> 0.01)
%> 3. The maximum number of terms (currently 40) has been reached.
%>
%> The examine option is intended for model developers to inspect the polynomials and make sure
%> they are good enough.  When regenerating polynomials for a model that has already been tested,
%> the examine and name options should not be used.
%>
%> Saves results at strrep(settings.osimfile, '.osim', '_momentarms.mat')
%>
%> @param obj                 Gait2d_osim class object
%> @param range_muscleMoment  Table: Ranges of DOFs which were used to fit polynomials
%> @param examine             (optional, default 0) when 1: plots results for each muscle and pauses for inspection
%> @param name                (optional) name of the muscle that should be processed (when missing, all are done)
%> @retval momentarm_model	  Struct: Resulting polynomials which are also saved at strrep(settings.osimfile, '.osim', '_momentarms.mat')
%======================================================================
function [momentarm_model] = computeMomentArms(obj, range_muscleMoment, examine, name)

try
    settings.osimfile = obj.osim.file;
catch
    error('No Opensim file stored. Please read Opensime file first calling Gait2d_osim.readOsim(filename).')
end

momentarmfile = strrep(which(settings.osimfile), '.osim', '_momentarms.mat');


warning('Opensim library let Matlab crash if IPOPT exits with an exception. Save Gait2d_osim model and open a new Malab session if running an optimization.');

fileID = fopen(which(obj.osim.file));
fileID_read =fread(fileID);
md = java.security.MessageDigest.getInstance('SHA-256');
osim_sha256 = typecast(md.digest(uint8(fileID_read))', 'uint8');


settings.steps = 25;				% the number of steps for the range of motion of each dof
settings.order = 4; 			% max. polynomial order
settings.relerror = 0.03;		% when relative error in moment arm is below this level, no more terms are added
settings.maxterms = 40;			% maximum number of polynomial terms
settings.minerrorchange = 0.01;	% when fit error improves less than this when a term is added, stop adding terms
if (nargin < 3)
    settings.examine = 0;
else
    settings.examine = examine;
end
%settings.examine = 1;
% do all muscles, unless a name was specified
fid = fopen('gait2d_osim_momentarms.log','w');
fprintf(fid,'muscle               MaxMomarm(mm) RMSfull(mm) Nterms RMSerr(mm) MAXerr(mm) mom.arm with max error    angles where max error occurs\n');
fprintf(fid,'-------------------- ------------- ----------- ------ ---------- ---------- ---------------------- --------------------------------\n');
fclose(fid);
muscle_names = obj.muscles.Properties.RowNames;
for imus = 1:obj.nMus
    mus = table2struct(obj.muscles(imus,:));
    mus.name = muscle_names{imus};
    if (nargin < 4) || strcmp(name,mus.name(imus))
        momentarm_model{imus} = muscle_poly(obj,mus,range_muscleMoment,settings);
        %         % obj.muscles{imus}.dof_rom = momentarm_model{imus}.rom; % this is the range of motion for each dof where the momentarm model is valid
        %         obj.muscles.lparam_count(imus) = momentarm_model{imus}.num_lparams; % the number of coeficients in the poly
        %         obj.muscles.lparams(imus) = momentarm_model{imus}.lparams; % the exponents for each coeficient - size is lparam_count x dof_count
        %         obj.muscles.lcoefs(imus) = momentarm_model{imus}.lcoef; % the array of coeficients - size is lparam_count x 1
    end
end

save(momentarmfile,'settings','momentarm_model','range_muscleMoment','osim_sha256');


%> @cond DO_NOT_DOCUMENT
%=========================================================
    function [poly] = muscle_poly(obj, mus, range_muscleMoment, settings)
        disp(['Getting Opensim moment arms for ', mus.name]);
        
        allvectxts = '';
        for idof=1:obj.nDofs
            if isempty(allvectxts)
                allvectxts = ['dof_values{' num2str(idof) '}'];
            else
                allvectxts = [allvectxts ',dof_values{' num2str(idof) '}'];
            end
        end
        
        % reset all dofs to the default value
        for idof=1:obj.nDofs
            dof_values{idof} = 0;
        end
        
        % for the dofs this muscle crosses, go through their range
        dof_range = zeros(mus.dof_count,2);
        for idof = 1:mus.dof_count
            dof = mus.dof_indexes(idof);
            range = range_muscleMoment{dof,:}; % with range_osim we do not get good polynomials
            values = linspace(range(1), range(2), settings.steps);
            
            dof_values{dof} = values;
            dof_range(idof,1:2) = range;
            disp(['  DOF: ', obj.dofs.Properties.RowNames(dof),' range:   ',num2str(range)]);
        end
        
        alljnts = nan; % Seems to be needed for Matlab 2021a
        eval(['alljnts = combvec(' allvectxts ')'';']);
        
        % allmomarms contains the moment arms of this muscle for all
        % combinations of joint angles in alljnts
        % rows: # combinations of angles
        % columns: # of dofs this muscle crosses
        allmomarms = zeros(size(alljnts,1),mus.dof_count);
        disp(['  number of poses: ',num2str(size(allmomarms))]);
        [allmomarms, mus_length] = opensim_get_momentarms(which(settings.osimfile), alljnts, mus.name, mus.dof_names);
        
        % now fit a polynomial to the Opensim results
        ndofs = mus.dof_count; % number of dofs spanned by this muscle
        num_data = size(allmomarms,1);
        
        % count how many parameters the polynomial model for muscle length has
        npar = prod(1:(ndofs+settings.order))/prod(1:ndofs)/prod(1:settings.order);
        fprintf(1,'Number of DOFs:   %d\n',ndofs);
        fprintf(1,'Polynomial order: %d\n',settings.order);
        fprintf(1,'Potential number of polynomial terms: %d\n',npar);
        tot_data = num_data*ndofs;	% total number of data points
        A = zeros(tot_data, npar);	% allocate memory space for A
        b = zeros(tot_data, 1);		% allocate memory space for b
        
        % get angle values for the dofs the muscle crosses
        musdof_indeces = mus.dof_indexes;
        ang = (alljnts(:,musdof_indeces) + 1e-6);	% protect against angle = 0.0
        
        for idof = 1:ndofs
            % read moment arm from allmomarms matrix and store in b (store dL/dq = -moment_arm)
            b((idof-1)*num_data+1:idof*num_data) = -allmomarms(:,idof);
            
            % generate the npar polynomial terms, and for each term, add a column to A
            polylist = zeros(npar,ndofs);
            expon = zeros(num_data,ndofs);	% start with all exponents zero
            for ii=1:npar
                polylist(ii,:) = expon(1,:);
                A((idof-1)*num_data+1:idof*num_data,ii) = expon(:,idof).*prod(ang.^expon,2)./ang(:,idof); % contribution of this term to moment arm idof
                % generate the next set of exponents, but remain within model order
                k = 1;
                while (1)
                    expon(:,k) = expon(:,k)+1;
                    if (sum(expon(1,:)) > settings.order && ii<npar)
                        expon(:,k)=0;
                        k = k+1;
                    else
                        break;
                    end
                end
            end     % done generating model terms
        end		% done reading all data files for this muscle
        fprintf('Total number of data points: %d\n',tot_data);
        
        % now we have all data for this muscle stored in A and b
        % remove from A the first column (constant term of polynomial)
        A = A(:,2:npar);
        polylist = polylist(2:npar,:);
        
        maxmomarm = max(abs(b))*1000; 	% maximum moment arm (in mm) for the muscle and dof
        
        % solve the full model with npar-1 terms
        p = A\b;		% compute coefficients of the best fitting model
        bpred = A*p;	% these are the moment arms predicted by this model
        res = bpred-b;	% residuals
        RMSfull = (sqrt(sum(res.^2)/tot_data)) * 1000;		% RMS of residuals, in mm
        fprintf('RMS fit error of the full model is: %f mm\n',RMSfull);
        fprintf('Maximum moment arm: %f mm\n',maxmomarm);
        
        % now do stepwise regression to select polynomial terms for a smaller model
        Aselected = [];
        polylist_selected = [];
        npar_selected = 0;
        % outer loop: successively add columns to Aselected
        for i = 1:npar-1
            [~, ncolumns] = size(A);
            % inner loop: find the column of A that causes most reduction in RMS when added to Aselected
            RMSnew = zeros(1,ncolumns); % this will store the RMS errors of each expanded model
            for j = 1:ncolumns
                % add column j from A to Anew
                Anew = [Aselected A(:,j)];
                % solve new p's
                pnew = Anew\b;
                % compute new RMS fit error
                RMSnew(j) = (norm(Anew*pnew - b)/(sqrt(tot_data))) * 1000; 	% convert to mm
            end
            % now determine which expanded model had the lowest RMS
            [RMSmin, col] = min(RMSnew);
            % if the change in error is less than a certain percentace, stop without adding this term
            if ( (i>1) && ((RMS - RMSmin)/RMS < settings.minerrorchange) )
                fprintf('Change in error: %3f. No more terms added.\n ',(RMS - RMSmin)/RMS);
                break;
            end
            if ( i>settings.maxterms )
                fprintf('Maximum number of terms (%d) reached.\n ', maxterms);
                break;
            end
            % otherwise add this column to Aselected
            Aselected = [Aselected A(:,col)];
            npar_selected = npar_selected + 1;
            p = Aselected\b;		% solve the coefficients (again)
            % compute some error measures
            [MAXERR, MAXERRindex] = max(abs(Aselected*p-b));		% maximum moment arm error
            SSE = sum((Aselected*p-b).^2); 				% summed squared errors (SSE)
            RMS = (sqrt(SSE/tot_data)) * 1000; 			% RMS error, is the same as what we had before
            GCV = SSE/((num_data - npar_selected)^2);    	% Generalized Cross Validation
            AIC = log(SSE) + 2*npar_selected/tot_data;   	% Akaike's Information Criterion
            % print what we just did, on screen and on output file
            fprintf('%3i: Added term ',i);
            fprintf('%7.4f ',p(end));
            fprintf('%i ',polylist(col,:));
            fprintf('-- RMS=%6.2f, GCV=%8.3e, AIC = %6.2f\n',RMS,GCV,AIC);
            
            % remember the exponents of this polynomial term
            polylist_selected = [polylist_selected ; polylist(col,:)];
            % remove this column from A and from polylist so it is not used again
            A = [A(:,1:(col-1)) A(:,(col+1):ncolumns)];
            polylist = [polylist(1:(col-1),:);polylist((col+1):ncolumns,:)];
            % stop adding terms if RMS error in fvectors is less than a percentage of max
            if ( RMS <= settings.relerror*maxmomarm )
                fprintf('Relative error in moment arm is less than %6.3f\n', settings.relerror);
                break;
            end
        end 		% and go find the next term
        
        % plot muscle length from Opensim and polynomial
        if (settings.examine)
            examine_momarms(mus, ang, allmomarms, -Aselected*p);
        end
        
        % save the polynomial in structure
        poly.name = mus.name;
        poly.num_lparams = npar_selected+1;
        poly.lparams = [zeros(1,ndofs);polylist_selected]; % size of lparams: num_lparams x num_dofs
        poly.lcoef = [mus_length; p];
        % we should store dof_range, if we ever make this smaller than the passive range of motion
        
        % write something on the log file
        % header: 'muscle               MaxMomarm(mm) RMSfull(mm) Nterms RMSerr(mm) MAXerr(mm) mom.arm with max error angles\n'
        log_name = strrep(which(settings.osimfile), '.osim', '_momentarms.log');
        fid1 = fopen(log_name,'a');
        fprintf(fid1,'%20s ', mus.name);
        fprintf(fid1,'%13.2f ', maxmomarm);
        fprintf(fid1,'%11.2f ', RMSfull);
        fprintf(fid1,'%6d', numel(poly.lcoef));
        fprintf(fid1,'%11.2f ', RMS);
        fprintf(fid1,'%10.2f ', 1000*MAXERR);
        fprintf(fid1,'%22s ', mus.dof_names{floor((MAXERRindex-1)/num_data)+1});
        for i=1:ndofs
            fprintf(fid1,'%10.2f ', 180/pi*ang(mod((MAXERRindex-1),num_data)+1,i));
        end
        fprintf(fid1,'\n');
        fclose(fid1);
        
    end
%=========================================================
    function [momArms,zerolength] = opensim_get_momentarms(osimfile, angles, Mus, Dofs)
        % function minusdLdq = opensim_get_dLdq(Mod, angles, Mus, Dofs)
        %
        % This function calculates -dL/dq of muscle "Mus" about dof set
        % "Dofs" at a given angle matrix "angles" of opensim model "Mod"
        %
        % The momArm output is from Opensim's built-in moment arm calculation, which was incorrect
        % in some versions of Opensim, but is correct in Opensim 3.3.
        %
        % "Angles" can be a vector (one hand position) or a matrix (one hand
        % position per row)
        % "Dofs" can a string (one dof) or a cell array (multiple dofs)
        %
        % Adapted from opensim_get_momarm.m
        % Dimitra Blana, March 2012
        % Patrick Kugler, March 2013
        %
        % 28/3/2012: Use setValue (with active constraints) instead of set state
        %
        import org.opensim.modeling.*;
        % load the model
        Mod = Model(osimfile);
        
        % initialize the system and get the initial state
        state = Mod.initSystem();
        
        % equilibrate all muscles
        Mod.equilibrateMuscles(state);
        
        % get the coordinates structure
        CoordSet = Mod.getCoordinateSet();
        DofNames = ArrayStr();
        CoordSet.getNames(DofNames);
        nDofs = CoordSet.getSize;
        for idof = 1:nDofs
            DofNames_str{idof} = char(DofNames.getitem(idof-1));
        end
        
        % How many moment arms do we want?
        if iscell(Dofs)
            num_request_dofs = length(Dofs);
        else
            num_request_dofs = 1;
        end
        
        % find indeces of dofs for moment arms
        if num_request_dofs==1
            dofindex = find(strcmp(Dofs,DofNames_str));
        else
            for idof = 1:num_request_dofs
                dofindex(idof) = find(strcmp(Dofs{idof},DofNames_str));
            end
        end
        
        % get the muscles
        MuscleSet = Mod.getMuscles();
        MuscleNames = ArrayStr();
        MuscleSet.getNames(MuscleNames);
        nMus = MuscleNames.getSize;
        for iMus = 1:nMus
            MuscleNames_str{iMus} = char(MuscleNames.getitem(iMus-1));
        end
        musindex = find(strcmp(Mus,MuscleNames_str));
        
        % angles matrix: one position per row
        [nrows,ncols] = size(angles);
        
        % initialise matrix for storing moment arms
        momArms = zeros(size(angles,1),num_request_dofs);
        
        % get zerolength (muscle-tendon length when all q's are zero)
        for idof = 1:nDofs
            currentDof = CoordSet.get(idof-1);
            currentDof.setValue(state,0,1);
        end
        zerolength = MuscleSet.get(musindex-1).getLength(state);
        
        for istep = 1:size(angles,1)
            if ~mod(istep,100)
                disp(['Muscle ',char(MuscleSet.get(musindex-1).getName()), ' - pose ',...
                    num2str(istep),' of ', num2str(size(angles,1))]);
            end
            
            % set dof values for this step
            idof_counter = 0;
            for idof = 1:nDofs
                currentDof = CoordSet.get(idof-1);
                if currentDof.getLocked(state)
                    continue;
                end
                idof_counter = idof_counter + 1;
                currentDof.setValue(state,angles(istep,idof_counter),1);
            end
            
            for idof = 1:num_request_dofs
                
                % find moment arm = -dL/dq with a central difference approximation
                % this was needed in a previous version of Opensim when they had a bug
                % in their computeMomentArm code
                % 	currentDof = CoordSet.get(dofindex(idof)-1);
                % 	currentDof.setValue(state,angles(istep,dofindex(idof))-0.0001,1);
                % 	L1 = MuscleSet.get(musindex-1).getLength(state);
                % 	currentDof.setValue(state,angles(istep,dofindex(idof))+0.0001,1);
                % 	L2 = MuscleSet.get(musindex-1).getLength(state);
                % 	momarm(istep,idof) = -(L2-L1)/0.0002;
                % set dof to original value
                % 	currentDof.setValue(state,angles(istep,dofindex(idof)),1);
                
                momArms(istep,idof) = MuscleSet.get(musindex-1).computeMomentArm(state,CoordSet.get(dofindex(idof)-1));
            end
        end
    end
%=========================================================
    function y = combvec(varargin)
        %COMBVEC Create all combinations of vectors.
        %
        %  <a href="matlab:doc combvec">combvec</a>(A1,A2,...) takes any number of inputs A, where each Ai has
        %  Ni columns, and return a matrix of (N1*N2*...) column vectors, where
        %  the columns consist of all combinations found by combining one column
        %  vector from each Ai.
        %
        %  For instance, here the four combinations of two 2-column matrices are
        %  found.
        %
        %    a1 = [1 2 3; 4 5 6];
        %    a2 = [7 8; 9 10];
        %    a3 = <a href="matlab:doc combvec">combvec</a>(a1,a2)
        %
        % Mark Beale, 12-15-93
        % Copyright 1992-2010 The MathWorks, Inc.
        % $Revision: 1.1.10.2 $  $Date: 2010/04/24 18:07:52 $
        %
        if length(varargin) == 0
            y = [];
        else
            y = varargin{1};
            for i=2:length(varargin)
                z = varargin{i};
                y = [copy_blocked(y,size(z,2)); copy_interleaved(z,size(y,2))];
            end
        end
    end
%=========================================================
    function b = copy_blocked(m,n)
        
        [mr,mc] = size(m);
        b = zeros(mr,mc*n);
        ind = 1:mc;
        for i=[0:(n-1)]*mc
            b(:,ind+i) = m;
        end
    end
%=========================================================
    function b = copy_interleaved(m,n)
        
        [mr,mc] = size(m);
        b = zeros(mr*n,mc);
        ind = 1:mr;
        for i=[0:(n-1)]*mr
            b(ind+i,:) = m;
        end
        b = reshape(b,mr,n*mc);
    end
%==========================================================
    function examine_momarms(mus, angles, ma_osim, ma_poly);
        [nposes, ndof] = size(ma_osim);
        ma_poly = reshape(ma_poly, nposes, ndof);
        ndof = mus.dof_count;
        figure(1);clf;
        for i=1:ndof
            subplot(ndof,2,2*i-1);
            plot(angles(:,i).*180/pi, [ma_osim(:,i) ma_poly(:,i) ma_osim(:,i)-ma_poly(:,i)], '.-');
            if (i==ndof)
                legend('Opensim','polynomial','difference');
            end
            ylabel('moment arm (m)');
            title(strrep([mus.dof_names{i} ' moment arm for ' mus.name], '_', '\_'));
            subplot(ndof,2,2*i);
            plot(angles(:,i)*180/pi, '.-');
            ylabel('angle (deg)');
            title(strrep(mus.dof_names{i}, '_', '\_'));
        end
        disp('Hit ENTER to continue');
        pause;
    end
%> @endcond


end


