%======================================================================
%> @file myvecnorm.m
%> @brief Function to compute the norm of a vector
%> @details
%> Details: myvecnorm()
%>
%> @author Marlies Nitschke
%> @date July, 2019
%======================================================================

%======================================================================
%> @brief Function to compute the norm of a vector
%>
%> @details
%> vecnorm() was introduced in Matlab 2017b. Since we are also working with
%> Matlab 2017, I reimplemented the functionality. 
%>
%> @param  A     Double matrix: Norm will be computed along the single columns (M x N)
%> @retval no    Double vector: Resulting norm (N x 1)
%======================================================================
function no = myvecnorm(A)

no = nan(size(A, 2), 1);
for iCol = 1 : size(A, 2)
    no(iCol) = norm(A(:, iCol));
end

end
