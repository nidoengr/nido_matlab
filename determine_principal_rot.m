function [flag,rotAng_rd,eigVecCorEigVal1] = determine_principal_rot(BN)
arguments
   BN (3,3) double = [[ 0.925417  0.0296956 0.377786 ]',...
                      [ 0.336824 -0.521281 -0.784102 ]',...
                      [ 0.173648  0.852869 -0.492404 ]'];
end

[eigVecs,eigVals] = eig(BN);

TOL = 1e-6;

eigValsReal = diag(real(eigVals));
eigValsImag = diag(imag(eigVals));

idxImagEq0 = find(eigValsImag==0);

isResultValid = ~isempty(idxImagEq0);

rotAng_rd = nan;
eigVecCorEigVal1 = nan(3,1);
flag = false;
if isResultValid
   
   isRealEigVal1 = any(abs( 1 - eigValsReal(idxImagEq0) ) < TOL); % Note I had to drop the tolerance
   
   if isRealEigVal1
      
      idxEigVal1 = idxImagEq0;%find(isEigVal1); % FIX ME! handle better for cases of more than 1

      eigVecCorEigVal1 = eigVecs(:,idxEigVal1)

      rotAng_rd = acos( 0.5 * ( BN(1,1) + BN(2,2) + BN(3,3) - 1 ) );

      eigVecCorEigVal1 = (1/(2*sin(rotAng_rd)))*[BN(2,3) - BN(3,2)
                                                 BN(3,1) - BN(1,3)
                                                 BN(1,2) - BN(2,1)];

      flag = true;

   end
else
   warning('All eigenvalues have imaginary parts.');
end