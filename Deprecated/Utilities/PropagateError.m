function var_f = PropagateError(A,B,var_A,var_B,var_AB,combinationType)
%Function for propagating the error of the nonlinear combination of A and B.

%Inputs:
%   A: mean value of A
%   B: mean value of B
%   var_A: squared uncertainty in A
%   var_B: squared uncertainty in B
%   var_AB: joint squared uncertainty in A and B
%   combinationType: function f(A,B)
%       'multiply': f = AB
%       'divide': f = A/B

%Outputs:
%   var_f: squared uncertainty of combination

if strcmp(combinationType,'multiply')
    f = @(x,y) x*y;
    var_f = f(A,B).^2 .* ((var_A./A).^2 + (var_B./B).^2 + 2.*var_AB./(A.*B));
elseif strcmp(combinationType,'divide')
    f = @(x,y) x./y;
    var_f = f(A,B).^2 .* ((var_A./A).^2 + (var_B./B).^2 - 2.*var_AB./(A.*B));
else
    error('Please input a correct function f = f(A,B)');
end

end