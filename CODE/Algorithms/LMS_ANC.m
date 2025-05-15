function [y_hat, x_hat, w] = LMS_ANC(x_input,z_input,mu,order,leakage)
% The least mean square (LMS) algorithm with optional gear shift
%   Inputs:  x_input = input signal x
%            z_input = real output signal z
%            mu_adaptgain = the adaptation gain 
%            order = order of the adaptive filter (Nw+1) 
%   Outputs: y_hat = the LMS estimate 
%            error = output the error vector e[n] 
%            x_hat = error
    
    arguments
        x_input {mustBeVector}
        z_input {mustBeVector,mustBeEqualSize(x_input,z_input)}
        mu (1,1) {mustBeNumeric}
        order (1,1) {mustBeInteger}
        leakage (1,1) double = 0;
    end

    N = size(x_input,2); % length(x)=length(z)
    w = zeros(order,N+1); % weights
    y_hat = zeros(N,1);
    x_hat = y_hat;
    delayed_z = zeros(order, N);

    for i = 1:order
        delayed_z(i, :) = [zeros(1, i - 1), z_input(1:N - (i - 1))];
    end

    for i = 1 : N 
        y_hat(i) = w(:, i)' * delayed_z(:, i);     % y_hat = weights * input;
        x_hat(i) = x_input(i) - y_hat(i) ;      % e[n] = z[n] − y_hat[n]

        % w_new = w(1-mu) + mu * x_n * error;
        w(:, i + 1) = (1 - mu * leakage) .* w(:, i) + mu * x_hat(i) * delayed_z(:, i);
    end 
    w = w(:, 2:end);
end

% Custom validation function
function mustBeEqualSize(a,b)
    % Test for equal size
    if ~isequal(size(a),size(b))
        eid = 'Size:notEqual';
        msg = 'Size of first input must equal size of second input.';
        error(eid,msg)
    end
end