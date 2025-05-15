function [y_hat, error, w] = LMS_ALE(x_input,mu,order,leakage,delay)
% The least mean square (LMS) algorithm with optional gear shift
%   Inputs:  x_input = input signal x
%            z_input = real output signal z
%            mu_adaptgain = the adaptation gain 
%            order = order of the adaptive filter (Nw+1) 
%   Outputs: y_hat = the LMS estimate 
%            error = output the error vector e[n] 
%            evol_weights = evolution of the adaptive weights in time, (Nw + 1) × N matrix
    
    arguments
        x_input {mustBeVector}
        mu (1,1) {mustBeNumeric}
        order (1,1) {mustBeInteger}
        leakage (1,1) double = 0;
        delay (1,1) double = 1;
    end

    N = size(x_input,2); % length(x)=length(z)
    w = zeros(order,N+1); % weights
    y_hat = zeros(1, N);
    error = zeros(1, N);
    delayed_x = zeros(order, N);

    for i = 1:order
        delayed_x(i, :) = [zeros(1, i + delay - 1), x_input(1:N - (i + delay - 1))];
    end

    for i = 1 : N 
        y_hat(i) = w(:, i)' * delayed_x(:, i);     % y_hat = weights * input;
        error(i) = x_input(i) - y_hat(i) ;      % e[n] = z[n] − y_hat[n]

        % w_new = w(1-mu) + mu * x_n * error;
        w(:, i + 1) = (1 - mu * leakage) * w(:, i) + mu * error(:, i) * delayed_x(:, i);
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