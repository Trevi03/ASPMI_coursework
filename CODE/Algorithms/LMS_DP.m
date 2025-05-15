function [y_hat, error, w, alphas] = LMS_DP(x_input,z_input,mu,order,w_init,bias,alpha,DC)
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
        z_input {mustBeVector,mustBeEqualSize(x_input,z_input)}
        mu (1,1) {mustBeNumeric}
        order (1,1) {mustBeInteger}
        w_init {mustBeVector}
        bias (1,1) logical
        alpha (1,1) double
        DC (1,1) {mustBeNumeric} = 1;
    end

    N = length(x_input); % length(x)=length(z)
    w = zeros(order+bias,N+1); % weights
    w(:,1) = w_init;
    y_hat = zeros(1, N);
    error = y_hat;

    alphas = zeros(1, N);

    delayed_x = zeros(order, N);
    for i = 1:order
        delayed_x(i, :) = [zeros(1, i), x_input(1:N -i)'];
    end

    if bias
        delayed_x = [DC*ones(1,N); delayed_x];
    end

    for i = 1 : N 
        y_hat(i) = alpha*tanh(w(:, i)' * delayed_x(:, i));
        error(i) = z_input(i) - y_hat(i);
        act_func = alpha*(1-(y_hat(:,i)/alpha)^2); % alpha[1 - tanh^2(w(:, i)' * delayed_x(:, i))]

        % w_new = w(1-mu) + mu * x_n * error;
        w(:, i + 1) = w(:, i) + (mu *act_func* error(i)) .* delayed_x(:, i);
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