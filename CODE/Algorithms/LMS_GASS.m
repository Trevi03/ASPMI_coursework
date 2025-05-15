function [y_hat, error, w, mus] = LMS_GASS(x_input,z_input,mu,rho, alpha, order, algorithm)
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
        rho (1,1) {mustBeNumeric}
        alpha (1,1) {mustBeNumeric}
        order (1,1) {mustBeInteger}
        algorithm (1,1) {mustBeInteger, mustBeLessThanOrEqual(algorithm,5)} % original, Benveniste, Ang & Farhang, Matthews & Xie
    end
    N = size(x_input,2); % length(x)=length(z)
    w = zeros(order,N+1); % weights
    y_hat = zeros(1, N);
    error = zeros(1, N);
    delayed_x = zeros(order, N);
    delay = 1;
    psi = w;
    leakage = 0;

    if algorithm <3
        mus = mu*ones(1,length(x_input)+1);
    else
        mus = zeros(1,length(x_input)+1);
        mus(1) = mu;
    end

    for i = 1:order
        delayed_x(i, :) = [zeros(1, i + delay - 1), x_input(1:N - (i + delay - 1))];
    end

    for i = 1 : N 
        y_hat(:, i) = w(:, i)' * delayed_x(:, i);     % y_hat = weights * input;
        error(i) = z_input(i) - y_hat(i) ;      % e[n] = z[n] − y_hat[n]
        
        w(:, i + 1) = (1 - mu*leakage) * w(:, i) + mus(i) * error(i) * delayed_x(:, i);

        switch algorithm
            case 1
                % don't do anything
            case 2
                % don't do anything
            case 3 % benveniste
                % now adding an update for mu
                mus(i+1)= mus(i) + rho*error(i)*delayed_x(:,i)'*psi(:,i);
                psi(:,i+1) = (eye(order,order) - (mus(i)*delayed_x(:,i)'*(delayed_x(:,i))))*psi(:,i) + error(i)*delayed_x(:,i);
            case 4 % ang and farhang
                % now adding an update for mu
                mus(i+1)= mus(i) + rho*error(i)*delayed_x(:,i)'*psi(:,i);
                psi(:,i+1) = alpha*psi(:,i) + error(i)*delayed_x(:,i);
            case 5 % Matthews & Xie
                % now adding an update for mu
                mus(i+1)= mus(i) + rho*error(i)*delayed_x(:,i)'*psi(:,i);
                psi(:,i+1) = error(i)*delayed_x(:,i);
        end
    end 
    w = w(:, 2:end);
    mus = mus(:,2:end);
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