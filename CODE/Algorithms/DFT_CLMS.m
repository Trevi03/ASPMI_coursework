function [w,error] = DFT_CLMS(x_input, z_input, mu, order, leakage)
    arguments
        x_input {mustBeNumeric}
        z_input {mustBeNumeric,mustBeEqualLength(x_input,z_input)}
        mu (1,1) {mustBeNumeric}
        order (1,1) {mustBeInteger}
        leakage (1,1) {mustBeNumeric} = 0;
    end

    N = length(x_input);
    w = complex(zeros(order,N));
    error = complex(zeros(1,N));
    x_hat = zeros(order,N);

    for k = 1: length(x_input)
        % calculate the prediction
        x_hat(k) = w(:,k)'*x_input(:,k);
        error(k) = z_input(k)-x_hat(k);

        % update
        w(:,k+1) = (1 - mu * leakage)* w(:,k) + mu*conj(error(k))*x_input(:,k);
    end
    w =  w(:,2:end);    
end


% Custom validation function
function mustBeEqualLength(a,b)
    % Test for equal size
    if ~isequal(size(a,2),size(b,2))
        eid = 'Size:notEqual';
        msg = 'Size of first input must equal size of second input.';
        error(eid,msg)
    end
end