function [h,g,error] = ACLMS(x_input, z_input, mu, order)
    arguments
        x_input {mustBeVector}
        z_input {mustBeVector,mustBeEqualSize(x_input,z_input)}
        mu (1,1) {mustBeNumeric}
        order (1,1) {mustBeInteger}
    end

    N = length(x_input);
    h = complex(zeros(order,N));
    g = h;
    error = complex(zeros(1,N));
    x_hat = complex(zeros(1,N));

    delayed_x = zeros(N+order-1,1);
    delayed_x(order:N+order-1) = x_input;
    
    for k = 1: length(x_input)
        % calculate the prediction
        x_hat(k) = h(:,k)'*delayed_x(k:k+order-1) + g(:,k)'*conj(delayed_x(k:k+order-1));
        error(k) = z_input(k)-x_hat(k);

        % update
        h(:,k+1) = h(:,k) + mu*conj(error(k))*delayed_x(k:k+order-1);
        g(:,k+1) = g(:,k) + mu*conj(error(k))*conj(delayed_x(k:k+order-1));
    end
    h = h(:,2:end);   
    g = g(:,2:end);
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