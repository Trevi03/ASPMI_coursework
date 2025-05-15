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
    x_hat = zeros(order,N);
    delayed_x = zeros(order,N);
    
    % creating two shifted vectors by i, i.e. the order length
    for i = 1: order
        delayed_x(i,:) = [zeros(1,i), x_input(1:N-i)]; 
    end
    
    for k = 1: length(x_input)
        % calculate the prediction
        x_hat(k) = h(:,k)'*delayed_x(:,k)+ g(:,k)'*conj(delayed_x(:,k));
        error(k) = z_input(k)-x_hat(k);

        % update
        h(:,k+1) = h(:,k) + mu*conj(error(k))*delayed_x(:,k);
        g(:,k+1) = g(:,k) + mu*conj(error(k))*conj(delayed_x(:,k));
    end
    h =  h(:,2:end);   
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