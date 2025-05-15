classdef FunctionCollection
    %ErrorCalculation collection of functions for error calculation
    properties
        Description
    end

    methods
        function obj = FunctionCollection()
            obj.Description = 'collection of functions for error calculation';
        end
        
        function error = SE(~,x1,x2)
            %Mean Squared Error
            error = (x1 - x2).^2;
        end

        function error = MSE(~,x1,x2)
            %Mean Squared Error
            error = mean((x1 - x2).^2,1);
        end

        function error = RMSE(~,x1,x2)
            %Root Mean Squared Error
            error = mean((x1 - x2).^2,1)^0.5;
        end

        function error = MAE(~,x1,x2)
            %Mean Absolute Error
            error = mean(abs(x1 - x2), 'all');
        end

        function t = calcTimeConst(u, l)
            t = -1 / log(abs(1 - 2 * u * l));
        end
    end
end