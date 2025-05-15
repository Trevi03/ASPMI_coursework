%%% Data generation for the 1D frequency decomposition CNN, ASPMI coursework
%%% Harry Davies and Danilo Mandic, 2024
%% 1D DCT - extended 8 point
N = 500;
Xtot = 10000;
sl_multiplier = 5;
t = ((2*pi)/N):((2*pi)/N):2*pi*sl_multiplier;
trainX = zeros(Xtot,N*sl_multiplier);
trainy = zeros(Xtot,8);
%disp('start')
for i = 1:Xtot
    disp(i) %to indicate progress if the script is running slowly
    noise = 0*randn(1)*randn(1,N*sl_multiplier);
    y_tot = randi([0,30],1);
    dct_ideal = zeros(1,500);
    if y_tot == 0
        y = ones(1,N*sl_multiplier);
        dct_ideal(1) = 1;
    else
        if y_tot > 2
            y_tot = 1; %majority of data generated with 1 sinusoid
        end
        y = zeros(1,N*sl_multiplier);
        for p = 1:y_tot 
            f2 = randi([1,7],1);
            if f2 == 31
                f2 = f2+randi([1,200],1);
            end

            dct_ideal(f2+1) = dct_ideal(f2+1) + 1;
            f = f2/2;
            time_add = (2*pi)/randi([1 50],1);
            y = y + cos(f*(t+time_add));
        end
        y = detrend(y);
    end
    y = y + noise;
    DCT_output = [];
    for u = 0:1:N-1
        if u == 0
            a = sqrt(1/N);
        else
            a = sqrt(2/N);
        end

        for x = 0:1:N-1
            temp_dct(x+1) = y(x+1)*cos((((2*x)+1)*u*pi)/(2*N));
        end
        DCT_output(u+1) = sum(temp_dct);

    end
    dct_ideal(1) = dct_ideal(1)*2;
    trainX(i,:) = y;
    trainy(i,:) = rescale(dct_ideal(1:8));
end

%%

save('example_training_n2500_1freq_rescale_k8.mat', 'trainX', 'trainy');