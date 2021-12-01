function [outputArg1,outputArg2,outputArg3] = opt_normalization_1st2ndstage_exnan(second_state, action, goals, ind)
% make exceptions for nan
% context 1,2 : specific
% context 1,2 : Low high
% context 3,4 : flexible
% context 3,4 : High low
token_value = [0 0 0 0 0 40 20 10 0];

switch ind
    case 1 % specific low
        token_value = zeros(1,9);
        token_value(goals) = 1;
    case 2 % specific high
        token_value = zeros(1,9);
        token_value(goals) = 1;
    case 3 % flexible low
        token_value = token_value/max(token_value);
    case 4 % flexible high
        token_value = token_value/max(token_value);        
end


expect_value = []; % expected value for all actions, regarding token value. L1, R1, L2, R2 respectively
if second_state == 1
    %action 1 : Left // 2 : Right
    switch ind
        case 1
            state2 = mean([0.9*token_value(7)+0.1*token_value(8) 0.9*token_value(8)+0.1*token_value(9)]);
            state3 = mean([0.9*token_value(8)+0.1*token_value(9) 0.9*token_value(7)+0.1*token_value(9)]);
            state4 = mean([0.9*token_value(7)+0.1*token_value(6) 0.9*token_value(6)+0.1*token_value(9)]);
            state5 = mean([0.9*token_value(7)+0.1*token_value(9) 0.9*token_value(9)+0.1*token_value(6)]);
            expect_value = [0.9*state2 + 0.1*state3, 0.9*state4 + 0.1*state5];
        case 2
            state2 = mean([0.5*token_value(7)+0.5*token_value(8) 0.5*token_value(8)+0.5*token_value(9)]);
            state3 = mean([0.5*token_value(8)+0.5*token_value(9) 0.5*token_value(7)+0.5*token_value(9)]);
            state4 = mean([0.5*token_value(7)+0.5*token_value(6) 0.5*token_value(6)+0.5*token_value(9)]);
            state5 = mean([0.5*token_value(7)+0.5*token_value(9) 0.5*token_value(9)+0.5*token_value(6)]);
            expect_value = [0.5*state2 + 0.5*state3, 0.5*state4 + 0.5*state5];
        case 4
            state2 = mean([0.9*token_value(7)+0.1*token_value(8) 0.9*token_value(8)+0.1*token_value(9)]);
            state3 = mean([0.9*token_value(8)+0.1*token_value(9) 0.9*token_value(7)+0.1*token_value(9)]);
            state4 = mean([0.9*token_value(7)+0.1*token_value(6) 0.9*token_value(6)+0.1*token_value(9)]);
            state5 = mean([0.9*token_value(7)+0.1*token_value(9) 0.9*token_value(9)+0.1*token_value(6)]);
            expect_value = [0.9*state2 + 0.1*state3, 0.9*state4 + 0.1*state5];
        case 3
            state2 = mean([0.5*token_value(7)+0.5*token_value(8) 0.5*token_value(8)+0.5*token_value(9)]);
            state3 = mean([0.5*token_value(8)+0.5*token_value(9) 0.5*token_value(7)+0.5*token_value(9)]);
            state4 = mean([0.5*token_value(7)+0.5*token_value(6) 0.5*token_value(6)+0.5*token_value(9)]);
            state5 = mean([0.5*token_value(7)+0.5*token_value(9) 0.5*token_value(9)+0.5*token_value(6)]);
            expect_value = [0.5*state2 + 0.5*state3, 0.5*state4 + 0.5*state5];
    end
    
    likelihood = expect_value ./ sum(expect_value);
    outputArg1 = likelihood(action);
    outputArg2 = find(likelihood == max(likelihood));
    
    
else
    switch ind
        case 1 % specific low
            switch second_state
                case 2
                    expect_value = [0.9*token_value(7)+0.1*token_value(8) 0.9*token_value(8)+0.1*token_value(9)];
                case 3
                    expect_value = [0.9*token_value(8)+0.1*token_value(9) 0.9*token_value(7)+0.1*token_value(9)];
                case 4
                    expect_value = [0.9*token_value(7)+0.1*token_value(6) 0.9*token_value(6)+0.1*token_value(9)];
                case 5
                    expect_value = [0.9*token_value(7)+0.1*token_value(9) 0.9*token_value(9)+0.1*token_value(6)];
            end            
        case 2 %HL
            switch second_state
                case 2
                    expect_value = [0.5*token_value(7)+0.5*token_value(8) 0.5*token_value(8)+0.5*token_value(9)];
                case 3
                    expect_value = [0.5*token_value(8)+0.5*token_value(9) 0.5*token_value(7)+0.5*token_value(9)];
                case 4
                    expect_value = [0.5*token_value(7)+0.5*token_value(6) 0.5*token_value(6)+0.5*token_value(9)];
                case 5
                    expect_value = [0.5*token_value(7)+0.5*token_value(9) 0.5*token_value(9)+0.5*token_value(6)];
            end
            
        case 4 % LH
            switch second_state
                case 2
                    expect_value = [0.9*token_value(7)+0.1*token_value(8) 0.9*token_value(8)+0.1*token_value(9)];
                case 3
                    expect_value = [0.9*token_value(8)+0.1*token_value(9) 0.9*token_value(7)+0.1*token_value(9)];
                case 4
                    expect_value = [0.9*token_value(7)+0.1*token_value(6) 0.9*token_value(6)+0.1*token_value(9)];
                case 5
                    expect_value = [0.9*token_value(7)+0.1*token_value(9) 0.9*token_value(9)+0.1*token_value(6)];
            end            
            
        case 3 %HH
            switch second_state
                case 2
                    expect_value = [0.5*token_value(7)+0.5*token_value(8) 0.5*token_value(8)+0.5*token_value(9)];
                case 3
                    expect_value = [0.5*token_value(8)+0.5*token_value(9) 0.5*token_value(7)+0.5*token_value(9)];
                case 4
                    expect_value = [0.5*token_value(7)+0.5*token_value(6) 0.5*token_value(6)+0.5*token_value(9)];
                case 5
                    expect_value = [0.5*token_value(7)+0.5*token_value(9) 0.5*token_value(9)+0.5*token_value(6)];
            end
            
    end
    
    likelihood = expect_value ./ sum(expect_value);
    
    outputArg1 = likelihood(action);
    outputArg2 = find(likelihood == max(likelihood));
    
    if (sum(expect_value) == 0)
        outputArg1 = ones(size(expect_value));
        outputArg1 = outputArg1/sum(outputArg1);
%         idx = find(outputArg1 == max(outputArg1));
        outputArg2 = find(outputArg1 == max(outputArg1));
        outputArg1 = outputArg1(action);
%         outputArg2 = idx(randi(length(idx)));
    end
    
end

end

