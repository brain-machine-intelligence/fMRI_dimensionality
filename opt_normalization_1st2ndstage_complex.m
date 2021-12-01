function [outputArg1,outputArg2] = opt_normalization_1st2ndstage_complex(second_state, action, token_value, ind, condi)


expect_value = []; % expected value for all actions, regarding token value. L1, R1, L2, R2 respectively
if second_state == 1
    %action 1 : Left // 2 : Right
    switch ind
        case 1
            state2 = mean([0.9*0 + 0.1*token_value(3),  0.9*token_value(1) + 0.1*token_value(2)]);
            state3 = mean([0.9*token_value(1) + 0.1*0,  0.9*token_value(2) + 0.1*token_value(3)]);
            expect_value = [0.9*state2 + 0.1*state3, 0.9*state3 + 0.1*state2];
        case 2
            state2 = mean([0.9*0 + 0.1*token_value(3),  0.9*token_value(1) + 0.1*token_value(2), 0.9*token_value(3) + 0.1*0, 0.9*token_value(2) + 0.1*token_value(1)]);
            state3 = mean([0.9*token_value(1) + 0.1*0,  0.9*token_value(2) + 0.1*token_value(3), 0.9*0 + 0.1*token_value(1), 0.9*token_value(3) + 0.1*token_value(2)]);
            expect_value = [0.9*state2 + 0.1*state3, 0.9*state3 + 0.1*state2];
        case 3
            state2 = mean([0.5*0 + 0.5*token_value(3),  0.5*token_value(1) + 0.5*token_value(2)]);
            state3 = mean([0.5*token_value(1) + 0.5*0,  0.5*token_value(2) + 0.5*token_value(3)]);
            expect_value = [0.5*state2 + 0.5*state3, 0.5*state3 + 0.5*state2];
        case 4
            state2 = mean([0.5*0 + 0.5*token_value(3),  0.5*token_value(1) + 0.5*token_value(2), 0.5*token_value(3) + 0.5*0, 0.5*token_value(2) + 0.5*token_value(1)]);
            state3 = mean([0.5*token_value(1) + 0.5*0,  0.5*token_value(2) + 0.5*token_value(3), 0.5*0 + 0.5*token_value(1), 0.5*token_value(3) + 0.5*token_value(2)]);
            expect_value = [0.5*state2 + 0.5*state3, 0.5*state3 + 0.5*state2];    
    end
    
    likelihood = expect_value ./ sum(expect_value);
    
    if strcmp (condi, 'likelihood')
        outputArg1 = likelihood(action);
        outputArg2 = find(likelihood == max(likelihood));
    else
        outputArg1 = ismember(action, find(likelihood == max(likelihood)));
    end
    
    
else
    switch ind
        case 1 % LL
            %                 opt(1,2) = opt(1,2)+1;
            if second_state == 2 % a1 should be L1. L1 : NOTHING, L2 : BLUE, R2 : RED, R1 : SILVER
                expect_value = [0.9*0 + 0.1*token_value(3),  0.9*token_value(1) + 0.1*token_value(2)];
            else % second_state == 3 L1 : SILVER, L2 : NOTHING, R2 : BLUE, R1 : RED
                expect_value = [0.9*token_value(1) + 0.1*0,  0.9*token_value(2) + 0.1*token_value(3)];
            end
            
        case 2 %HL
            %                 opt(2,2) = opt(2,2)+1;
            if second_state == 2
                expect_value = [0.9*0 + 0.1*token_value(3),  0.9*token_value(1) + 0.1*token_value(2), 0.9*token_value(3) + 0.1*0, 0.9*token_value(2) + 0.1*token_value(1)];
            else
                expect_value = [0.9*token_value(1) + 0.1*0,  0.9*token_value(2) + 0.1*token_value(3), 0.9*0 + 0.1*token_value(1), 0.9*token_value(3) + 0.1*token_value(2)];
            end
            
        case 3 % LH
            %                 opt(3,2) = opt(3,2)+1;
            if second_state == 2 % a1 should be L1. L1 : NOTHING, L2 : BLUE, R2 : RED, R1 : SILVER
                expect_value = [0.5*0 + 0.5*token_value(3),  0.5*token_value(1) + 0.5*token_value(2)];
            else % second_state == 3 L1 : SILVER, L2 : NOTHING, R2 : BLUE, R1 : RED
                expect_value = [0.5*token_value(1) + 0.5*0,  0.5*token_value(2) + 0.5*token_value(3)];
            end
            
        case 4 %HH
            %                 opt(4,2) = opt(4,2)+1;
            if second_state == 2
                expect_value = [0.5*0 + 0.5*token_value(3),  0.5*token_value(1) + 0.5*token_value(2), 0.5*token_value(3) + 0.5*0, 0.5*token_value(2) + 0.5*token_value(1)];
                
                %             I = find(expect_value == max(expect_value));
                %             if ~isempty(find(sbj_data(4) == I))
                %                 cc_HH(sub,2) = cc_HH(sub,2) + 1;
                %                 if isequal(I,[2,4])
                %                     cc_HH(sub,1) = cc_HH(sub,1) + 1;
                %                 end
                %             end
            else
                
                expect_value = [0.5*token_value(1) + 0.5*0,  0.5*token_value(2) + 0.5*token_value(3), 0.5*0 + 0.5*token_value(1), 0.5*token_value(3) + 0.5*token_value(2)];
            end
            
    end
    
    likelihood = expect_value ./ sum(expect_value);
    if (sum(expect_value) == 0)
        likelihood = ones(size(expect_value));
        likelihood = likelihood/sum(likelihood);
    end
        
    if strcmp (condi, 'likelihood')
        outputArg1 = likelihood(action);
        if (ind == 2) || (ind == 4)
            outputArg1 = 2*likelihood(action);
            likelihood = likelihood*2;
        end
        outputArg2 = find(likelihood == max(likelihood));
    elseif strcmp (condi, 'uncertainty')
        expect_expect = [];
        
        switch ind
            case 1
                expect_expect = [0.9*mean([0.9*0 + 0.1*token_value(3),  0.9*token_value(1) + 0.1*token_value(2)]) + 0.1*mean([0.9*token_value(1) + 0.1*0,  0.9*token_value(2) + 0.1*token_value(3)]),...
                    0.1*mean([0.9*0 + 0.1*token_value(3),  0.9*token_value(1) + 0.1*token_value(2)]) + 0.9*mean([0.9*token_value(1) + 0.1*0,  0.9*token_value(2) + 0.1*token_value(3)])];
            case 2
                expect_expect = [0.9*mean([0.9*0 + 0.1*token_value(3),  0.9*token_value(1) + 0.1*token_value(2), 0.9*token_value(3) + 0.1*0, 0.9*token_value(2) + 0.1*token_value(1)]) + 0.1*mean([0.9*token_value(1) + 0.1*0,  0.9*token_value(2) + 0.1*token_value(3), 0.9*0 + 0.1*token_value(1), 0.9*token_value(3) + 0.1*token_value(2)]),...
                    0.1*mean([0.9*0 + 0.1*token_value(3),  0.9*token_value(1) + 0.1*token_value(2), 0.9*token_value(3) + 0.1*0, 0.9*token_value(2) + 0.1*token_value(1)]) + 0.9*mean([0.9*token_value(1) + 0.1*0,  0.9*token_value(2) + 0.1*token_value(3), 0.9*0 + 0.1*token_value(1), 0.9*token_value(3) + 0.1*token_value(2)])];
            case 3
                expect_expect = [0.5*mean([0.5*0 + 0.5*token_value(3),  0.5*token_value(1) + 0.5*token_value(2)]) + 0.5*mean([0.5*token_value(1) + 0.5*0,  0.5*token_value(2) + 0.5*token_value(3)]),...
                    0.5*mean([0.5*0 + 0.5*token_value(3),  0.5*token_value(1) + 0.5*token_value(2)]) + 0.5*mean([0.5*token_value(1) + 0.5*0,  0.5*token_value(2) + 0.5*token_value(3)])];
            case 4
                expect_expect = [0.5*mean([0.5*0 + 0.5*token_value(3),  0.5*token_value(1) + 0.5*token_value(2), 0.5*token_value(3) + 0.5*0, 0.5*token_value(2) + 0.5*token_value(1)]) + 0.5*mean([0.5*token_value(1) + 0.5*0,  0.5*token_value(2) + 0.5*token_value(3), 0.5*0 + 0.5*token_value(1), 0.5*token_value(3) + 0.5*token_value(2)]),...
                    0.5*mean([0.5*0 + 0.5*token_value(3),  0.5*token_value(1) + 0.5*token_value(2), 0.5*token_value(3) + 0.5*0, 0.5*token_value(2) + 0.5*token_value(1)]) + 0.5*mean([0.5*token_value(1) + 0.5*0,  0.5*token_value(2) + 0.5*token_value(3), 0.5*0 + 0.5*token_value(1), 0.5*token_value(3) + 0.5*token_value(2)])];
        end
        likelihood = expect_expect ./ sum(expect_expect);
        outputArg1 = likelihood(action);
        
    else
        outputArg1 = ismember(action, find(likelihood == max(likelihood)));
    end
end

end

