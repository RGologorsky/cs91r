% precompute avg. payoff, cooperation freq., time in game1
function [] = precompute_stats()
    
    %% Constants
    
    % game1 & game2 values. 
    game2 = [1.2 1];
    game1 = [1.8 1];
    
    % error rate
    eps = 0.01;
    
    num_strategies = 2^12;
    payoffs      = zeros(num_strategies, num_strategies);
    frac_coops   = zeros(num_strategies, num_strategies);
    frac_game1s  = zeros(num_strategies, num_strategies);
    
    %% Helper Functions - convert between strategy/binary representation
    function strategy = num_to_strat(num)
        coop_defect = de2bi(num);
        if size(coop_defect,2) < 12
            coop_defect(12) = 0; % extend strategy to be a 12 component vector
        end
        
        % add in error: change 0 -> epsilon, and 1 -> 1 - epsilon
        strategy = (1 - 2*eps)*coop_defect + eps;
    end
    
    function num = strat_to_num(strategy)
        % undo error: change epsilon -> 0, and 1 - epsilon -> 1
        binary = (strategy - eps)/(1-2*eps);
        num = bi2de(binary);
    end
    
    %% Main Loop
    for i = 1:num_strategies
        for j = i:num_strategies
            strat_i = num_to_strat(i);
            strat_j = num_to_strat(j);
            
            [payoffs, frac_coops, frac_game1] = ...
                get_stats(game1, game2, strat_i, strat_j);

            payoff_i    = payoffs(1);
            payoff_j    = payoffs(2); 
            
            frac_coop_i = frac_coops(1);
            frac_coop_j = frac_coops(2);

            payoffs(i,j)    = payoff_i;
            frac_coops(i,j) = frac_coop_i;
            frac_game1s(i,j)= frac_game1;

            payoffs(j,i)    = payoff_j;
            frac_coops(j,i) = frac_coop_j;
            frac_game1s(j,i)= frac_game1; % save 2x space w/i > j only
        end
    end
    
    assignin('base', 'precomputed_payoffs', payoffs);
    assignin('base', 'precomputed_frac_coops', frac_coops);
    assignin('base', 'precomputed_frac_game1s', frac_game1s);
end