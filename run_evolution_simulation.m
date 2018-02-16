function [] ...
        = run_evolution_simulation()
    %% Constants

    % poplution size
    N = 100;

    % number of time steps
    T = 10^5;
    time_vec = 1:T;

    % probability of mutation
    mu = 0.01;
    
    % selection pressure
    beta = 1;

    % game1 & game2 values. INTRODUCE SCALING PARAMETER?
    game2 = [1.2 1];
    game1 = [1.9 1];
    
    % error rate
    eps = 0.01;

    % E(new strategies) = mu*T
    max_num_strategies = 2*mu*T; 

  
    
    %% Setup variables

    % debugging
    num_times_sampled_same = 0;
    
    % list of all strategies that have evolved
    strategies_list = zeros(max_num_strategies,12);

    % keeps track of the total number of strategies evolved.
    total_num_strategies = 1; 
    
    % keeps track of the number of strategies in the current populaion
    curr_num_strategies = 1;
    
    % initial strategy = pure ALLD, w/ error rate eps
    all_d_eps = [eps eps eps eps ...
                 eps eps eps eps ...
                 eps eps eps eps];
           
    strategies_list(1,:) = all_d_eps;

    % stores the composition of the current population: specifically,
    % stores strategy index in stategies_list, strategy frequency
    curr_pop = [[1, N]]; % N indivuals start w/ ALLD.

    % payoffs, coop rate, & Game1 freq of current_pop strategy encounters
    curr_payoffs = zeros(curr_num_strategies, curr_num_strategies);
    curr_coops   = zeros(curr_num_strategies, curr_num_strategies);
    curr_game1   = zeros(curr_num_strategies, curr_num_strategies);

    % below time data vectors store strategy population data at each time step
    avg_coop_data   = time_vec; % avg. coop rate stored per time step
    avg_payoff_data = time_vec; % avg. payoff stored per time step
    avg_game1_data  = time_vec; % avg. freq players in Game1 per time step

    % how often each strategy was played over course of full simulation
    cum_strategy_counts = zeros(max_num_strategies,1);
    
    %% Helper Functions

    % sigmoid function
    sigmoid = @(x) 1/(1.0 + exp(x));

    % sigmoid imitation probability
    get_imitation_prob = @(pi_r, pi_l) sigmoid(-beta * (pi_r - pi_l));

    % returns True = 1 with input probability p
    coin_toss = @(p) rand() <= p; 

    % generates a memory-1 strategy $\in [0,1]^12$ uniformly at random
    generate_mem1_strategy = @() [rand() rand() rand() rand() ...
                                 rand() rand() rand() rand() ...
                                 rand() rand() rand() rand()];
                             
    % generates a pure memory-1 strategy, with error rate epsilon.
    function strategy = generate_pure_mem1_strategy()
        coop_defect = [coin_toss(0.5) coin_toss(0.5), coin_toss(0.5), ...
                       coin_toss(0.5) coin_toss(0.5), coin_toss(0.5), ...
                       coin_toss(0.5), coin_toss(0.5), coin_toss(0.5), ...
                       coin_toss(0.5), coin_toss(0.5), coin_toss(0.5)];
        % add in error: change 0 -> epsilon, and 1 -> 1 - epsilon
        strategy = (1 - 2*eps)*coop_defect + eps;
    end

    % strategy's avg payoff against N - 1 other individuals in current pop
    % does not count individual w/ strategy playing against himself 
    function avg_payoff = get_avg_strategy_payoff(strategy_index)
        all_payoff = curr_payoffs(strategy_index,:) * curr_pop(:,2);
        avg_payoff = (all_payoff - curr_payoffs(strategy_index,strategy_index))/(N-1);
    end
    
    % updates strategy list, current pop, payoffs, #strategies, etc.
    function add_strategy(new_strategy)
            [already_listed,strategy_index] = ...
                ismember(new_strategy,strategies_list(1:total_num_strategies,:),'rows');
            
            already_in_pop = false;
            
            if already_listed % check if curr pop has the random strategy 
                [already_in_pop, pop_index] = ...
                    ismember(strategy_index,curr_pop(1:curr_num_strategies, 1));
                
                if already_in_pop
                    curr_pop(pop_index,2) = curr_pop(pop_index,2) + 1;
                else
                    % update current pop #strategies(also used as index)
                    curr_num_strategies = curr_num_strategies + 1;
                    curr_pop(curr_num_strategies,:) = [strategy_index 1];
                end
                
            else % strategy not listed, so not in current pop either.
            
                % total num strategies also used as an index 
                total_num_strategies = total_num_strategies + 1;
                curr_num_strategies  = curr_num_strategies  + 1;
                                
                strategies_list(total_num_strategies,:)   = new_strategy;
                cum_strategy_counts(total_num_strategies) = 0;
                
                curr_pop(curr_num_strategies,:) = [total_num_strategies 1];
            end
            
            if ~already_in_pop
                % calculate payoffs for newly added strategy
                expand_curr_payoffs()
            end
    end
    
    function expand_curr_payoffs()
        new_strategy_index = curr_num_strategies;
        new_strat = strategies_list(curr_pop(new_strategy_index,1),:);
        
        for j = 1:curr_num_strategies
            
            strategy_j = strategies_list(curr_pop(j, 1),:);
            
            [payoffs, frac_coops, frac_game1] = ...
                get_stats(game1, game2, new_strat, strategy_j);

            payoff_i    = payoffs(1);
            payoff_j    = payoffs(2); 
            
            frac_coop_i = frac_coops(1);
            frac_coop_j = frac_coops(2);

            curr_payoffs(new_strategy_index,j) = payoff_i;
            curr_coops(new_strategy_index,j)   = frac_coop_i;
            curr_game1(new_strategy_index,j)   = frac_game1;

            curr_payoffs(j,new_strategy_index) = payoff_j;
            curr_coops(j,new_strategy_index)   = frac_coop_j;
            curr_game1(j,new_strategy_index)   = frac_game1;
        end
    end
                
    
    % deletes a strategy (updates current pop, #strategies, payoffs, etc)
    function delete_strategy(strategy_index)
        curr_pop(strategy_index,:) = []; % delete associated row
        shrink_curr_payoffs(strategy_index)
        % curr_num_strategies also serves as index
        curr_num_strategies = curr_num_strategies - 1; 
    end
    
    function shrink_curr_payoffs(index)
        curr_payoffs(index,:) = []; % delete associated row
        curr_payoffs(:,index) = []; % delete associated column
        
        curr_coops(index,:) = []; % delete associated row
        curr_coops(:,index) = []; % delete associated column
        
        curr_game1(index,:) = []; % delete associated row
        curr_game1(:,index) = []; % delete associated column 
    end
    
    % sample strategy weighted by strategy prevalence
    sample_strategy= @() ...
        randsample(curr_num_strategies, 1, true, curr_pop(:,2)/N);
    
    function record_timestep_data(timestep)
        % update strategy cumulative totals
        
        indices_in_play = curr_pop(:,1);
        % assignin('base', 'curr_pop', curr_pop)
        % assignin('base', 'cum_strategy_counts', cum_strategy_counts)
        
        cum_strategy_counts(indices_in_play) = ...
            cum_strategy_counts(indices_in_play) + curr_pop(:,2);
        
        % each game contributes 2 values to the overall pool of results; 
        % (N choose 2) games.
        num_contribs = N * (N - 1);
                
        % for each strategy, avg the sum of the contributions of N-1 games
        % (discard 1 game of playing yourself)
        
        % total strategy payoff/coop/game1 in matchup w/other N-1 opponents
        strategy_payoffs = curr_payoffs * curr_pop(:,2) - diag(curr_payoffs); 
        strategy_coops   = curr_coops   * curr_pop(:,2) - diag(curr_coops);
        strategy_game1s  = curr_game1   * curr_pop(:,2) - diag(curr_game1);
        
        % weight by each strategy's freq; sum for overall payoff/coop/game1
        avg_payoff = sum(strategy_payoffs .* curr_pop(:,2)/num_contribs);
        avg_coop   = sum(strategy_coops   .* curr_pop(:,2)/num_contribs);
        avg_game1  = sum(strategy_game1s  .* curr_pop(:,2)/num_contribs);
       
        avg_payoff_data(timestep) = avg_payoff;
        avg_coop_data(timestep)   = avg_coop;
        avg_game1_data(timestep)  = avg_game1;
    end
    
    function plot_timestep_data()
        
        figure(1);
        subplot(2,2,1);
        plot(time_vec, avg_coop_data); 
        title('Evolution of Avg. Cooperation Frequency');
        xlabel('Timestep');
        ylabel('Avg. Fraction of Cooperation');

        subplot(2,2,2); 
        
        % plot payoff for reward in game 2, R = b - c.
        game2_reward_payoff = game2(1) - game2(2);
        % plot(time_vec,ones(size(T))*game2_reward_payoff,'--rc);
        % plot([1 T], [game2_reward_payoff game2_reward_payoff],'--r');
    
        plot(time_vec, avg_payoff_data);
        title('Evolution of Overall Avg. Payoff');
        xlabel('Timestep');
        ylabel('Avg. Payoff');
        
        subplot(2,2,3); 
        plot(time_vec, avg_game1_data);  
        title('Evolution of Avg. Game1 Frequency');
        xlabel('Timestep');
        ylabel('Avg. Fraction of Time in Game1');
    end

    %% Main Evolution Loop
    for timestep = time_vec

        % INDIVIDUALS UPDATE STRATEGIES

        % w/probability mu, someone invents a new strategy;else can "learn"
        % (imitate) strategy from randomly chosen role model, w/probability
        % determined by selection pressure & relative strategy performance

        if coin_toss(mu) 
            % sample strategy weighted by strategy prevalence
            old_strategy = randsample(curr_num_strategies, 1, true, curr_pop(:,2)/N);

            new_strategy = generate_pure_mem1_strategy();
            
            % old_strategy lost adherent
            curr_pop(old_strategy,2) = curr_pop(old_strategy,2) - 1; 
            
            if curr_pop(old_strategy,2) == 0
                % old_strategy went extinct in current population
                delete_strategy(old_strategy)
            end
            
            % updates current pop,strategy list,etc - if not already exists
            add_strategy(new_strategy);

        else % randomly sample "Learner" & "Rolemodel" (w/replacement)

            learner   = randsample(curr_num_strategies, 1, true, curr_pop(:,2)/N);
            rolemodel = randsample(curr_num_strategies, 1, true, curr_pop(:,2)/N);
            
            if learner == rolemodel
                num_times_sampled_same = num_times_sampled_same + 1;
            end
            
            

            % avg payoff agaisnt N - 1 other individuals in the population
            pi_learner   = get_avg_strategy_payoff(learner);
            pi_rolemodel = get_avg_strategy_payoff(rolemodel);

            imitation_prob = get_imitation_prob(pi_rolemodel, pi_learner);
 
            
            if coin_toss(imitation_prob) 
                % learner switches to rolemodel strategy
                curr_pop(learner, 2)   = curr_pop(learner, 2)   - 1;
                curr_pop(rolemodel, 2) = curr_pop(rolemodel, 2) + 1;
                
                if curr_pop(learner, 2) == 0
                    % learner strategy went extinct
                    delete_strategy(learner);
                end
            end
        end
        
        record_timestep_data(timestep);
    end
    
    num_times_sampled_same
    
    % PLOT GRAPHS
    plot_timestep_data();
    
    %% Assign variables to base for later examination
    
    % strategy cumulative totals
    cum_strategy_counts = cum_strategy_counts(1:total_num_strategies)/sum(cum_strategy_counts(1:total_num_strategies));
    
    assignin('base', 'curr_pop', curr_pop);
    assignin('base', 'strategies_list', strategies_list);
    assignin('base', 'cum_strategy_counts', cum_strategy_counts);
    
    %[most_freq_strategy_counts, indices] = ...
    %    sort(cum_strategy_counts, 'descend');
    
    % most_freq_strategies = strategies_list(indices,:);
    % assignin('base', 'top_10_strategy_counts', most_freq_strategy_counts(1:10));
    % assignin('base', 'top_10_strategies_by_count', most_freq_strategies(1:10,:));
     
end
