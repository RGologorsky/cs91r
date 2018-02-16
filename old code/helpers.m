% get_payoffs.m helper functions

% converts b, c values into full payoff matrix values
player1_game_vals = @(b, c) [(b-c); (-c); b; 0];
player2_game_vals = @(b, c) [b-c; b; -c; 0];

% computes payoff = dot product of stationary eigenvector and game values
get_payoff = @(evec, game_vals) evec * game_vals;

% evolution.m helper functions

% sigmoid function
sigmoid = @(x) 1/(1.0 + exp(x));

% sigmoid imitation probability
get_imitation_prob = @(pi_r, pi_l) sigmoid(-beta * (pi_r - pi_l));

% returns True = 1 with input probability p
coin_toss = @(p) rand() <= p; 

% generates a new strategy $\in [0,1]^12$ uniformly at random
generate_new_strategy = @() [rand() rand() rand() rand() ...
                             rand() rand() rand() rand() ...
                             rand() rand() rand() rand()];
