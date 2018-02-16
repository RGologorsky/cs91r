# class Player(object):
#     def __init__(f, (p1_cc, p1_cd, p1_dc, p1_dd),
#                     (p2_cc, p2_cd, p2_dc, p2_dd),
#                     (x_cc, x_cd, x_dc, x_dd))):
#         self.action_strategy = [
#             1: [
#                 (C,C): p1_cc,
#                 (C,D): p1_cd,
#                 (D,C): p1_dc,
#                 (D,D): p1_dd
#                 }
#             2: [
#                 (C,C): p1_cc,
#                 (C,D): p1_cd,
#                 (D,C): p1_dc,
#                 (D,D): p1_dd
#                 }
#         }
#
#         self.transition_strategy = [
#             (C,C): x_cc,
#             (C,D): x_cd,
#             (D,C): x_dc,
#             (D,D): x_dd
#         }
#
#     def prob_next_game_number(prev_round):
#         return self.transition_strategy(prev_round)
#
#     def prob_next_state(prev_round, next_state):
#         (next_action1, _), next_game_number = next_state
#
#         prob_next_C = self.action_strategy.next_game_number.prev_round
#         if next_action1 == C:
#             return prob_next_C
#         return 1 - prob_next_C
#
#
#
# class Game(object):
#     def __init__(f, b, c, p1, p2):
#
#         # reward structure
#         self.g = (b, -b,
#              b - c, 0)
#
#         # transition function between game states
#         self.f = lambda x,y : x * y
#
#         # player1, player2 (with their strategies)
#         self.p1 = p1
#         self.p2 = p2
#
#     def switch_perpective(round):
#         (action1, action2) = round
#         return (action2, action1)
#
#     def prob_transition(prev_state, game_number):
#         x = self.p1.prob_game1(prev_round)
#         y = self.p2.prob_game1(prev_round)
#
#         if game_number == 1:
#             return f(x, y)
#         return 1 - f(X, y)
#
#     # probability of transitioning from previous state -> next state
#     # state = (round, game number), where round = (action1, action2)
#     # action = 1 is cooperate, action = 0 is defect.
#     def prob(prev_state, next_state):
#         prev_state = (prev_action1, prev_action2, prev_game_number)
#         next_state = (next_round, next_game_number)
#
#         next_action1, next_action2 = next_round
#
#         prob_next_game_number = prob_transition(prev_state, next_game_number)
#         prob_next_action1 = self.p1.prob_transition(prev_state, next_action1)
#         prob_next_action2 = self.p2.prob_transition(self.switch_perspective(prev_round),
#                                                         next_action2)
#
#
#         return prob_next_game_number * prob_next_action1 * prob_next_action2
#
# # defines the probability of transitioning to Game State 1 ("better" game state)
# def f(x, y, epsilon = 0):
#   return x * y
#
# # Probability of choosing to switch to Game State 1 given previous round outcome
# def transition_strategy(x_cc, x_cd, x_dc, x_dd):
#     if p1_cc + p1_cd + p1_dc + p1_dd != 1:
#         print("Error: Probabilities for Strategy 1 should sum to 1.")
#     return ["p1_cc": p1_cc,
#             "p1_cd": p1_cd,
#             "p1_dc": p1_dc,
#             "p1_dd": p1_dd}
#
# # Within Game State 1, probability of cooperating given previous round outocme
# def strategy_1(p1_cc, p1_cd, p1_dc, p1_dd, epsilon = 0):
#     if p1_cc + p1_cd + p1_dc + p1_dd != 1:
#         print("Error: Probabilities for Strategy 1 should sum to 1.")
#     return ["p1_cc": p1_cc,
#             "p1_cd": p1_cd,
#             "p1_dc": p1_dc,
#             "p1_dd": p1_dd}
#
# # Within Game State 2, probability of cooperating given previous round outocme
# def strategy_2(p2_cc, p2_cd, p2_dc, p2_dd, epsilon = 0):
#     if p2_cc + p2_cd + p2_dc + p2_dd != 1:
#         print("Error: Probabilities for Strategy 2 should sum to 1.")
#     return ["p2_cc": p2_cc,
#             "p2_cd": p2_cd,
#             "p2_dc": p2_dc,
#             "p2_dd": p2_dd}
#
# def main():
