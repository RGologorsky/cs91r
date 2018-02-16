#!/usr/bin/python

import numpy as np
import scipy.sparse.linalg as sla
from ast import literal_eval as make_tuple
import sys

import matlab

def input_to_tuple(str):
    return make_tuple("(" + str + ")")

def main():
  f = lambda a, b: a * b

  # print("Player 1")
  # p1 = get_user_input("p1cc, p1cd, p1dc, p2dd")
  # p2 = get_user_input("p2")
  # x =  get_user_input("x")
 
  # print("Player 2")
  # q1 = get_user_input("q1")
  # q2 = get_user_input("q2")
  # y =  get_user_input("y")
  
  # game_vals1 = raw_input("b1, c1\n")
  # game_vals2 = raw_input("b2, c2\n")

  p1 = (1,0,0,1) 
  p2 = (1,0,0,1) 
  x =  (1,0,0,0)

  q1 = (0,0,0,0)
  q2 = (0,0,0,0)
  y  = (1,0,0,0)

  game_vals1 = (10, 1)
  game_vals2 = (2, 1)

  calculate_payoff(p1, p2, x, q1, q2, y, f, game_vals1, game_vals2)


# Helper functions

# converts b,c values to reward tuple from Player 1's perspective
def p1_game_vals_full(game_vals):
    b, c = game_vals
    return (b - c, -c, b, 0)

# converts b,c values to reward tuple from Player 2's perspective
def p2_game_vals_full(game_vals):
    b, c = game_vals
    return (b - c, b, -c, 0)

# calculates Player 1 and Player 2's payoff, given input parameters
def calculate_payoff(p1, p2, x, q1, q2, y, f, game_vals1, game_vals2):

    evector = get_steady_state(p1, p2, x, q1, q2, y, f)
    p1_game_vals = np.asarray(p1_game_vals_full(game_vals1) +
                              p1_game_vals_full(game_vals2))
    p2_game_vals = np.asarray(p2_game_vals_full(game_vals1) +
                              p2_game_vals_full(game_vals2))
    p1_payoff = np.dot(evector, p1_game_vals)
    p2_payoff = np.dot(evector, p2_game_vals)

    eng = matlab.engine.start_matlab()
    print("Steady state")
    print(evector)
    print("P1 Payoff: %f" % (p1_payoff))
    print("P2 Payoff: %f" % (p2_payoff))

    return (p1_payoff, p2_payoff)

# gets the steady state = eigenvector w/ eigenvalue 1 of Q, given input param.
def get_steady_state(p1, p2, x, q1, q2, y, f):

    # calculate Q = transition matrix
    Q = get_Q(p1, p2, x, q1, q2, y, f)
    # Q = Q.T

    # eng = matlab.engine.start_matlab()
    # [V,D] = eng.eigs(Q,1)

    # print(V)

    evals, evecs = sla.eigs(Q.T, k=1, which='LM', maxiter=5000) # largest eigenvalue/vector
    evec = evecs[:,0].real

    # normalize L1
    evec = evec/np.linalg.norm(evec, ord=1)
    return evec

def convert_to_float(int_tuple):
  return tuple(map(lambda x: float(x), int_tuple))

def get_Q(p1, p2, x, q1, q2, y, f):

    # convert inputs to floats
    p1 = convert_to_float(p1)
    p2 = convert_to_float(p2)
    x  = convert_to_float(x)
    
    q1 = convert_to_float(q1)
    q2 = convert_to_float(q2)
    y  = convert_to_float(y)

    # decompose inputs into probabilities
    (p1cc, p1cd, p1dc, p1dd) = p1
    (p2cc, p2cd, p2dc, p2dd) = p2
    (xcc, xcd, xdc, xdd)     = x
    (q1cc, q1cd, q1dc, q1dd) = q1
    (q2cc, q2cd, q2dc, q2dd) = q2
    (ycc, ycd, ydc, ydd)     = y

    return np.matrix([
     [
      f(xcc, ycc)*p1cc*q1cc,
      f(xcc, ycc)*p1cc*(1 - q1cc),
      f(xcc, ycc)*(1 - p1cc)*q1cc,
      f(xcc, ycc)*(1 - p1cc)*(1 - q1cc),

      (1 - f(xcc, ycc))*p2cc* q2cc,
      (1 - f(xcc, ycc))*p2cc * (1 - q2cc),
      (1 - f(xcc, ycc)) * (1 - p2cc)*q2cc,
      (1 - f(xcc, ycc))* (1 - p2cc) * (1 - q2cc)
     ],

     [
      f(xcd, ydc)*p1cd*q1dc,
      f(xcd, ydc)*p1cd*(1 - q1dc),
      f(xcd, ydc)*(1 - p1cd)*q1dc,
      f(xcd, ydc)*(1 - p1cd)*(1 - q1dc),

      (1 - f(xcd, ydc))*p2cd*q2dc,
      (1 - f(xcd, ydc))*p2cd*(1 - q2dc),
      (1 - f(xcd, ydc))*(1 - p2cd)*q2dc,
      (1 - f(xcd, ydc))*(1 - p2cd)*(1 - q2dc)
     ],

     [
      f(xdc, ycd)*p1dc*q1cd,
      f(xdc, ycd)*p1dc*(1 - q1cd),
      f(xdc, ycd)* (1 - p1dc)*q1cd,
      f(xdc, ycd)*(1 - p1dc)*(1 - q1cd),

      (1 - f(xdc, ycd))*p2dc*q2cd,
      (1 - f(xdc, ycd))*p2dc*(1 - q2cd),
      (1 - f(xdc, ycd))*(1 - p2dc)* q2cd,
      (1 - f(xdc, ycd))*(1 - p2dc)*(1 - q2cd)
     ],

     [
      f(xdd, ydd)*p1dd*q1dd,
      f(xdd, ydd)*p1dd*(1 - q1dd),
      f(xdd, ydd)*(1 - p1dd)*q1dd,
      f(xdd, ydd)*(1 - p1dd)*(1 - q1dd),

      (1 - f(xdd, ydd))*p2dd*q2dd,
      (1 - f(xdd, ydd))*p2dd*(1 - q2dd),
      (1 - f(xdd, ydd))*(1 - p2dd)*q2dd,
      (1 - f(xdd, ydd))*(1 - p2dd)*(1 - q2dd)
     ],

     [
      f(xcc, ycc)*p1cc*q1cc,
      f(xcc, ycc)*p1cc*(1 - q1cc),
      f(xcc, ycc)*(1 - p1cc)*q1cc,
      f(xcc, ycc)*(1 - p1cc)*(1 - q1cc),

      (1 - f(xcc, ycc))*p2cc* q2cc,
      (1 - f(xcc, ycc))*p2cc * (1 - q2cc),
      (1 - f(xcc, ycc)) * (1 - p2cc)*q2cc,
      (1 - f(xcc, ycc))* (1 - p2cc) * (1 - q2cc)
     ],

     [
      f(xcd, ydc)*p1cd*q1dc,
      f(xcd, ydc)*p1cd*(1 - q1dc),
      f(xcd, ydc)*(1 - p1cd)*q1dc,
      f(xcd, ydc)*(1 - p1cd)*(1 - q1dc),

      (1 - f(xcd, ydc))*p2cd*q2dc,
      (1 - f(xcd, ydc))*p2cd*(1 - q2dc),
      (1 - f(xcd, ydc))*(1 - p2cd)*q2dc,
      (1 - f(xcd, ydc))*(1 - p2cd)*(1 - q2dc)
     ],

     [
      f(xdc, ycd)*p1dc*q1cd,
      f(xdc, ycd)*p1dc*(1 - q1cd),
      f(xdc, ycd)* (1 - p1dc)*q1cd,
      f(xdc, ycd)*(1 - p1dc)*(1 - q1cd),

      (1 - f(xdc, ycd))*p2dc*q2cd,
      (1 - f(xdc, ycd))*p2dc*(1 - q2cd),
      (1 - f(xdc, ycd))*(1 - p2dc)* q2cd,
      (1 - f(xdc, ycd))*(1 - p2dc)*(1 - q2cd)
     ],

     [
      f(xdd, ydd)*p1dd*q1dd,
      f(xdd, ydd)*p1dd*(1 - q1dd),
      f(xdd, ydd)*(1 - p1dd)*q1dd,
      f(xdd, ydd)*(1 - p1dd)*(1 - q1dd),

      (1 - f(xdd, ydd))*p2dd*q2dd,
      (1 - f(xdd, ydd))*p2dd*(1 - q2dd),
      (1 - f(xdd, ydd))*(1 - p2dd)*q2dd,
      (1 - f(xdd, ydd))*(1 - p2dd)*(1 - q2dd)
     ],
    ])

def get_user_input(string):
    str_q = "Please enter " + string + " as a comma seperated sequence of values:"
    q = raw_input(str_q)
    print("\n")
    return input_to_tuple(q)


if __name__ == '__main__':
    main()
