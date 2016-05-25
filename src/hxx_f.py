#! /usr/bin/python

import numpy as np
import random

class artxx_f ():
    def __init__ (self, raw, lambdas_in):
        self.lambdas = lambdas_in
        self.previous_prediction = None
        self.num_predictions = 0
        self.num_correct_predictions = 0
        self.alpha_hat = {}

        lambdas = self.lambdas

        for key in lambdas:
            split = lambdas[key]['split']
            A = np.array (lambdas[key]['A'])
            B = np.array (lambdas[key]['B'])
            Pi = np.array (lambdas[key]['Pi'])

            N = Pi.shape[0]
            M = B.shape[1]
            T = len (raw)

            O = np.zeros (T, dtype=np.int8)
            np.place (O, raw >= split, 1)

            alpha_bar = np.dot (Pi, B[:, O[0]])
            alpha_hat = alpha_bar / alpha_bar.sum ()

            for t in range (0, T - 1):
                alpha_bar = np.multiply (alpha_hat, np.dot (A, B[:, O[t+1]]))

                alpha_hat = alpha_bar / alpha_bar.sum ()

            self.alpha_hat[key] = alpha_hat

    def update (self, raw):
        if self.previous_prediction != None:
            self.num_predictions = self.num_predictions + 1

            if raw < 0. and self.previous_prediction['decision'] == 'short':
                self.num_correct_predictions = self.num_correct_predictions + 1
            elif raw > 0. and self.previous_prediction['decision'] == 'long':
                self.num_correct_predictions = self.num_correct_predictions + 1
 
        lambdas = self.lambdas

        short_votes = 0
        total_votes = 0

        for key in lambdas:
            split = lambdas[key]['split']
            A = np.array (lambdas[key]['A'])
            B = np.array (lambdas[key]['B'])
            Pi = np.array (lambdas[key]['Pi'])

            N = Pi.shape[0]
            M = B.shape[1]

            alpha_hat = self.alpha_hat[key]
            
            if raw < split:
                O = 0
            else:
                O = 1

            alpha_bar = np.multiply (alpha_hat, np.dot (A, B[:, O]))

            self.alpha_hat[key] = alpha_bar / alpha_bar.sum ()

            print "alpha_bar", alpha_bar

            if alpha_bar[0] > alpha_bar[1]:
                short_votes = short_votes + 1
            elif alpha_bar[0] < alpha_bar[1]:
                pass
            else:
                continue # No vote.

            total_votes = total_votes + 1

        long_votes = total_votes - short_votes

        if short_votes > long_votes:
            new_prediction = { 'decision': 'short',
                               'p': float(short_votes) / total_votes }
        elif short_votes < long_votes:
            new_prediction = { 'decision': 'long',
                               'p': float(long_votes) / total_votes }
        else:
            new_prediction = { 'decision': 'uncertain' }

        self.previous_prediction = new_prediction

        return new_prediction

    def summary (self):
        print ("Predictions: {0} correct out of {1}".format (self.num_correct_predictions, self.num_predictions))

if __name__ == '__main__':
    lambdas = { 'first': { 'split': 0.,
                           'A': np.array ([[.6, .4], [.5, .5]]),
                           'B': np.array ([[.74, .26], [.16, .84]]),
                           'Pi': np.array ([0.3, 0.7])} }

    artxx = artxx_f (np.array([-.1, .05, .01, .6]), lambdas)
    for i in range (0, 100):
        rnd = random.uniform (-1., 1.)
        print rnd
        print artxx.update (rnd)

    artxx.summary ()
