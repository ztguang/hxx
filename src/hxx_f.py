#! /usr/bin/python

import numpy as np
import random

class artxx_f ():
    def __init__ (self, raw, lambdas_in):
        self.lambdas = lambdas_in
        self.previous_prediction = None
        self.num_obs_fall = 0
        self.num_obs_rise = 0
        self.num_predictions = 0
        self.num_correct_short_predictions = 0
        self.num_correct_long_predictions = 0
        self.num_uncertain_predictions = 0
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

            alpha_bar = np.multiply (Pi, B[:, O[0]])
            alpha_hat = alpha_bar / alpha_bar.sum ()

            for t in range (0, T - 1):
                # multiply is elementwise, dot is matrix multiplication
                alpha_bar = np.multiply (np.dot (alpha_hat, A), B[:, O[t+1]])
                alpha_hat = alpha_bar / alpha_bar.sum ()

            self.alpha_hat[key] = alpha_hat

    def update (self, raw):
        if self.previous_prediction != None:
            self.num_predictions = self.num_predictions + 1

            if raw < 0.:
                self.num_obs_fall = self.num_obs_fall + 1
                if self.previous_prediction['decision'] == 'short':
                    self.num_correct_short_predictions = self.num_correct_short_predictions + 1
            elif raw > 0.:
                self.num_obs_rise = self.num_obs_rise + 1
                if self.previous_prediction['decision'] == 'long':
                    self.num_correct_long_predictions = self.num_correct_long_predictions + 1
            elif self.previous_prediction['decision'] == 'uncertain':
                self.num_uncertain_predictions = self.num_uncertain_predictions + 1
 
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

            # multiply is elementwise, dot is matrix multiplication
            alpha_bar = np.multiply (np.dot (alpha_hat, A), B[:, O])
            #alpha_bar = np.multiply (alpha_hat, np.dot (A, B[:, O]))

            self.alpha_hat[key] = alpha_bar / alpha_bar.sum ()

            if alpha_bar[0] > alpha_bar[1]:
                short_votes = short_votes + 1
            elif alpha_bar[0] < alpha_bar[1]:
                pass
            else:
                continue # No vote.

            total_votes = total_votes + 1

        long_votes = total_votes - short_votes

        print ("short_votes ({0}), long_votes ({1})".format(short_votes,long_votes))

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
        num_correct_predictions = self.num_correct_short_predictions + self.num_correct_long_predictions
        print ("Predictions: {0} correct out of {1} ({2})".format (num_correct_predictions, self.num_predictions, float(num_correct_predictions)/self.num_predictions))
        print ("Falls and predicted: {0} predicted out of {1} ({2})".format(self.num_correct_short_predictions, self.num_obs_fall, float(self.num_correct_short_predictions)/self.num_obs_fall))
        print ("Rises and predicted: {0} predicted out of {1} ({2})".format(self.num_correct_long_predictions, self.num_obs_rise, float(self.num_correct_long_predictions)/self.num_obs_rise))

if __name__ == '__main__':
    lambdas = {
        'k-0.0134492':{'A':[[0.999999,1e-06,],[0.0207499,0.97925,],],'B':[[1e-06,0.999999,],[0.829887,0.170113,],],'Pi':[1e-06,0.999999,],'split':-0.0134492},
        'k-0.0116564':{'A':[[0.39264,0.60736,],[0.597101,0.402899,],],'B':[[0.181444,0.818556,],[0.190973,0.809027,],],'Pi':[0.936416,0.0635843,],'split':-0.0116564},
        'k0.00110254':{'A':[[0.709008,0.290992,],[0.881111,0.118889,],],'B':[[0.532422,0.467578,],[0.556633,0.443367,],],'Pi':[0.410067,0.589933,],'split':0.00110254},
        'k-0.00484202':{'A':[[0.991692,0.00830823,],[1e-06,0.999999,],],'B':[[0.731118,0.268882,],[1e-06,0.999999,],],'Pi':[0.999999,1e-06,],'split':-0.00484202},
        'k-0.00915241':{'A':[[0.979181,0.0208188,],[1e-06,0.999999,],],'B':[[0.999208,0.00079205,],[0.0402307,0.959769,],],'Pi':[0.999999,1e-06,],'split':-0.00915241},
        'k0.0212183':{'A':[[0.312113,0.687887,],[0.420737,0.579263,],],'B':[[0.897634,0.102366,],[0.984316,0.0156844,],],'Pi':[0.33922,0.66078,],'split':0.0212183},
        'k-0.0044398':{'A':[[0.878221,0.121779,],[0.0699217,0.930078,],],'B':[[0.859974,0.140026,],[0.0559965,0.944003,],],'Pi':[0.999972,2.83995e-05,],'split':-0.0044398},
        'k0.00106572':{'A':[[0.999999,1e-06,],[0.00793998,0.99206,],],'B':[[0.0500032,0.949997,],[0.999999,1e-06,],],'Pi':[1e-06,0.999999,],'split':0.00106572},
        'k0.0134992':{'A':[[0.998747,0.00125298,],[0.00687784,0.993122,],],'B':[[0.629155,0.370845,],[0.993401,0.00659939,],],'Pi':[1e-06,0.999999,],'split':0.0134992},
        'k-0.000516974':{'A':[[0.991667,0.00833275,],[1e-06,0.999999,],],'B':[[0.991597,0.00840277,],[1e-06,0.999999,],],'Pi':[0.999999,1e-06,],'split':-0.000516974},
        'k0.0480553':{'A':[[0.712883,0.287117,],[0.627739,0.372261,],],'B':[[0.996326,0.00367384,],[0.99513,0.00487042,],],'Pi':[0.829617,0.170383,],'split':0.0480553},
        'k0.00954846':{'A':[[0.310494,0.689506,],[0.158019,0.841981,],],'B':[[0.576628,0.423372,],[0.823206,0.176794,],],'Pi':[0.13044,0.86956,],'split':0.00954846},
        'k-0.0216141':{'A':[[0.684116,0.315884,],[0.910812,0.0891883,],],'B':[[0.0373603,0.96264,],[0.11217,0.88783,],],'Pi':[0.513366,0.486634,],'split':-0.0216141},
        'k-0.00708661':{'A':[[0.96353,0.0364699,],[0.0104285,0.989571,],],'B':[[0.797684,0.202316,],[0.0434204,0.95658,],],'Pi':[0.999999,1e-06,],'split':-0.00708661},
        'k0.00155099':{'A':[[0.624376,0.375624,],[0.64916,0.35084,],],'B':[[0.520075,0.479925,],[0.602863,0.397137,],],'Pi':[0.0597071,0.940293,],'split':0.00155099},
        'k-0.00612136':{'A':[[0.3344,0.6656,],[0.232507,0.767493,],],'B':[[0.230627,0.769373,],[0.339964,0.660036,],],'Pi':[0.0949019,0.905098,],'split':-0.00612136},
        'k-0.0092462':{'A':[[0.567925,0.432075,],[0.617362,0.382638,],],'B':[[0.242799,0.757201,],[0.183859,0.816141,],],'Pi':[0.991618,0.00838155,],'split':-0.0092462},
        'k0.00261666':{'A':[[0.992054,0.0079456,],[1e-06,0.999999,],],'B':[[0.999999,1e-06,],[0.12501,0.87499,],],'Pi':[0.999999,1e-06,],'split':0.00261666},
        'k-0.0146226':{'A':[[0.999999,1e-06,],[0.0207052,0.979295,],],'B':[[1e-06,0.999999,],[0.766091,0.233909,],],'Pi':[1e-06,0.999999,],'split':-0.0146226},
        'k0.00459199':{'A':[[0.882353,0.117647,],[0.49575,0.50425,],],'B':[[0.741469,0.258531,],[0.205791,0.794209,],],'Pi':[0.997588,0.00241244,],'split':0.00459199},
                           }

    artxx = artxx_f (np.array([-.1, .05, .01, .6]), lambdas)
    for i in range (0, 100):
        rnd = random.uniform (-1., 1.)
        print rnd
        print artxx.update (rnd)

    artxx.summary ()
