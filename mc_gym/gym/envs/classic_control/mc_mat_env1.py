import gym
from gym import spaces
from gym.utils import seeding
import numpy as np
import pandas as pd
import random

dataset = pd.read_csv('train5000.csv')

def onehot(observation,v_min,v_max,n=10):
    my_onehot = np.zeros((n))
    value = int(observation//((v_max-v_min)/float(n)))
    if value <= n-1:
        my_onehot[value] = 1.0
    else: my_onehot[int(n-1)] = 1.0
    return my_onehot

class mc_mat1(gym.Env):

    def __init__(self):

        self.action_space = spaces.Discrete(4)
        # self.observation_space = spaces.Discrete(4)
        self.round = 1
        self.esb_log_file = None

        self.NROW = dataset.shape[0]-1
        self.nrow = random.randint(0,self.NROW)

        self.LAM = dataset.loc[self.nrow,'LAM']
        self.D = dataset.loc[self.nrow,'D']
        self.N_0 = dataset.loc[self.nrow,'N_0']
        self.N_n = dataset.loc[self.nrow,'N_n']
        # self.observation = np.array([self.LAM, self.D, self.N_0, self.N_n])#.reshape(1,4)
        self.observation = (self.LAM, self.D, self.N_0, self.N_n)#.reshape(1,4)
        self.observation = tuple(np.concatenate((onehot(self.LAM,0.1,20), onehot(self.D,5,40), 
                            onehot(self.N_0,5,20), onehot(self.N_n,5,20),np.array([self.LAM,self.D,self.N_0,self.N_n])),axis=0))

        self.esb = map(float, dataset.loc[self.nrow,'esb'][1:-1].split(", "))

        self.guess_count = 0
        self.guess_max = 10
        self.log_row = []
        self.log_esb = []
        self.log_esb_tmp = []
        # self._reset()

    def _step(self, action):
        assert self.action_space.contains(action)

        # res = self.mlab.run('/Users/zonemercy/Documents/MATLAB/mc_gym/Copy_of_oraginal.m', 
        #                {'arg1': 1, 'arg2': self.LAM, 'arg3': self.D, 'arg4': self.N_0, 'arg5': self.N_n})

        # esb = res['result']
        esb = [ round(elm, self.round) for elm in self.esb ]
        result = np.where(esb == np.min(esb))[0]

        reward = 0
        done = False
        self.log_esb_tmp.append(esb[action])

        if action in result:
            reward = 10
            done = True
            self.log_row.append(self.nrow)
            self.log_esb.append(self.log_esb_tmp)

        self.guess_count += 1
        if self.guess_count >= self.guess_max:
            done = True
            self.log_row.append(self.nrow)
            self.log_esb.append(self.log_esb_tmp)

        # return self.observation, reward, done, {"action": action, "result": result,'esb':esb}
        return self.observation, reward, done, {}


    def _reset(self):
        self.nrow = random.randint(0,self.NROW)

        self.LAM = dataset.loc[self.nrow,'LAM']
        self.D = dataset.loc[self.nrow,'D']
        self.N_0 = dataset.loc[self.nrow,'N_0']
        self.N_n = dataset.loc[self.nrow,'N_n']
        # self.observation = np.array([self.LAM, self.D, self.N_0, self.N_n])#.reshape(1,4)
        self.observation = (self.LAM, self.D, self.N_0, self.N_n)
        self.observation = tuple(np.concatenate((onehot(self.LAM,0.1,20), onehot(self.D,5,40), 
                            onehot(self.N_0,5,20), onehot(self.N_n,5,20),np.array([self.LAM,self.D,self.N_0,self.N_n])),axis=0))

        self.esb = map(float, dataset.loc[self.nrow,'esb'][1:-1].split(", "))
        self.guess_count = 0
        self.log_esb_tmp = []

        return self.observation

    def _configure(self, nround=3, esb_log=None, save_esb=False, close_mat=False):

        self.round = nround
        self.esb_log_file = esb_log
        if save_esb == True:
            # df = pd.DataFrame()
            df = pd.DataFrame({ 'log_row':self.log_row,'log_esb':self.log_esb},
                  columns=['log_row', 'log_esb'])
            df.to_csv(self.esb_log_file,index=False)
        # if close_mat == True:
        #     self.mlab.stop()
