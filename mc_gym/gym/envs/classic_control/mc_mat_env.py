import gym
from gym import spaces
from gym.utils import seeding
import numpy as np
from pymatbridge import Matlab
import random


class mc_mat(gym.Env):

    def __init__(self):
        self.mlab = Matlab(matlab='/Applications/MATLAB_R2014a.app/bin/matlab')
        self.mlab.stop()
        self.mlab.start()

        self.action_space = spaces.Discrete(4)
        self.observation_space = spaces.Discrete(4)
        self.round = 2

        self.LAM = round(random.uniform(0.1, 20),3)
        self.D = round(np.random.uniform(5,40),2)
        self.N_n = np.random.randint(5,20)
        self.N_0 = round(random.uniform(5, 20),3)
        # self.observation = np.array([self.LAM, self.D, self.N_n, self.N_0])#.reshape(1,4)
        self.observation = (self.LAM, self.D, self.N_n, self.N_0)#.reshape(1,4)

        self._reset()

    def _step(self, action):
        assert self.action_space.contains(action)


        res = self.mlab.run('/Users/zonemercy/Documents/MATLAB/mc_gym/Copy_of_oraginal.m', 
                       {'arg1': 1, 'arg2': self.LAM, 'arg3': self.D, 'arg4': self.N_0, 'arg5': self.N_n})

        # esb = res['result']
        esb = [ round(elm, self.round) for elm in res['result'] ]
        result = np.where(esb == np.min(esb))[0]

        if action in result:
            reward = 10
            done = True
        else: 
            reward = 0
            done = False

        # return self.observation, reward, done, {"action": action, "result": result,'esb':esb}
        return self.observation, reward, done, {}


    def _reset(self):
        self.LAM = round(random.uniform(0.1, 20),3)
        self.D = round(np.random.uniform(5,40),2)
        self.N_n = np.random.randint(5,20)
        self.N_0 = round(random.uniform(5, 20),3)
        # self.observation = np.array([self.LAM, self.D, self.N_n, self.N_0])#.reshape(1,4)
        self.observation = (self.LAM, self.D, self.N_n, self.N_0)

        return self.observation

    def _configure(self, nround=3, close_mat=False):
        self.round = nround
        if close_mat == True:
            self.mlab.stop()
