import sys
import os

import numpy as np
from scipy.sparse import csr_matrix
from scipy.special import logsumexp
import pandas as pd
from sklearn import mixture
from sklearn.naive_bayes import GaussianNB
import util
import time

class MarkovRandomField:
    def __init__(self, n, edges, obs, args, n_state=10, edge_potential=None, label=None, label_prob=None, n_iter=10, tol=1e-2, covariance_type='full'):
        self.n = n
        self.obs = obs
        self.edges = edges
        self.edges_matrix = csr_matrix((np.ones(np.size(edges,0)), (edges[:, 0], edges[:, 1])), shape=(n, n))
        self.n_edges = np.size(edges,0)
        self.reverse_index = np.zeros(self.n_edges)
        for i in range(int(self.n_edges/2)):
            self.reverse_index[i]=i+int(self.n_edges/2)
            self.reverse_index[i+int(self.n_edges/2)]=i
        
        self.n_state = args.n
        self.edge_potential = np.zeros((n_state, n_state))
        self.label = np.zeros(n_state)
        self.label_prob = np.zeros(n)
        self.n_iter = n_iter
        self.tol = tol
        self.cov = covariance_type
        self.last = []
        self.args = args

    def get_state(self, i):
        return self.label[i]

    def get_neighbor(self, i):
        return self.edges[i].nonzero()

    def check_edge(self, i, j):
        return (self.edges_matrix[i,j]>0)

    def init_gmm(self):
        gmm_model = mixture.GaussianMixture(n_components=self.n_state, covariance_type=self.cov)
        gmm_model.fit(self.obs)
        label = gmm_model.predict(self.obs)
        label_prob = gmm_model.predict_proba(self.obs)
        self.label = label
        
        self.label_prob = np.zeros(shape=label_prob.shape)
        for i in range(self.n_state):
            self.label_prob[:,i] = util.get_multi_normal_pdf(gmm_model.means_[i,:],gmm_model.covariances_[i,:,:], self.obs)
        
        #self.label_prob = util.log_transform(label_prob)
        self.gaussian = self.label_prob.copy()

    def init_trans(self):
        # trans_matrix = np.ones((self.n_state, self.n_state))
        #edges = self.edges.toarray()
        #edges = self.edges
        #edges[edges > 0] = 1

        #label_one_hot = np.zeros((self.n_state, self.n))
        #label_one_hot[self.label, np.arange(self.n)] = 1
        
        label_one_hot = csr_matrix((np.ones(self.n), (self.label, np.arange(self.n))), shape=(self.n_state, self.n))
        
        label_one_hot_t = np.transpose(label_one_hot)
        
        trans_matrix = label_one_hot.dot(self.edges_matrix).dot(label_one_hot_t)
        
        #self.edge_potential = util.log_transform(trans_matrix.toarray())

        trans_matrix = util.log_transform(trans_matrix.toarray())
        trans_matrix = trans_matrix - logsumexp(trans_matrix)
        self.edge_potential = trans_matrix
        
        
    def predict(self):
        self.label = np.array(np.argmax(np.transpose(self.label_prob), axis=0))[0]
    
    def check_converge(self):
        if (len(self.last) >= 2 and abs(self.last[1] - self.last[0]) < self.tol):
            return True
        else:
            return False

    def loop_belief_propagation(self, max_iter=100):
        #check
    
        #reverse index/neighbor
        reverse_matrix =  csr_matrix((np.ones(self.n_edges), (np.arange(self.n_edges),self.reverse_index)), shape=(self.n_edges, self.n_edges))
        #message
        #message = np.ones((self.n_state,self.n_edges))
        message = csr_matrix(util.log_transform(np.ones((self.n_state,self.n_edges))))

        #receive t
        t = csr_matrix((np.ones(self.n_edges), (np.arange(self.n_edges),self.edges[:,1])), shape=(self.n_edges, self.n))
        #send f
        f = csr_matrix((np.ones(self.n_edges), (np.arange(self.n_edges),self.edges[:,0])), shape=(self.n_edges, self.n))
        #print("lbp.")
        util.print_log(time.ctime() + " loopy belief propagation.", self.args.o + "/log.txt")
        self.last = []
        for i in range(max_iter):

            belief_new = csr_matrix((np.transpose(self.gaussian))) + message.dot(t)
            #normalization
            belief = csr_matrix(belief_new - util.sparse_logsumexp_row(belief_new))
            #print(logsumexp(belief_new))
            
            source_b = belief.dot(np.transpose(f))
            source_message = message.dot(reverse_matrix)
            pairwise = np.repeat(self.edge_potential, self.n_edges).reshape(self.n_state,self.n_state,self.n_edges)

            message = pairwise + source_b.toarray() - source_message.toarray()
            message = csr_matrix(logsumexp(message,axis=1))

            #message = np.ones(belief_new.shape).dot(logsumexp(np.transpose(self.label_prob)) + belief.dot(np.transpose(f)) - message.dot(reverse_matrix))
            self.last.append(np.sum(belief))
            if len(self.last)>2:
                tmp = self.last.pop(0)
            if self.check_converge():
                break
                #pass
        util.print_log(time.ctime() + " Number of iteration: " + str(i), self.args.o + "/log.txt")

        self.label_prob = np.transpose(belief)
        self.label = np.array(np.argmax(belief, axis=0))[0]

        #print(np.transpose(np.argmax(self.label_prob, axis=1)))

    def solve(self, max_iter=10):
        last = []
        for i in range(max_iter):
            self.loop_belief_propagation()
            util.print_log(time.ctime() + " EM iteration " + str(i), self.args.o + "/log.txt")
            #em algorithm
            
            #print(self.gaussian)
            
            label_prev = csr_matrix((np.ones(self.n), (np.arange(self.n),self.label)), shape=(self.n, self.n_state))
            neighbor_likehood = self.edges_matrix.dot(label_prev).dot(self.edge_potential)
            posterior = (self.gaussian +  neighbor_likehood)
            posterior = posterior - logsumexp(posterior, axis=0)
            #posterior = logsumexp(self.gaussian +  neighbor_likehood, axis=1)           
            estimate_label = np.array(np.argmax(posterior, axis=1))
            target = csr_matrix(np.transpose(self.gaussian)).dot(self.edges_matrix).dot(label_prev)
            
            #clf = GaussianNB()
            #clf.fit(self.obs, estimate_label)
            #self.gaussian = clf.predict_log_proba(self.obs)
            self.gaussian = np.zeros(shape=self.gaussian.shape)
            for i in range(self.n_state):
                #self.gaussian[:,i] = util.get_multi_normal_pdf(gmm_model.means_[i,:],gmm_model.covariances_[i,:,:], self.obs)
                data = self.obs[estimate_label==i]

                if len(data)>0:
                    self.gaussian[:,i] = util.get_multi_normal_pdf(np.mean(data, axis=0),np.cov(data, rowvar=0), self.obs)
            
            label_one_hot = csr_matrix((np.ones(self.n), (estimate_label, np.arange(self.n))), shape=(self.n_state, self.n))
            label_one_hot_t = np.transpose(label_one_hot)
            trans_matrix = label_one_hot.dot(self.edges_matrix).dot(label_one_hot_t)
            trans_matrix = util.log_transform(trans_matrix.toarray())
            trans_matrix = trans_matrix - logsumexp(trans_matrix)
            self.edge_potential = trans_matrix
            
            if np.array_equal(last, estimate_label):
                break
            last = estimate_label

    def __del__(self):
        pass
