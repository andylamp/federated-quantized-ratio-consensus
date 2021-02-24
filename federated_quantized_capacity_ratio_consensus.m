%% Federated quantized capacity ratio consensus /w finite time termination
%
% This is the main stub that performs the experiments in our paper. We test
% across a variety of parameters; concretely we define a parameter tuple as
% the node count, delay and run it across a number of trials - in the paper
% we used 10 trials.
%
% For more details, please see either the README.md or our paper which can
% be found here: https://arxiv.org/abs/xxxx.xxxx
%
% Authors:
%
%  - Andreas A. Grammenos (ag926@cl.cam.ac.uk)
%  - Apostolos Rikos (rikos@kth.se)
%  - Themistoklis Charalambous (themistoklis.charalambous@aalto.fi)
%
% License: GPLv3
%

%% Initialisation
%
clear; clc; close all;

% for reproducibility, fix the seed
rng("default")

% the maximum iterations to go through
max_iter = 100;
% the graph connectivity target
graph_connectivity = 0.5;
% the number of nodes
nodes_to_test = [20, 50]; % , 100, 150, 200, 400, 600
% nodes for "large scale" testing
% nodes_to_test = [20, 200, 500, 1000, 5000, 10000];
node_to_test_len = length(nodes_to_test);
% set the node length as a convenience variable
node_len_array = 1:node_to_test_len;

% trials for regular testing
trials = 2;
% trials for large scale testing
% trials = 5;
trials_arr = 1:trials;

% the quantisation step to use
quantization_step = 100;

for t=trials_arr
    fprintf("\n** Running trial %d\n", t);
    
    for n=node_to_test_len
      nodes = nodes_to_test(n);
      
      % generate the capacity for z0
      z0 = gen_workload(nodes);
      % now generate the capacity for y0
      y0_opt.multiplier = quantization_step;
      y0_opt.distribution_type = "randomised";
      y0 = gen_workload(nodes, y0_opt);
      % firstly, lets generate the graph
      [~, diameter, nodes, adjMatrix] = gen_graph(nodes);
    end
    fprintf("** Finished trial %d\n", t);
end

