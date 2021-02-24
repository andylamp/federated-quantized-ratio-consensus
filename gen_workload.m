function [node_capacities] = gen_workload(nodes, opt)
%GEN_NODE_CAPACITY Summary of this function goes here
%   Detailed explanation goes here
%
% Strategies:
%
% - default: odd nodes get 10, even get 20
% - randomised: nodes get a randomised capacity between a range.
%

  % let the user know of the mishap.
  if nodes < 1
    error("we requite at least 1 node")
  end
  
  % if opt was not supplied initialise an empty struct
  if nargin < 2
    opt = struct;
  end  
  
  % check if we want to use randomised values
  if ~isfield(opt, "distribution_type")
    opt.distribution_type = "fixed";
  end
  
  % check if we have a multiplier, and if not use the default.
  if ~isfield(opt, "multiplier")
    opt.multiplier = 1;
  end
  
  % now check if we set fixed values
  if opt.distribution_type == "fixed"
    % check if we have variables for the array
    if ~isfield(opt, "fixed_values")
      opt.fixed_values = [10, 20];
    end
    % generate the node capacity vector
    node_capacities = ones(nodes, 1); 
    % odd nodes are the first bit of the slice.
    node_capacities(1:2:end) = opt.fixed_values(1);
    % even nodes are the second bit of the slice.
    node_capacities(2:2:end) = opt.fixed_values(2);
    
  % our distribution is randomised
  elseif opt.distribution_type == "randomised"
    
    % cap it from 1
    if ~isfield(opt, "min_cap")
      opt.min_cap = 1;
    end
    
    % to hero!
    if ~isfield(opt, "max_cap")
      opt.max_cap = nodes;
    end
    
    % sanity check
    if opt.min_cap > opt.max_cap
      error("max cap cannot be lower than min cap")
    end
    
    % assign a random distribution based on the nodes
    node_capacities = opt.multiplier * randi([opt.min_cap, opt.max_cap], nodes, 1);
  else
    error("unknown option supplied for distribution type");
  end
  
end

