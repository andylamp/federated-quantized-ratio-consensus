close all
clear;
clc

rng('default');

%largeconstant=10^10                %%% WORK ALSO FOR RING
connectivity=0.5;                   
nodes=200;
iter=100;

Ad = zeros(nodes,nodes);

for i=1:nodes
    for j=1:nodes
        if rand()>connectivity
            Ad(i,j)=1;              % AD - Adjacency Martix
        end
    end
end

for i=1:nodes
    Ad(i,i)=1;
end


diam=1;
AAA = Ad;
while (find(AAA==0)>0)
    AAA = AAA*Ad;
    diam=diam+1;
end


quant_step=1000;   % we can multiply by max of z0
y0=quant_step*randi(100,nodes,1);
z0=2*ones(nodes,1);
for j=1:nodes
    if (mod(j,2)==1)
        z0(j) = 100;
    else
        z0(j) = 300;
    end
end
y=y0
z=z0




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   - RING
% Ad = zeros(nodes,nodes);
% 
% for i=1:nodes-1
%     Ad(i+1,i)=1;
%     Ad(i,i)=1;
% end
% Ad(1,nodes)=1;
% 
% y0 = [4; 8; 10; 3; 5; 1];
% z0 = [1; 1; 1; 1; 1; 1];
% y=y0;
% z=z0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Ad = [0 1 1 0 0 0 0; 1 0 0 1 0 0 0; 0 0 0 0 0 1 0; 0 0 0 0 0 0 1; 1 1 1 1 0 0 0; 0 0 0 0 1 0 1; 0 0 0 0 1 0 0];   % CDC paper plot
%  nodes = 7;
%  y0 = [15; 5; 11; 4; 3; 13; 9]
%  z0 = [1; 1; 1; 1; 1; 1; 1]
%  y=y0;
%  z=z0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=1:nodes
%     Ad(i,i)=1;
% end



transm_z = zeros(nodes,nodes);
transm_y = zeros(nodes,nodes);


initial_aver = sum(y0)/sum(z0);

first = [];
first2 = [];

nodes_states_y = y0;
nodes_states_z = z0;



y
z
initial_aver

first_max = [];
first_min = [];


max_plot = [];
min_plot = [];

nonz = zeros(1,nodes);


   

% for k=1:iter
%     
%     for j=1:nodes
%         if (z(j) >0)
%             nodes_states_y(j) = ceil(y(j)/z(j));
%         end
%     end
%     
%     first2 = [first2; (nodes_states_y)'];
%     nodes_states_y = nodes_states_y.*z0;
%     first = [first; (nodes_states_y)'];
%     
%    
%     
%     if (mod(k,diam)==1)
%         states_for_vote = ((y)./(z));
%         max_vot = ceil(states_for_vote);
%         min_vot = floor(states_for_vote);
%     end
%     
%     
%     for j=1:nodes
%     %j = randi(nodes);   % pick a random node
%     
%         while (z(j) > 1)    % transmission process  - we send equal amounts of y
%             
%             edg = randi(nodes);
%         
%             if (Ad(edg,j) > 0)    % this is an outgoing edge
%            
%                 transm_y(edg,j) = transm_y(edg,j) + floor(y(j)/z(j));
%                 transm_z(edg,j) = transm_z(edg,j) + 1;
%                 y(j) = y(j) - floor(y(j)/z(j));
%                 z(j) = z(j) - 1;
%             
%             end
%         
%         end
%     
%     end
%     
%     
%     
%     for j=1:nodes           % sum incoming
%         
%         y(j) = y(j) + sum(transm_y(j,:));
%         z(j) = z(j) + sum(transm_z(j,:));
%     end
%     
%     
%     transm_z = zeros(nodes,nodes);
%     transm_y = zeros(nodes,nodes);
%         
%     
% %     for j=1:nodes
% %         if (z(j) > 0)
% %             nonz(j) = y(j) / z(j);
% %         end
% %     end
%     
%     first_max = [first_max; max(nonz)];
%     first_min = [first_min; min(nonz(nonz>0))];
%     
%     nonz = zeros(1,nodes);
%     
%     
%     
%     for i=1:nodes
%         for j=1:nodes
%             if (Ad(i,j)==1)
%                 % check for node i the incoming
%                 max_vot(i) = max(max_vot(i),max_vot(j));
%                 min_vot(i) = min(min_vot(i),min_vot(j));
%             end
%         end
%     end
%     
%     flag_to_stop=0;
%     if (mod(k,diam)==0)
%         for akka = 1:diam
%             max_plot = [max_plot; max_vot(1)];
%             min_plot = [min_plot; min_vot(1)];
%         end
%         termination_var = 1;
%         for i=1:nodes
%             if ((max_vot(i) - min_vot(i) > 1))
%                 termination_var = 0;
%             end
%         end
%         if (termination_var==1)
%             flag_to_stop=1;
%             break
%         end
%     end
%     
%     if (flag_to_stop==1)
%         break
%     end
%        
%     
% end

for k=1:max_iter
  % transmission buffers for the nodes
  transm_z = zeros(nodes, nodes);
  transm_y = zeros(nodes, nodes);

  % get the indices that are greater than zero
  z_idcs = find(z > 0);
  % now get the ceil for the y-states based on the predicate above
  y_states(z_idcs) = ceil(y(z_idcs)./z(z_idcs));

  % get current y states of the nodes
  y_states = y_states .* z0;

  % check if the diameter is right for us to check
  if mod(k, diameter) == 1
    vote_states = y ./ z;
    max_votes = ceil(vote_states);
    min_votes = floor(vote_states); 
  end

  % get the nodes that need to be adjusted
  n_idcs = find(z > 1);


  % -- transmission process -- %

  % outgoing transmission
  for idc=1:length(n_idcs)
    out_edge = n_idcs(idc);
    while z(out_edge) > 1
      % get a random node out of the available pool
      node = randi(nodes);

      % now adjust the transmitted weights
      if adjMatrix(node, out_edge) > 0
        flr = floor(y(out_edge) / z(out_edge));
        transm_y(node, out_edge) = transm_y(node, out_edge) + flr;
        transm_z(node, out_edge) = transm_z(node, out_edge) + 1;

        y(out_edge) = y(out_edge) - flr;
        z(out_edge) = z(out_edge) - 1;
      end

    end

  end

  % sum with the incoming transmitted values
  y = y + sum(transm_y, 2);
  z = z + sum(transm_z, 2);

  % -- end transmission process -- %

  % find the min/max votes after values settle
  for h=node_len_array
    [~, cols] = find(adjMatrix(h, :) == 1);
    max_votes(h) = max(max_votes(h), max(max_votes(cols)));
    min_votes(h) = min(min_votes(h), min(min_votes(cols)));
  end

  % check if the termination condition is met
  if mod(k, diameter) == 0
    vote_index = max_votes - min_votes <= epsilon;
    nodes_converged = find(vote_index == 1);
    nodes_diverged = find(vote_index ~= 0);

    a = 1;  
    node_stats(nodes_converged, 1) = min(node_stats(nodes_converged, 1), k);
    node_stats(nodes_converged, 2) = max(node_stats(nodes_converged, 2), k);
    %min(node_stats(:, 1))

    % check if we can terminate, which is when all votes are raised
    if vote_index
      break;
    end
  end
end

y
z

    for j=1:nodes
        if (z(j) >0)
            nodes_states_y(j) = ceil(y(j)/z(j));
        end
    end

% figure
% hold
% title('Quantized Average Consensus via Mass Equalization - Envelope','FontWeight','Normal')
% plot(first_max)
% plot(first_min)
% hold


%first = (first'*z0)'

figure
hold
stairs(ceil(first/quant_step))
title('Load per Node According to Processing Capacity', 'FontWeight','Normal')
ylabel('Node State Variables ($q_j^s[k]$)', 'interpreter', "latex", "fontsize", 18)
xlabel('Number of Iterations ($k$)', 'interpreter', "latex", "fontsize", 18)

C_average = nodes_states_y

initial_aver 


figure
hold
stairs(first2)
title('Load per Processing Cycle + Max/Min plots (Dashed)')
ylabel('Node State Variables ($q_j^s[k]$)', 'interpreter', "latex", "fontsize", 18)
xlabel('Number of Iterations ($k$)', 'interpreter', "latex", "fontsize", 18)
stairs(max_plot, '--')
stairs(min_plot, '--')


