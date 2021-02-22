close all
clear all
clc

%largeconstant=10^10                %%% WORK ALSO FOR RING
connectivity=0.5;                   
nodes=100;
iter=200;

Ad = zeros(nodes,nodes);

for i=1:nodes
    for j=1:nodes
        if rand()<connectivity
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



%%% CONSTRAINT FOR WORKING: SUM(Y0) > SUM(Z0)
quant_step=10000;   % we can multiply by max of z0

y0=quant_step*randi(50,nodes,1);
z0=ones(nodes,1);      % z0 is the total capacity of the network
for i=1:nodes
    if (i<10)
        z0(i)=100;
    else
        z0(i)=500;
    end
end
y=y0;
z=z0;



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




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Ad = [1 0 1 0; 1 1 1 0; 1 0 1 1; 0 1 0 1];   % CDC paper plot
%  nodes = 4;
%  y0 = [5; 3; 7; 2];
%  z0 = [1; 1; 1; 1];
%  y=y0;
%  z=z0;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


passes = ones(nodes,1);


first_max = [];
first_min = [];

max_plot = [];
min_plot = [];


node_states_no_oscill = y0./z0;   % WE PLOT THESE STATES FOR NO OSCILLATIONS


states_for_vote = ((y+node_states_no_oscill)./(z+1));
max_vot = ceil(states_for_vote)';
min_vot = floor(states_for_vote)';


nonz = zeros(1,nodes);

flagg=0;

for k=1:iter
    
    if (mod(k,diam)==1)
        states_for_vote = ((y+node_states_no_oscill)./(z+1));
        max_vot = ceil(states_for_vote)';
        min_vot = floor(states_for_vote)';
    end
    
    
    if flagg==0
        flagg=1;
    else
        node_states_no_oscill2 = node_states_no_oscill.*z0;
        first2 = [first2; (node_states_no_oscill2)'];
        first = [first; (node_states_no_oscill)'];
    end
    % node_states_no_oscill --- THIS IS ONLY THE MASS THAT IS KEPT INSIDE
    
    for j=1:nodes
    %j = randi(nodes);   % pick a random node
    
    if (z(j) > 0)
        y(j) = y(j) + node_states_no_oscill(j);
        aaa = ceil((y(j))/(1+z(j)));
        node_states_no_oscill(j) = aaa;
        y(j) = y(j) - aaa;
        clear aaa
    end
    
    while (z(j) > 0)    % transmission process  - we send equal amounts of y
        
        edg = randi(nodes);
        
        if (Ad(edg,j) > 0)    % this is an outgoing edge
           
            transm_y(edg,j) = transm_y(edg,j) + floor(y(j)/z(j));
            transm_z(edg,j) = transm_z(edg,j) + 1;
            y(j) = y(j) - floor(y(j)/z(j));
            z(j) = z(j) - 1;
            
        end
    end
    
    end
    
    
    for j=1:nodes           % sum incoming
        
        y(j) = sum(transm_y(j,:));
        z(j) = sum(transm_z(j,:));
    end
    
    
    transm_z = zeros(nodes,nodes);
    transm_y = zeros(nodes,nodes);
        
    
    for j=1:nodes
        if (z(j) > 0)
            nonz(j) = y(j) / z(j);
        end
    end
    
    first_max = [first_max; max(nonz)];
    first_min = [first_min; min(nonz(nonz>0))];
    
    nonz = zeros(1,nodes);
    
    
    for j=1:nodes
        if (z(j) >0)
            nodes_states_y(j) =  round(y(j)/z(j));
        end
    end
    
    
    for i=1:nodes
        for j=1:nodes
            if (Ad(i,j)==1)
                % check for node i the incoming
                max_vot(i) = max(max_vot(i),max_vot(j));
                min_vot(i) = min(min_vot(i),min_vot(j));
            end
        end
    end
    
    flag_to_stop=0;
    if (mod(k,diam)==0)
        for akka = 1:diam
            max_plot = [max_plot; max_vot(1)];
            min_plot = [min_plot; min_vot(1)];
        end
        termination_var = 1;
        for i=1:nodes
            if ((max_vot(i) - min_vot(i) > 1))
                termination_var = 0;
            end
        end
        if (termination_var==1)
            flag_to_stop=1;
            break
        end
    end
    
    if (flag_to_stop==1)
        break
    end
            
%      nodes_states_y
    % first = [first; (y./z)'];
    % first = [first; (nodes_states_y)'];
    if abs(sum(y+node_states_no_oscill) - sum(y0)) >= 1e-10
        k
        error("this should not happen");
    end
end


% figure
% hold
% title('Quantized Average Consensus via Mass Equalization - Envelope','FontWeight','Normal')
% plot(first_max)
% plot(first_min)
% hold

node_states_no_oscill2 = node_states_no_oscill.*z0;
first2 = [first2; (node_states_no_oscill2)'];
first = [first; (node_states_no_oscill)'];


C_average = node_states_no_oscill

initial_aver


first2=first2/quant_step;
% first2 = round(first2);   % because sometimes the final values are real
% because we increase quantization level and we divide. 

figure
hold
stairs(first)
title('Load per Processing Cycle + Max/Min plots (Dashed)','FontWeight','Normal')
ylabel('Node State Variables (q_j^s[k])')
xlabel('Number of Iterations (k)')
stairs(max_plot, '--')
stairs(min_plot, '--')


figure
hold
stairs(first2)
title('Load per Node According to Processing Capacity','FontWeight','Normal')
ylabel('Node State Variables (q_j^s[k])')
xlabel('Number of Iterations (k)')

