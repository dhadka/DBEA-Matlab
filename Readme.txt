% This is a code for constrained and unconstrained multi/many objective
% optimization: Authors Md. Asafuddoula, Tapabrata Ray and Ruhul Sarker

% Generates reference directions via systematic sampling.
% Population size is defined based on the number of weight vectors
% Uses SBX and polynomial mutation to create two children, one of which is
% considered.

% The objective functions should be positive and remain positive throughout
% the search domain.

% Intercept update rule: if the norm of (ideal point in current generation
% minus ideal point in the last generation is less than tolerance of 1e-7,
% the intercepts can reduce).

% eps_obj is set to 0: Allign the solutions
% eps_obj can be adaptively choosen in which case there will be compromise
% in the allignment while it may offer a faster convergence for many
% objective optimization problems.

% For details, Please see the problem definition file for creation and addition of new problems 

%Please forecited the following two papers for your reference

%[1] Asafuddoula, A., Ray, T. and Sarker, R., “A decomposition based evolutionary algorithm for many objective optimization 
%    with systematic sampling and adaptive epsilon control,” in Proceedings of the Seventh International Conference on 
%    Evolutionary Multi-Criterion Optimization, vol. 7811 Lecture Notes in Computer Science, pp. 417-427, Springer, 2013. 

%[2] Asafuddoula, M., Ray, T. and Sarker, R., “A decomposition based evolutionary algorithm for many objective optimization” 
%    IEEE Transactions on Evolutionary Computation, In Press, (Accepted 10/06/2014) 


%%USAGE%%

% Add the required path of the nessesary files i.e. path of the core code, path of the problem , for multiobjective mex path

%run the multirun_DBEA.m file.
%for a specific problem set the problem name(s).
multirun_DBEA
%Datafile
 
%1. all datafile contains
        %for i=1 to size of population 
        %alldata=[alldata; gen i pop(i).x pop(i).c pop(i).Fi pop(i).f];
%2. def.stats_all =1 will produce the data for all generations.





