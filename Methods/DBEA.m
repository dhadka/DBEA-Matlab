
%% This is a code for constrained and unconstrained multi/many objective
% optimization: Authors Md. Asafuddoula and Tapabrata Ray

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
%%
function [pop]=DBEA(def)
rand('twister', def.seed);
% Get the problem definition
prob=load_problem_definition(def);
def.nf=prob.nf;
% 2-layer weight direction generation
[def.popSize,w]=direction_vector(def.s1, def.s2, def.nf);
global alldata;
alldata=[];
% Initialize the population of solutions
[pop]=initialize_pop(prob,def);
% Evaluate the performance of these solutions
[pop]=compute_objective(pop,def);
%initialize the ideal and intercept points
idealpoint=ones(1,def.nf)*1e20;intercept=ones(1,def.nf)*1e-20;
for i=1:size(pop,2)
    if(pop(i).Fi==0)
        idealpoint=min(idealpoint,pop(i).f);
        intercept=max(intercept,pop(i).f);
    end
    if(i<=def.popSize)
        w(i,:)=w(i,:)./norm(w(i,:));
    end
end
gen=1;
if def.displayIterativePlot
    show_progress(pop, w ,1, gen,1);
end
D1=[];D2=[];

while(gen<def.gen)
    % Generating the child solutions
    p1=randperm(def.popSize);
    for i=1:def.popSize
        %[child]=evolve(pop(i).alpha,range,prob.nx,ql,def.inc,def);
        [child]=evolve_real(pop,def,i,p1(i),prob);
        [dom_flag]=check_domination(child.f,child.Fi,pop,prob);
        if(dom_flag==1)
            % Computing the intercept and the idealpoint for pop+child
            if(child.Fi==0)
                idealpoint=min(idealpoint,child.f);
                [intercept]=update_intercept(child,pop,intercept);
            end
            % Computing the epsilon threshold for constraints
            eps_con=constraint_approach(pop);
            eps_obj=0;
            % Updating the pop from pop+child
            [pop]=update_population(child,pop,w,intercept,idealpoint,eps_con,eps_obj,prob,def);
        end
    end
 
    if(def.distanceConvergencePlot)
        % calculate the distance measures and store for the final plot
        d1=0;d2=0;
        for j=1:def.popSize
            F1=pop(j).f;
            %distance along the line
            d1=d1+F1*w(j,:)';
            %distance to the segment
            d2=d2+sqrt(norm(F1)^2-pop(j).d1_parent^2);
        end
        D1=[D1 d1/def.popSize];
        D2=[D2 d2/def.popSize];
    end
    if(def.displayIterativePlot)
        clc;
        disp(strcat(num2str(round((gen/def.gen)*100)),'%'));
        
%       show_progress(pop,w,2,gen,2);  % iterative plot with tagged solution
        show_progress(pop,w,2,gen,1);
    else
        clc;
        disp(strcat(num2str(round((gen/def.gen)*100)),'%'));
    end
    
    % Writing the population information to a file
    write_progress(pop,gen,def);
    gen=gen+1;
end
all_f=reshape([pop.f],prob.nf,size(pop,2))';
f=all_f(1:def.popSize,:);
ndffeas=nd_rank_one(f);
dlmwrite('ndfeasibleobj.dat',ndffeas,'precision','%10.4f','delimiter',' ');
% if you wish to see the final front
if(def.displayFinalFront)
    show_progress(pop,w,1,gen,1);
    show_progress(pop,w,2,gen,2);
end
% if you wish to see the d1, d2 convergence
if(def.distanceConvergencePlot)
    figure(3);
    plot(1:gen-1,D1(1:gen-1),'k.-','markersize',8,'linewidth',2); hold on;
    plot(1:gen-1,D2(1:gen-1),'m.--','markersize',8,'linewidth',2); hold on;
    legend('d_1','d_2');
end
return

% Computing the intercept along the objective axes
function [intercept]=update_intercept(child,pop,intercept)
I_max=1e20;
M=size(child.f,2);
f=reshape([pop.f],M,size(pop,2))';
f=[f(reshape([pop.Fi],1,size(pop,2))==0,:);child.f];
% identify M extreme solutions using corner-sort
r=sort_corner(f);
f=f(r(1:2*M),:);
[f_nd,~]=nd_rank_one(f);
m=size(f,2);
if ~isempty(f)
    for i=1:m
        [a(i), b(i)]=max(f_nd(:,i));
    end
    point=f_nd(b,:);
    point=unique(point,'rows');
    if(size(point,1)<m )
        I_max=a;
    else
        % compute the intercepts of the identified hyperplane
        numerator=det(point);right_side=ones(m,1)*numerator;
        normal=point\right_side; % inv(point)*right_side
        I_max=zeros(1,m);
        for i=1:m
            I_max(i)=numerator/normal(i);
        end
        % condition if intercepts fail
        for i=1:m
            if(I_max(i)<0 || isnan(I_max(i))||isinf(I_max(i)))
                I_max(i)=a(i);
            end
        end
    end
end
intercept=I_max;

return


% Updating the population based on a child pop and a pop.
function [pop]=update_population(child,pop,w,intercept,idealpoint,eps_con,eps,prob,def)
for i=1:size(pop,2)
    F1=(pop(i).f-idealpoint)./(intercept-idealpoint);
    pop(i).F=F1;
end
success=0;
% normalize this solution is a dominated solution
F2=(child.f-idealpoint)./(intercept-idealpoint);
%random shuffled order to give the fair chance to compare all
order=randperm(def.popSize);
for k=1:def.popSize
    j=order(k);
    F1=(pop(j).F);
    if(child.Fi<eps_con && pop(j).Fi<eps_con || (child.Fi == pop(j).Fi))
        
        %distance along the direction vector
        d1_parent=F1*w(j,:)';
        d1_child=F2*w(j,:)';
        %distance perpendicular to the direction vector
        d2_parent=sqrt(norm(F1)^2-d1_parent^2);
        d2_child=sqrt(norm(F2)^2-d1_child^2);
        pop(j).d1_parent=d1_parent;
        pop(j).d2_parent=d2_parent;
        
        
        [b]=compare_objective(d1_child,d2_child,d1_parent,d2_parent,eps);
        if(b==1)
            %Child winner
            pop(j).x=child.x;
            pop(j).f=child.f;
            pop(j).c=child.c;
            pop(j).Fi=child.Fi;
            pop(j).d1_parent=d1_child;
            pop(j).d2_parent=d2_child;
            pop(j).F=F2;
            success=success+1;
        end
    end
    if(child.Fi < pop(j).Fi)
        %Child winner
        pop(j).x=child.x;
        pop(j).f=child.f;
        pop(j).c=child.c;
        pop(j).Fi=child.Fi;
        pop(j).F=F2;
        success=success+1;
    end
    if(success==1)
        break;
    end
end

return


function [b]=compare_objective(d1_child,d2_child,d1_parent,d2_parent,eps)
b=0;
if(d2_child==d2_parent) || ((d2_child<eps) && (d2_parent<eps))
    if(d1_child<d1_parent)
        %Child winner
        b=1;
    end
elseif (d2_child<d2_parent)
    %Child winner
    b=1;
end
return

% Computing the threshold for constrained problems
function [e]=constraint_approach(pop)
feasible=1;cv_max=0;popsize=size(pop,1);
for i=1:popsize
    if(pop(i).Fi==0)
        feasible=feasible+1;
    end
    cv_max=cv_max+pop(i).Fi;
end
cv_max=cv_max/popsize;
feasibility=feasible/popsize;
e=feasibility*cv_max;
return


function [child]=evolve_real(pop,def,p1,p2,prob)
lb=prob.range(:,1)';ub=prob.range(:,2)';
[childpop]=realbinarycrossover(pop,p1,p2,ub,lb);
rate=1/prob.nx;
[child.x]=realmutate(childpop(1).x,rate,prob.range);
func=str2func(def.problem_name);
[f,c,x] = func(prob.nf,child.x);
c=max(c,0);    % max of violation and 0 is equality value
child.f=f;
child.c=c;
child.x=x;
child.Fi=sum(c);
return

% Evaluating a population of solutions
function [pop]=compute_objective(pop,def)
func=str2func(def.problem_name);
for i=1:size(pop,2)
    [f,c,x]=func(def.nf,pop(i).x);
    c=max(c,0);    % max of violation and 0 is equality value
    pop(i).f=f;
    pop(i).c=c;
    pop(i).x=x;
    pop(i).Fi=sum(c);
    pop(i).d1_parent=0;
    pop(i).d2_parent=0;
    pop(i).F=f;
end

return

% Initializing a population of solutions
function [pop]=initialize_pop(prob,def)
lb=prob.range(:,1)';ub=prob.range(:,2)';
x = bsxfun(@plus,lb,bsxfun(@times,lhsdesign(def.popSize,prob.nx),(ub-lb)));
for i=1:def.popSize
    pop(i).x=x(i,:);
end
return


% Generating the direction vectors
%Function to generate uniformly distributed weight vectors
function [n,w] = direction_vector(s1,s2,m)
[n,w] = gen_weight(s1,m);
if s2 > 0
    [n2,w2] = gen_weight(s2,m);
    n = n+n2;
    w = [w;w2*0.5+(1 - 0.5)/(m)];
end
return

% Generating the direction vectors
function [n,w]=gen_weight(p,M)
NumPoints((1:M-1))=p;
% partitioned into the first objective
Beta1 = (0:(NumPoints(1)))'/(NumPoints(1));
Beta = Beta1;
for i = 2:M-1
    ki = round((1-sum(Beta1,2))*(NumPoints(i)));
    Beta = [];
    for j =1:size(Beta1,1)
        BetaVec = (0:ki(j))'/(NumPoints(i));
        Beta = [Beta; [repmat(Beta1(j,:), length(BetaVec),1) BetaVec] ];
    end
    Beta1 = Beta;
end
w= [Beta (1 - sum(Beta1,2))]; 
n=size(w,1);
return



% Function to write the all.dat file
function write_progress(pop,gen,def)
global alldata;
if(def.stats_all==1)
    for i=1:size(pop,2)
        alldata=[alldata; gen i pop(i).x pop(i).c pop(i).Fi pop(i).f];
    end
elseif(gen==def.gen-1)
    for i=1:size(pop,2)
        alldata=[alldata; gen i pop(i).x pop(i).c pop(i).Fi pop(i).f];
    end
end

if(gen==def.gen-1)
    save(strcat(def.problem_name,'_DBEAp_all','.mat'),'alldata');
end
return

% Function to show the solutions over generations
function show_progress(pop, w, fig_id, gen, flag)
drawnow;
figure(fig_id);
m=size(pop(end).f,2);
feasible=[];
infeasible=[];
%scale=max(reshape([pop.f],m,size(pop,2))');
for i=1:size(pop,2)
    if(pop(i).Fi==0)
        if(flag==1)
            feasible=[feasible;pop(i).f];
        else
            feasible=[feasible;pop(i).F];
        end
    else
        if(flag==1)
            infeasible=[infeasible;pop(i).f];
        else
            infeasible=[infeasible;pop(i).F];
        end
    end
end
if(flag==2)
    for i=1:size(pop,2)
        
        if(m==2)
            O= [0 0];
        elseif(m==3)
            O=[0 0 0];
        end
        P=w(i,:);
        d1=pop(i).F*w(i,:)';
        cord=d1*w(i,:);
        if(m==2)
            plot([pop(i).F(1,1) cord(1,1)],[pop(i).F(1,2) cord(1,2)],'color','m','Marker','.','LineStyle','-','linewidth',2);hold on;
            plot([O(1,1) P(1,1)],[O(1,2) P(1,2)],'color','k','Marker','.','LineStyle','-','linewidth',2);
        elseif(m==3)
            plot3([pop(i).F(1,1) cord(1,1)],[pop(i).F(1,2) cord(1,2)],[pop(i).F(1,3) cord(1,3)],'color','m','Marker','.','LineStyle','-');hold on;
            plot3([O(1,1) P(1,1)],[O(1,2) P(1,2)],[O(1,3) P(1,3)],'color','k','Marker','.','LineStyle','-');
            view(111,28);
        end
        
    end
    hold on;
end
if(~isempty(feasible) && ~isempty(infeasible))
    if(m==2)
        plot(infeasible(:,1), infeasible(:,2), 'r*', ...
            feasible(:,1), feasible(:,2), 'b*');
        grid on;
    elseif(m==3)
        plot3(infeasible(:,1), infeasible(:,2), infeasible(:,3),'r*', ...
            feasible(:,1), feasible(:,2), feasible(:,3),'b*');
        grid on;
    end
elseif(~isempty(infeasible))
    if(m==2)
        plot(infeasible(:,1), infeasible(:,2), 'r*');
        grid on
    elseif(m==3)
        plot3(infeasible(:,1), infeasible(:,2),infeasible(:,3), 'r*');
        grid on;
    end
elseif(~isempty(feasible))
    if(m==2)
        plot(feasible(:,1), feasible(:,2), 'bo','linewidth',2);
        grid on;
    elseif(m==3)
        plot3(feasible(:,1), feasible(:,2),feasible(:,3), 'b*');
        grid on;
    end
end
title(gen);
hold off;
xlabel('f_1');
ylabel('f_2');
return

% Checking if this child solution is dominated by any solution in the
% population.
function [dom_flag]=check_domination(f,Fi,pop,prob)
dom_flag=1;
if(Fi==0)
    for i=1:size(pop,2)
        if(pop(i).Fi==0)
            diff=pop(i).f<f;
            if sum(diff)==prob.nf
                dom_flag=0;
                break;
            end
        end
    end
end
return







