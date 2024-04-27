function [childpop]=realbinarycrossover(pop,p1,p2,ub,lb)
% parameters related to real SBX crossover
eta_c = 30;
crossoverprob=1.0;
%%
EPS = 1.2e-7;
numVariables = size(pop(p1).x,2);
if (rand <= crossoverprob)
    for i=1:numVariables
        if (rand <= 0.5)
            if (abs(pop(p1).x(i) - pop(p2).x(i)) > EPS)
                if (pop(p1).x(i) < pop(p2).x(i))
                    
                    y1 = pop(p1).x(i);
                    y2 = pop(p2).x(i);
                    
                else
                    y1 = pop(p2).x(i);
                    y2 = pop(p1).x(i);
                    
                end
                yl = lb(i);
                yu = ub(i);
                r = rand;
                beta = 1.0 + (2.0*(y1-yl)/(y2-y1));
                alpha = 2.0 - beta^(-(eta_c+1.0));
                if (r <= (1.0/alpha))
                    
                    betaq =(r*alpha)^(1.0/(eta_c+1.0));
                    
                else
                    
                    betaq = (1.0 / (2.0 - r * alpha))^ (1.0 / (eta_c + 1.0));
                end
                c1 = 0.5*((y1+y2)-betaq*(y2-y1));
                beta = 1.0 + (2.0*(yu-y2)/(y2-y1));
                alpha = 2.0 - beta^(-(eta_c + 1.0));
                if (r <= (1.0/alpha))
                    
                    betaq = (r * alpha)^ (1.0 / (eta_c + 1.0));
                    
                else
                    
                    betaq = (1.0 / (2.0 - r * alpha))^ (1.0 / (eta_c + 1.0));
                end
                c2 = 0.5*((y1+y2)+betaq*(y2-y1));
                if (c1<yl)
                    c1=yl;
                end
                if (c2<yl)
                    c2=yl;
                end
                if (c1>yu)
                    c1=yu;
                end
                if (c2>yu)
                    c2=yu;
                end
                if (rand<=0.5)
                    
                    childpop(1).x(i) = c2;
                    childpop(2).x(i) = c1;
                    
                else
                    
                    childpop(1).x(i) = c1;
                    childpop(2).x(i) = c2;
                    
                end
            else
                
                childpop(1).x(i) = pop(p1).x(i);
                childpop(2).x(i) = pop(p2).x(i);
                
            end
        else
            
            childpop(1).x(i) = pop(p1).x(i);
            childpop(2).x(i) = pop(p2).x(i);
            
        end
    end
else
    for i=1:numVariables
        
        childpop(1).x(i) = pop(p1).x(i);
        childpop(2).x(i) = pop(p2).x(i);
        
    end
end
return