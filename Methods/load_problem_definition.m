function [prob]=load_problem_definition(def)
funh=str2func(def.problem_name);
try
prob=funh(def.nf);
catch
    error('Problem does not exist...');
end
return