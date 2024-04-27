addpath(pwd)
addpath(strcat(pwd, filesep, 'Methods'));
addpath(strcat(pwd, filesep, 'Problems'));

def.gen = 400;
def.seed = 99;
def.displayIterativePlot = 0;
def.stats_all = 0;
def.distanceConvergencePlot = 0;
def.displayFinalFront = 0;
def.s1 = 12;
def.s2 = 0;
def.nf = 3;
def.problem_name = 'dtlz1';

DBEA(def)