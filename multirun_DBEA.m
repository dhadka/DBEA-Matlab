function multirun_DBEA
path=pwd;
cd(path);
addpath(pwd);
addpath(strcat(pwd,filesep,'Methods'));
addpath(strcat(pwd,filesep,'Problems'));
path = pwd;
path11=strcat(path,filesep);
dir=strcat(path11,'Results');
mkdir(dir);
% please specify the number of runs  
nruns=1;
% please set the problem name here
%pb={'dtlz1','dtlz2'};
pb={'dtlz1'};
% number of objectives to be considered
obj=[3, 3];
% number of generations to evolve
gen=[400, 400];
% spacing in upper and lower level
s1=[12, 12]; % upper level spacing
s2=[0, 0]; % lower level spacing
% initial seed to be considered
seed=99;count=1;
for j=1:length(pb)
    dir1=strcat(dir,filesep,pb{j});
    mkdir(dir1);
    cd(dir1);
    for i=1:1:nruns
        dir2=strcat(dir1,filesep,'run-',num2str(i));
        mkdir(dir2);
        cd(dir2);
        def.gen=gen(j);
        seed=seed+j;
        def.seed=seed;
        % if one wish to see the progressive plot in each iteration
        def.displayIterativePlot=0;
        def.stats_all=0;
        % if one wish to see the d1, d2 convergene plot
        def.distanceConvergencePlot=0;
        % if one wish to see the final front
        def.displayFinalFront=0;
        % direction
        def.s1=s1(j);
        def.s2=s2(j);
        def.nf=obj(j);
        def.problem_name=pb{j};
        save('Params.mat', 'def');
        path1{count} = dir2; 
        count = count+1;
        cd ..;
    end
     cd(path);
end

for i=1:length(path1)
    cd(path1{i});
    disp(strcat('Running -> ',path1{i}));
    param=load('Params.mat');
    tic; 
    DBEA(param.def); 
    toc;
    cd(dir);
end
cd(path);
return

