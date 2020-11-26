clc; clear 
p_av=[];
n_av=[];
t_av=[];
ran=16;

for name=["whole","localcomp", "reduce", "getdata"]
for i=2:64
fname=name+i+".out";
data=load(fname);
N=size(data(:,1));
p_av(i) = data(1,1);
n_av(i) = data(1,2);
t_av(i) = mean(data(:,3));
end
plot(p_av(2:ran), t_av(2:ran),'*-');
hold on;
end
%------Set the axes names:--------------
lx="Number of processors";
ly="Average Elapsed time (s)";
%-------Set Plot Title:------------------
tit= "Average Elapsed time (s) as a function of the number of processors. N=1000";
%------------------Legend-----------
legend("whole","localcomp", "reduce", "getdata")


%errorbar(p_av(2:16), t_av(2:16), stderror(2:16), 'LineStyle','none', 'Color', 'k','linewidth', 1)
%legend('Whole Program','Time for adding Local Summations')
xlabel(lx), ylabel(ly)
title(tit)