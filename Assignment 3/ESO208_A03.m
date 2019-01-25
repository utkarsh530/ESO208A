clc;
%Problem 9 of PS 7
h = input('Grid size h?: ');
n = (2/h)+1;
a = @(x) 1-((-(x+3)/(x+1))*h)/2;
b = @(x) (h^2)*((x+3)/(x+1)^2)-2;
c = @(x) 1+((-(x+3)/(x+1))*h)/2;
f  = @(x) (2*(x+1)+3*((x+3)/(x+1)^2))*(h^2);
prompt = input('1.2nd order Backward Difference\n2.2nd order central difference with Ghost Node\nEnter Your Option(1-2): ');
if(prompt == 2)
    
end    