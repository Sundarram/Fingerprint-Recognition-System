function y=minutie(x)
i=ceil(size(x)/2);
if x(i,i)==0;
    y=0;
else
    y=sum(x(:)) - 1;
end