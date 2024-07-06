function r=Bernoulli(p)

u = rand(1,1);
if (u < p) 
	r=1;
else
	r=0;
end


