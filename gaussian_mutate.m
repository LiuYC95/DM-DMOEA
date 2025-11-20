function newpoint = gaussian_mutate(newpoint,V,min_range,max_range)
   prob=1/V;
   sigma = (max_range-min_range)./20;
   newparam = min(max(normrnd(newpoint,sigma),min_range),max_range); 
   C = rand(1,V)<prob;
   newpoint(C) = newparam(C);
end
