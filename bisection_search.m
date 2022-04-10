function idx=bisection_search_mod(init,g,xq,Ng)
   dis = 0.0;
   %check if the query value is actually on the grid i.e. less than first
   %value, keep first value; greater than last value, keep last value
   if xq <= g(1)
        idx = 1;
    elseif xq >= g(Ng)
        idx = Ng;
   else
       %initializing search with low and high values
        low = init;
        high = Ng;
        found = false;
        %while it's not found, stay in loop (once found, exit)
        while ~found 
            %mid-index on grid
            id = floor((low+high)/2);
            %distance of query point and value of mid-index
            dis=xq-g(id);
            
            %distance is positive, and point is in between id and id+1, so
            %we've found lower and upper bound
            if dis<=(g(id+1)-g(id)) && dis>=0.0
                idx=id;
                found=true;
            
            %distance is negative, and point is in between id-1 and id, so
            %we've found lower and upper bound 
            elseif dis>=(g(id-1)-g(id)) && dis<=0.0 
                idx=id-1;
                found=true;
                
            %distance is positive, and we're above prior midpoint so search
            %up
            elseif dis>(g(id+1)-g(id)) && dis>0.0
                low=id+1;
                
            %distance is negative, and we're below prior midpoint, so
            %search below
            elseif dis<(g(id-1)-g(id)) && dis<0.0
                high=id-1;
            end
        end
    end
end