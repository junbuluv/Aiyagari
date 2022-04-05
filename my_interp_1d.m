function val=my_interp_1d(g,xq,f,N_x)
    init=1;
    % See below for the bisection_search() function
    idx=bisection_search(init,g,xq,N_x);
    x_L=g(idx);
    x_H=g(idx+1);
    y_L=f(idx);
    y_H=f(idx+1);
    val=((xq-x_L)*y_H+(x_H-xq)*y_L)/(x_H-x_L);
end

% I recommend adding comments to this to make sure you understand what's
% going on
function idx=bisection_search(init,g,xq,Ng)
   dis=0.0;
   if xq<=g(1)
        idx=1;
    elseif xq>=g(Ng-1)
        idx=Ng-1;
   else
        low=init;
        high=Ng;
        found=false;
        while ~found 
            id = floor((low+high)/2);
         
            dis=xq-g(id);
            if dis<=(g(id+1)-g(id)) && dis>=0.0
                idx=id;
                found=true;
            elseif dis>=(g(id-1)-g(id)) && dis<=0.0 
                idx=id-1;
                found=true;
            elseif dis>(g(id+1)-g(id)) && dis>0.0
                low=id+1;
            elseif dis<(g(id-1)-g(id)) && dis<0.0
                high=id-1;
            end
        end
    end
end