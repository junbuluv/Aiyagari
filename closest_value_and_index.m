function [x, ibest] = closest_value_and_index(xgrid, val)


ibest = closest_index(xgrid, val);

x = xgrid(ibest);
end

