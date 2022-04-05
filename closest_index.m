function val =  closest_index(a,x)

if isempty(a) == true
    error("xGrid is empty in function closest_index.")
end
if isnan(a) == true
    error("xGrid is Nan")
end

idx = find(a,x,"first");
if idx == 1
    val = idx;
end
if idx > length(a)
    val = length(a);
end
if (a(idx) ==x)
    val = idx;
end
if abs(a(idx)-1) < abs(a(idx-1)-x)
    val = idx;
else
    val = idx-1;

end
end
