function [newpoint]=fixnew(newpoint,min_range,max_range)
newpoint=max(newpoint,min_range);
newpoint=min(newpoint,max_range);
end

  