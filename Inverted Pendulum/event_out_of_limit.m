function [value, isterminate, direction] = event_out_of_limit(t, x, xlim)
% Terminate the solving if the x exceeds the xlimit
value = 1;
for i=1:length(x)
    if abs(x(i)) > xlim(i)
        value = 0;
    end
end
isterminate = 1;
direction = 0;
end

