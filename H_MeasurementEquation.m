function [z] = H_MeasurementEquation(state)
    %Only valid when the sensor does not move with respect to the ECI
    range = state(1:3);
    range_rate = state(4:6);

    range_i = range(1);
    range_j = range(2);
    range_k = range(3);

    range_rate_i = range_rate(1);
    range_rate_j = range_rate(2);
    range_rate_k = range_rate(3);

    if range_i == 0 && range_j == 0
        alpha = atan2(range_rate_j,range_rate_i);
    else
        alpha = atan2(range_j,range_i);
    end
    delta = asin(range_k/norm(range));
    z = [alpha; delta];

end