function retVal = arccos( argument )
        if abs(argument-1) < 1e-12
            argument = 1;
        elseif abs(argument+1) < 1e-12
            argument = -1;
        end
        retVal = acos(argument);
end
