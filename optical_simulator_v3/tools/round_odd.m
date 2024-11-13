function o_num = round_odd(i_num)
    i_num = floor(i_num);

    if mod(i_num , 2) == 0
        o_num = i_num + 1;
    else
        o_num = i_num;
    end  
end

