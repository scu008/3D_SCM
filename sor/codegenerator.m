function [codevector] = codegenerator(rate, cl)
if rate == 1/2
    if cl == 3
        codevector=[7 5];
    elseif cl == 4
        codevector=[17 15];
    elseif cl == 5
        codevector=[35 23];
    elseif cl == 6
        codevector=[75 53];
    elseif cl == 7
        codevector=[171 133];
    elseif cl == 8
        codevector=[171 247];
    elseif cl == 9
        codevector=[753 561];
    end
elseif rate == 1/3
    if cl == 3
        codevector=[7 7 5];
    elseif cl == 4
        codevector=[17 15 13];
    elseif cl == 5
        codevector=[37 33 25];
    elseif cl == 6
        codevector=[75 53 47];
    elseif cl == 7
        codevector=[171 165 133];
    elseif cl == 8
        codevector=[367 331 225];
    end
end
