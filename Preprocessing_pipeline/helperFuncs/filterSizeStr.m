function [filtertag]=filterSizeStr(filter_lo)
    %%Takes in filter cutoff and returns it in string form without period
    
    ones=floor(filter_lo);
    tenths=floor((filter_lo-ones)*10);
    hundredths=floor(((filter_lo-ones)*10-tenths)*10);
    filtertag=[num2str(ones) 'pt' num2str(tenths) num2str(hundredths)];