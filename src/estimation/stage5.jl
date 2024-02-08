# add to data frame so we can use a different status variable
function kidmoms_data(d) #<-?
    m = @chain d begin
        @subset :AGE.<=18
        @subset :AGE.>=3
        #@transform :AGE = round(:AGE ./ 4)
        groupby([:DIV,:AGE])
        @combine :S = mean(skipmissing(:AP_raw))
    end
    return m
end
