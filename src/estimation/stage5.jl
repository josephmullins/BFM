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


# TODO:
# finish and test code to calculate test score moments in data (above)
# write code to update function with age profiles
# write code to fetch the simulated moments
# write code to get the objective
# write code to minimize and see if it's enough
