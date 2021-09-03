# extremal_index
### Extremal index estimators: 
    -Intervals
    -Sliding
    -Northrop
    -runs   : n^.*            
    -blocks : n^.*   
    -stable : n^.*/bl      , * in {0.6,0.7}
    
Main function: 
    -plot_estimators(  )  : plots implemented estimators)
    -exIndex(  )              : returns a data.frame with all estimators as a function of k (block length)
        plot_estimators  plots implemented estimators
        
Comments:
     run: plot_estimators(th0=exIndex)    -> If the return of exIndex has been computed previously and you want to plot results.
     run: plot_estimators(sample0,alpha0)  If estimators + plot needs to be computed.
            set ei0 = true value of the extremal index (if known).


