#!/usr/local/bin/Rscript

# this is ..../artistry/papers/current/ann_stapleton_ms/our_parts/simulations/modified_eqn.r
#
# derived from  ../artistry/papers/current/ann_stapleton_ms/our_parts/simulations/flipper.r and
# ../artistry/papers/current/ann_stapleton_ms/our_parts/simulations/sweep.r
#
# this sweeps through parameter space using the simpler version of the
# elliptical paraboloid's eqn:
#
# z = c(ax^2 + by^2) + dx + ey
#
# It also writes out the row and column names for p.csv and *results*
# correctly.
#
# It relies on plot_dir.r to generate the plots.
#
# Kazic, 21.2.2016


# call is ./modified_eqn.r  SWEEP_NUM SWING
#
# e.g., ./modified_eqn.r 54 0.01 &




# need to go back and figure out how to do the rotation matrix for Mo17,
# since the mesh points for the central axis are computed in a completely
# ad hoc way.
#
# Kazic, 24.1.2017



# hard-wire the contour_set variable here
#
# Kazic, 24.1.2017

# contour_set <- "mo"



# clear all variables, just in case, and load packages

rm(list=ls(all=TRUE))





# hard-wired now
#
# Kazic, 7.8.2016

# setwd("/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts/simulations")
setwd("~/me/artistry/papers/current/ann_stapleton_ms/our_parts/simulations")  
print(paste0("working dir is: ",getwd()))
# print(paste0("all fcns before sweep_fcns are: ",ls()))



# source("/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts/simulations/sweep_fcns.r")
source("~/me/artistry/papers/current/ann_stapleton_ms/our_parts/simulations/sweep_fcns.r")
print(paste0("sourcing sweep_fcns.r"))
# print(paste0("all fcns after sweep_fcns are: ",ls()))


# source("/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts/simulations/analysis_fcns.r")







# grab variables passed in from the command line; 0.01 =< swing =< 0.05


args = commandArgs(trailingOnly=TRUE)

# args <- c(65,0.01)


if ( length(args) < 2 ) { stop("modified.r requires the sweep number and the value of the swing.\n") 
} else {
	sweep_num <- args[1]
	swing  <- as.numeric(args[2])
        }








# after all re-running of current sweeps done, think about having initial
# location of the parm_file be the sweep's directory.  Leave as is for
# right now.
#
# Kazic, 3.12.2015


# modified to use partition toggling
#
# Kazic, 7.8.2016



local <- toggle_partition(sweep_num)
wd <- build_wd(sweep_num)




parm_stem <- paste0("parm",sweep_num,".r")
parm_file <- file.path(paste0(local,"/results/parameter_sweeps/",parm_stem))






# read in parameter file, then move it to its sweep dir
#
# see src "test data and basic loops" in paraboloid.org for different structures of varying size

source(parm_file)
file.rename(parm_file,file.path(wd,parm_stem))







# compute number of steps in our sweep = num columns in proxies_matrix;
# write this information to the parm file so that the p.csv doesn't yield
# rows of NAs on reloading

steps <- length(a) * length(b) * length(c) * length(d) * length(e)
date <- date()

write(paste("\n\n\n# date = ",date,"\n# local = ",local,"\n# sweep_num = ",sweep_num,"\n# swing = ",swing,"\n# steps = ",steps),file.path(wd,parm_stem),append=TRUE)
write(paste0("\n\n\n# start sweep: ",Sys.time()),file.path(wd,parm_stem),append=TRUE)






# initialize the proxies_matrix; current rownames are in sweep_fcns.r
#
# added in a matrix to detect rotation of the surface by examining the
# slope and intercept of the central axis of each surface
#
# Stapleton and Kazic, 21.12.2016


# added an extra row for the serial position of the peak in matrix to the
# rot_matrix
#
# Kazic, 12.1.2017

proxies_matrix = matrix(NA,nrow=length(proxies),ncol=steps)
# rot_matrix = matrix(NA,nrow=3,ncol=steps)







# generate parameters matrix and write to the output directory, wd;
# include useful lines from other variables as well

p <- t(expand.grid("a"=a,"b"=b,"c"=c,"d"=d,"e"=e))
rownames(p) <- c("a","b","c","d","e")
colnames(p) <- paste0("s",sweep_num,"z",seq(1,ncol(p)))
write.csv(p,file.path(wd,"p.csv"))
save(p,file=file.path(wd,"p.RData"))


cat(paste0("now computing surfaces for ",steps," steps\n"))






# now, loop over p, compute z, and fill out the proxies_matrix
#
# added in computation of intercepts and slopes of the central axis for
# each surface, writing these out to a separate matrix so as not to break
# any other code.
#
# Stapleton and Kazic, 21.12.2016
#
# added in the serial location of the peak using
# http://stackoverflow.com/questions/6920441/index-values-from-a-matrix-using-row-col-indicies
#
# and
# http://eamoncaddigan.net/r/programming/2015/10/22/indexing-matrices/



for ( r in 1:steps ) {
        print(paste0("r: ",r))
	rows <- length(x)
        z=matrix(NA,nrow=rows,ncol=length(y))
        q <- p[,r]
        for ( i in seq_along(x) ) { for ( j in seq_along(y) ) { z[i,j]= q[3]*(q[1] * (x[i]**2) + q[2] * (y[j]**2))  + (q[4]*x[i] + q[5]*y[j]) } }


        z_proxies <- compute_proxies(z,length(proxies),swing)
        proxies_matrix[,r] <- z_proxies


#        rot_list <- grab_surface_rotation_data(x,y,z,contour_set)
	serial_peak <- z_proxies["row",1] + rows*(z_proxies["col",1] - 1 )
#        rot_matrix[,r] <- c(rot_list[[1]],rot_list[[2]],serial_peak)
                            
        zfile <- paste0("/z",r,".RData",collapse="")
        save(z,file=file.path(wd,zfile))

        rm(z)

        if ( identical(as.numeric(steps %% 500000),0) ) {
                safety_proxies_file <- paste0(wd,"safety_proxies_",steps,".rda")
#                safety_rot_file <- paste0(wd,"safety_rot_",steps,".rda")		
                save(list=c(proxies_matrix),file=safety_proxies_file)
#                save(list=c(rot_matrix),file=safety_rot_file)		
                }


        }




rownames(proxies_matrix) <- rownames(z_proxies)
colnames(proxies_matrix) <- paste0("s",sweep_num,"z",seq(1,steps))
summaries <- summary_statistics(proxies_matrix,wd)




# matrix of intercepts and slopes of central axes; added serial position of
# peak in matrix

# rownames(rot_matrix) <- c("new_intcpt","slope","serial_pk")
# colnames(rot_matrix) <- colnames(proxies_matrix)
# save(rot_matrix,file=file.path(wd,"rot_matrix.rda"))







# plot the 2d plots of the results with plotrix; may want to
# further tune colors or dynamic ranges; good thing we're not red-green color blind
#
# Commented out as it's just a blur any more with so many surfaces!
#
# Kazic, 4.8.2016


# results_title <- paste("sweep ",sweep_num,", proxies' values for all parameter combinatns",sep="",collapse="")
# tresults_title <- paste("sweep ",sweep_num,", mags of proxies' values for all parameter combinatns",sep="",collapse="")


# plot_results(wd,"/results_matrix.png",summaries$results,results_title)
# plot_results(wd,"/tresults_matrix.png",summaries$tresults,tresults_title)



write(paste0("\n\n\n# finish job: ",Sys.time()),file.path(wd,parm_stem),append=TRUE)



# exit nicely when called from the command line!

quit(save="no",status=0)
