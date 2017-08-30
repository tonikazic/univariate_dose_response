# this is ../artistry/papers/current/ann_stapleton_ms/our_parts/simulations/sweep_fcns.r

# a subroutine library for sweep.r
#
# Kazic, 2.11.2015
#
# added plot_em to speed up recalculations
#
# Kazic, 25.11.2015

# rm(list=ls(all=TRUE))


# source("/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts/simulations/analysis_fcns.r")
# source("~/me/artistry/papers/current/ann_stapleton_ms/our_parts/simulations/analysis_fcns.r")
# print(paste0("sourcing analysis_fcns.r"))
# print(paste0("all fcns after analysis_fcns are: ",ls()))





############ directory and file management ##############

# not sure if setting or calling it is faster, opt for setting

# modified to include toggling of partitions.  This will break some older
# subroutines.
#
# Kazic, 6.8.2016


build_wd <- function(sweep_num) {

        local <- toggle_partition(sweep_num)

        sweep_dir <- paste("sweep",sweep_num,sep="",collapse="")
        wd <- paste(local,"results/parameter_sweeps",sweep_dir,sep="/",collapse="/")
	if ( ( !file_test("-d",wd) && ( !file_test("-f",wd) ) ) ) { dir.create(file.path(wd)) }

        return(wd)
	}



# a hack, since it seems R doesn't look at the function's 
# arity when looking for functions.  Since analysis_fcns.r is
# sourced after sweep_fcns.r in modified_eqn.r, calls to build_wd/2
# might confuse the name space?
#
# Kazic, 26.9.2016

build_wd_dbl <- function(local,sweep_num) {

#        local <- toggle_partition(sweep_num)

        sweep_dir <- paste("sweep",sweep_num,sep="",collapse="")
        wd <- paste(local,"results/parameter_sweeps",sweep_dir,sep="/",collapse="/")
	if ( ( !file_test("-d",wd) && ( !file_test("-f",wd) ) ) ) { dir.create(file.path(wd)) }

        return(wd)
	}



















############## proxy functions #############


# list the proxies!
#
# current as of 30.11.2015
#
# added the values at the corners
#
# Kazic, 17.1.2016

proxies <- c("reflectn_min","area_min","curvature_min",
             "reflectn_max","area_max","curvature_max",
#
             "level0.25","level0.50","level0.75",
#
             "surface_vol",
#	     
	     "displ(b)","displ(t)","displ(r)","displ(l)",
             "slope(b)","slope(t)","slope(r)","slope(l)",
#	     
     	     "z(lr)","z(ll)","z(ur)","z(ul)",
#
             "minn","maxn","min_z_range","max_range",
             "z-lr","z-ll","z-ur","z-ul",
             "col","row","z_min","z_max")







# maxima, minima, and ranges
#
# hmm, found bug: s9z705 has:
#
# Browse[2]> corners
#           [,1]      [,2]      [,3]      [,4]
# [1,] -44.20917 -43.98417 -19.00917 -18.78417
#
# Browse[2]> min(z)
# [1] -44.20917
# Browse[2]> max(z)
# [1] 2.879688
#
# but e is
#
# Browse[2]> e
#                   [,1]
# z_max         2.879688
# z_min       -44.209170
# row         123.000000
# col          16.000000
# ul_z         47.088857
# ur_z         46.863857
# ll_z         21.888857
# lr_z         21.663857
# max_range    47.088857
# min_z_range  21.663857
# maxn          1.000000
# minn          1.000000
#
#
# hmmmmm! ul_z, ur_z, ll_z, and lr_z are the displacements of those
# corners relative to the peak!  So these are correctly calculated, but not
# what I expect by now.  So just report the corner values separately,
# increasing the number of proxies; and leave these ranges alone.
#
# Kazic, 17.1.2016



extrema <- function(z,maxl,corners) {

# minima and maximum
# for the maximum value, take the middlish one in case of multiples

        if ( identical(nrow(maxl),as.integer(1),num.eq=TRUE) ) {
	        z_max <- max(z)
		midmax <- which(z == max(z), arr.ind = TRUE)
	} else {
	        m = floor(nrow(maxl)/2)              # middle of peak if a very flat one
                midmax <- t(as.matrix(maxl[m,]))
                z_max <- z[midmax]
                }          

        
        minl = which(z == min(z), arr.ind = TRUE)       # returns [row,col] of global minimum value
        z_min = min(z)                            # minimum value; take first in case of multiples





# maximum position and number of minima (the global minimum, should be 1, but it may have
# very flat curvature)

        row = midmax[1,1]
        col = midmax[1,2]
        maxn = length(maxl)/2
        minn = length(minl)/2




# ranges

        max_range = abs(max(range(z)) - min(range(z)))
        all_z_ranges = z_max - corners;
        min_z_range = min(all_z_ranges)


        e <- matrix(NA,12,1)


# as.integer not coercing here; order correct
#
# Kazic, 17.1.2016

        e[,1] <- c(z_max,z_min,as.integer(row),as.integer(col),
                   all_z_ranges[1],all_z_ranges[2],all_z_ranges[3],all_z_ranges[4],
                   max_range,min_z_range,as.integer(maxn),as.integer(minn))

		    
        rownames(e) <- c("z_max","z_min","row","col",
	                 "z-ul","z-ur","z-ll","z-lr",
                         "max_range","min_z_range","maxn","minn")

        return(e)
        }








# absolute value of the displacements of peak, projected to the four edges,
# relative to the slope of each edge at the peak's position

displcmts <- function(z,row,col,rows,cols,z_max,corners) {


# hmmm, this is the silly step.  I'm just rescaling z_max by its relative
# position in the surface! and the perfect parallels come from repeating
# the entries . . . ah me.
#
# Kazic, 17.1.2016
#
#        projectns <- (c(col_pos, col_pos, row_pos, row_pos)*z_max)


# slopes of the bounding edges
#
# these are correct
#
# Kazic, 17.1.2016

        displmts <- matrix(NA,8,1)
        displmts[1,1] <- as.numeric(corners[1] - corners[3])             # left 
        displmts[2,1] <- as.numeric(corners[2] - corners[4])             # right
        displmts[3,1] <- as.numeric(corners[1] - corners[2])             # top
        displmts[4,1] <- as.numeric(corners[3] - corners[4])             # bottom



# displacemnts of projected peaks over each edge        
#
# ah, the other half of the silliness:
#
# displmts[5:8,1] <- abs(displmts[1:4,1] - projectns)
#
# what I want is the value of z at each of those positions on the edges ---
#
#  -----------------------------------------------------
#  |                (1,col)                            |
#  |(row,1)         (row,col)            (row,cols)    |
#  |                (rows,col)                         |
#  -----------------------------------------------------
#
# (hooray for ascii art)
#
# z[row,1]   z[row,cols]   z[1,col]    z[rows,col]
#   left        right         top        bottom
#
# Kazic, 17.1.2016

        displmts[5,1] <- as.numeric( abs(z_max - z[row,1]))
        displmts[6,1] <- as.numeric( abs(z_max - z[row,cols]))
        displmts[7,1] <- as.numeric( abs(z_max - z[1,col]))
        displmts[8,1] <- as.numeric( abs(z_max - z[rows,col]))


        rownames(displmts) <- c("slope(l)","slope(r)","slope(t)","slope(b)",
	                        "displ(l)","displ(r)","displ(t)","displ(b)")

        return(displmts)
        }






# level sets
#
# after a long detour into connected components
# (see rmk "level sets, counting, and convex covers" in paraboloid.org),
# we just count the cells for now, maybe someday do convex covers
#
# count from the maximum downward
#
# allow a little swing for now, figure out a nice way to adjust so we get more of a curve, not a volume
#
# ca. December, 2015
#
#
# Hmmm, it would be better to calculate the level_sets relative to an
# absolute measure --- the maximum of B73, which is the highest overall ---
# rather than the z_max of each surface. See
# /Volumes/b/artistry/papers/current/ann_stapleton_ms/first_manuscript/working_draft/stuff_for_figs.org
# at "looking at the category 4 guys in sweep 46" for the thinking here.
#
# Kazic, 10.4.2016
#

level_sets <- function(z,z_max,swing) {

#        mask0.75 <- z > z_max*(0.75 - swing) & z < z_max*(0.75 + swing)
#        mask0.50 <- z > z_max*(0.50 - swing) & z < z_max*(0.50 + swing)
#        mask0.25 <- z > z_max*(0.25 - swing) & z < z_max*(0.25 + swing)

        mask0.75 <- z > 37.615*(0.75 - swing) & z < 37.615*(0.75 + swing)
        mask0.50 <- z > 37.615*(0.50 - swing) & z < 37.615*(0.50 + swing)
        mask0.25 <- z > 37.615*(0.25 - swing) & z < 37.615*(0.25 + swing)

        level0.75 <- length(which(mask0.75,arr.ind = FALSE,useNames = FALSE))
        level0.50 <- length(which(mask0.50,arr.ind = FALSE,useNames = FALSE))
        level0.25 <- length(which(mask0.25,arr.ind = FALSE,useNames = FALSE))


        return(t(data.frame("level0.75"=as.integer(level0.75),"level0.50"=as.integer(level0.50),"level0.25"=as.integer(level0.25))))
        }










# curvature, areas, and reflections around the peak


# per sullivan2006 and sullivan2007, discrete curvature is 2*pi - \Sigma_i^n \theta_i,
# where \theta is the angle at the vertex of the peak.  This is the curvature at the point (the peak),
# NOT the Gaussian curvature over the whole surface (if I correctly understand Sullivan; not 
# sure I do).
#
# We might compute the area of all the triangles that subtend the peak, and then divide the
# curavture by that.  I'm not sure that would be very useful, though.  But why not.
#
# in a crude approximation of the middle of the peak's area for now,
# floor rounds down to nearest integer, which is 0 if nrow == 1

# called to also find the curvature about the min
#
# Kazic, 20.1.2016

find_curvature <- function(point,z) {
        r <- recenter_matrix(point,z);
        angle_area_mat <- compute_angles(r$local_mat);
        curvature <- 2*pi - sum(angle_area_mat[1,]);
        areas <- sum(angle_area_mat[2,]);

        return(t(data.frame("curvature"=curvature,"areas"=areas,"reflectn"=r$reflectn)));
        }





# if the peak is not in the interior of the z matrix, then there are
# two cases:  the peak is smack on a corner, or it is displaced one cell
# from a corner.  If on a corner, then we must also reflect the cell at the opposite
# diagonal corner.  Otherwise, just reflecting the row or column is enough.
#
# tested using z <- B in the test data above, and the peak at [3,3]
# all rewriting checks out! results are in [[file:recentering_test.org]].
#
# Kazic, 21.10.2015


# just rechecked, coordinates of first two rows of local_mat correct
#
# Kazic, 20.1.2016

recenter_matrix <- function(point,z) {

        local_mat <- rbind(c(0, -1,-1,-1, 0,1, 1,1,   0, -1), 
                          c(0, -1,0,1,   1,1, 0,-1, -1, -1), 
                          c(NA,NA,NA,NA,NA,NA,NA,NA,NA,NA));

        crow = as.numeric(point[1,1])
        ccol = as.numeric(point[1,2])

        local_mat[3,1] = z[point];


# in the interior, no reflectn

        if ( ( crow > 1 ) && ( crow < nrow(z) ) && ( ccol > 1 ) && ( ccol < ncol(z) ) ) { 
                local_mat[3,2:4] = z[crow-1,c(ccol-1,ccol,ccol+1)];
                local_mat[3,5:6] = z[c(crow,crow+1),ccol+1];
                local_mat[3,7:8] = z[crow+1,c(ccol,ccol-1)];
                local_mat[3,9] = z[crow,ccol-1];
                reflectn = 0;


# in a corner

        } else if ( ( ( crow == 1 ) || ( crow == nrow(z) ) ) && ( ( ccol == 1 ) || ( ccol == ncol(z) ) ) ) { 
                r <- reflect_corner(crow,ccol,local_mat,z);
		local_mat <- r$local_mat;
		reflectn <- r$reflectn;


# on an edge

        } else if ( ( ( ( crow == 1 ) || ( crow == nrow(z) ) ) && ( ( ccol > 1 ) && ( ccol < ncol(z) ) ) )
                  || ( ( ( crow > 1 ) && ( crow < nrow(z) ) ) && ( ( ccol == 1 ) || ( ccol == ncol(z) ) ) ) ) {
                r <- reflect_edge(crow,ccol,local_mat,z);
		local_mat <- r$local_mat;
		reflectn <- r$reflectn;
        } 

        local_mat[3,10] = local_mat[3,2];

        return(list("local_mat"=local_mat,"reflectn"=reflectn));
        }






######## reflections ############

reflect_corner <- function(crow,ccol,local_mat,z) {
        if ( crow == 1 ) {
                if ( ccol == 1 ) {
                        local_mat[3,2]   = z[crow+1,ccol+1];
                        local_mat[3,3:4] = z[crow+1,c(ccol,ccol+1)];
                        local_mat[3,5:6] = z[c(crow,crow+1),ccol+1];
                        local_mat[3,7]   = z[crow+1,ccol];
                        local_mat[3,9:8] = z[c(crow,crow+1),ccol+1];
                } else {
                        local_mat[3,4]   = z[crow+1,ccol-1];
                        local_mat[3,2:3] = z[crow+1,c(ccol-1,ccol)];
                        local_mat[3,5:6] = z[c(crow,crow+1),ccol-1];
                        local_mat[3,8:7] = z[crow+1,c(ccol-1,ccol)];
                        local_mat[3,9]   = z[crow,ccol-1];
                        }
         } else {
                if ( ccol == 1 ) {
                        local_mat[3,8]   = z[crow-1,ccol+1];
                        local_mat[3,2]   = z[crow-1,ccol+1];
                        local_mat[3,3:4] = z[crow-1,c(ccol,ccol+1)];
                        local_mat[3,5]   = z[crow,ccol+1];
                        local_mat[3,7:6] = z[crow-1,c(ccol,ccol+1)];
                        local_mat[3,9]   = z[crow,ccol+1];
                } else {
                        local_mat[3,6]   = z[crow-1,ccol-1];
                        local_mat[3,2:3] = z[crow-1,c(ccol-1,ccol)];
                        local_mat[3,4:5] = z[c(crow-1,crow),ccol-1];
                        local_mat[3,8:7] = z[crow-1,c(ccol-1,ccol)];
                        local_mat[3,9]   = z[crow,ccol-1];
                        }
         }
         reflectn = -10;

         return(list("local_mat"=local_mat,"reflectn"=reflectn));
         }







reflect_edge <- function(crow,ccol,local_mat,z) {
        if ( crow == 1 ) {
                local_mat[3,2:4] = z[crow+1,c(ccol-1,ccol,ccol+1)];
                local_mat[3,5:6] = z[c(crow,crow+1),ccol+1];
                local_mat[3,8:7] = z[crow+1,c(ccol-1,ccol)];
                local_mat[3,9]   = z[crow,ccol-1];

        } else if ( crow == nrow(z) ) {
                local_mat[3,2:4] = z[crow-1,c(ccol-1,ccol,ccol+1)];
                local_mat[3,5]   = z[crow,ccol+1];
                local_mat[3,8:6] = z[crow-1,c(ccol-1,ccol,ccol+1)];
                local_mat[3,9]   = z[crow,ccol-1];

        } else if ( ccol == 1 ) {
                local_mat[3,2]   = z[crow-1,ccol+1];
                local_mat[3,3:4] = z[crow-1,c(ccol,ccol+1)];
                local_mat[3,5:6] = z[c(crow,crow+1),ccol+1];
                local_mat[3,7]   = z[crow+1,ccol];
                local_mat[3,9:8] = z[c(crow,crow+1),ccol+1];

        } else if ( ccol == ncol(z) ) {
                local_mat[3,2:3] = z[crow-1,c(ccol-1,ccol)];
                local_mat[3,4:6] = z[c(crow-1,crow,crow+1),ccol-1];
                local_mat[3,7]   = z[crow+1,ccol];
                local_mat[3,9:8] = z[c(crow,crow+1),ccol-1];
        }
        reflectn = -1;

        return(list("local_mat"=local_mat,"reflectn"=reflectn));
        }








# Now first shift the local_mat down in z by the amount of the peak, making the peak (0,0,0).
#
# Then for each adjacent column after the first, compute the 
# dot product of each column, the dot product between columns, and then the arccos
# of their ratio, as follows (recall that $ vec \cdot vec = ||vec||$):
#
# theta = arccos[(col(x) . col(x+1) ) / ( (col(x) . col(x) ) \times ( col(x+1) . col(x+1) )]
# acos == arccos
# %*% == dot product: matrix[,2] %*% matrix[,3] == dot product of cols 1 and 3
# * == scalar multiplication
#
# Just for giggles, also compute the area of each triangle, which will be:
# area = 1/2 * ||col(x)|| * ||col(x+1|| * sin(theta)
#
#
# Final results go to a 2x8 matrix, angle_area_mat, with construction of 
# an intermediate matrix, dot_prods:
#
# dot_prods row 1 == col(x) . col(x)      = ||col(x)||
# dot_prods row 2 == col(x+1) . col(x+1)  = ||col(x+1)||
# dot_prods row 3 == col(x) . col(x+1)
# angle_area_mat row 1 == acos(row3/(row1 * row2) = theta
# angle_area_mat row 2 == 0.5 * row1 * row2 * sin(row4)
#
#
# Note one can sum along a row by sum(matrix[row,]).

# Error in shifted[, i] : subscript out of bounds

compute_angles <- function(local_mat) {
        shifted <- shift_z(local_mat);
        dot_prods <- dot_prod_cols(shifted);
        theta <- acos(dot_prods[3,]/(dot_prods[1,]*dot_prods[2,]));
        area <- 0.5 * (dot_prods[1,]*dot_prods[2,]) * sin(theta);
        angle_area_mat <- matrix(rbind(theta,area),nrow=2,ncol=8);

        return(angle_area_mat);
        }








# deparse.level prevents carry-over of row names

shift_z <- function(local_mat) {
        coors <- local_mat[1:2,];
        p <- local_mat[3,1];
        o <- rep(p,10);
        ov <- local_mat[3,] - o;
        shifted <- rbind(coors,ov, deparse.level = 0);

        return(shifted);
        }







# shifted =      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
# dot_prods =            [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#
# checks out using
# shifted <- matrix(1:10,nrow=3,ncol=10,byrow=TRUE)
# but shifted not passed in

dot_prod_cols <- function(shifted) {
        dot_prods = matrix(NA,nrow=3,ncol=8);
        for ( i in 2:9 ) {
                dot_prods[1,i-1] = shifted[,i] %*% shifted[,i];   
                dot_prods[2,i-1] = shifted[,i+1] %*% shifted[,i+1];
                dot_prods[3,i-1] = shifted[,i] %*% shifted[,i+1];
                }

        return(dot_prods);
        }






# now the whole proxy shebang
#
# hard-wired the corner minima, since we know the function; 
# BUT this isn't correct if the peak is in a corner . . . what to do?
# 
# if the peak is at a corner, then when we compute the all_z_ranges and min_z_range,
# one of those values will be 0.  That's fine.
#
# corners are: ul = upper left; ur = upper right; ll = lower left; 
# lr = lower right
#
# that is, c(ul, ur, ll, lr)
#             1   2   3   4
#
#
#
# Modified to improve the shape proxies that are calculated.
#
# Kazic, 10.4.2016


compute_proxies <- function(z,num_proxies,swing) {

        rows <- nrow(z)
	cols <- ncol(z)


# returns [row,col] of maximum value; should be of length 1 (may be spread)

        maxl = which(z == max(z), arr.ind = TRUE)
        minl = which(z == min(z), arr.ind = TRUE)   	


# corners:
#
# ul   ur    ll      lr
# 1,1  1,c   r,1     r,c
#
# slopes:
# ul - ll = left = left water
# ur - lr = right = right water
# ul - ur = top = top n
# ll - lr = bottom = bottom n

        corners <- cbind(z[1,1], z[1,cols], z[rows,1], z[rows,cols])
        rep_corners <- t(corners)
	rownames(rep_corners) <- c("z(ul)","z(ur)","z(ll)","z(lr)")

        e <- extrema(z,maxl,corners)


        displ <- displcmts(z,e["row",1],e["col",1],rows,cols,e["z_max",1],corners)


# aargh!  this is incorrect!  it should be the sum(abs(z) - abs(z_min))!
# see "TODO recalculate the surface volumes <2016-01-04 Mon> :toni:" in
# paraboloid.org
#
# fixed and recalculate!
#
# Kazic, 11.1.2016
#
# no,
#
# surface_vol = sum(abs(z) - abs(e["z_min",1]))
#
# is still not right!  if z_max = 2.879688e+00 and z_min = -4.420917e+01,
# then 3 - 44 < 0!  so want instead,
#
#
# Kazic, 17.1.2016
#
#
# Hmmm, it would be better to calculate a surface_vol relative to an
# absolute measure --- the minimum of Mo17, which is the lowest overall ---
# rather than the z_min of each surface. See
# /Volumes/b/artistry/papers/current/ann_stapleton_ms/first_manuscript/working_draft/stuff_for_figs.org
# at "looking at the category 4 guys in sweep 46" for the thinking here.
#
# Kazic, 10.4.2016
#
#        surface_vol = sum(abs(z - e["z_min",1]))
#
        surface_vol = sum(abs(z - -0.701))

        l <- level_sets(z,as.numeric(e["z_max",1]),swing)


        if ( identical(e["maxn",1],1) ) { c <- find_curvature(maxl,z) }
        else { c <- find_curvature(t(as.matrix(c(e["row",],e["col",]))),z) }


# added a call to curvature for the minimum
#
# Kazic, 20.1.2016


        if ( identical(e["minn",1],1) ) { cm <- find_curvature(minl,z) }
        else { cm <- find_curvature(t(as.matrix(c(minl[1,1],minl[1,2]))),z) } 


        z_proxies <- matrix(NA,num_proxies,1)
        z_proxies <- rbind(e,rep_corners,displ,surface_vol,l,c,cm,deparse.level=1)
	rownames(z_proxies) <- rev(proxies)
	
        return(z_proxies)
        }


























#################### summary statistics ######################




# do a log10 of each value, but preserve zeroes and negative numbers;
# no mixing of two scales (then fcn not 1:1); coerce NAs and NaNs to 0.
#
# Kazic, 27.11.2015

mags <- function(value) {
        if ( identical(value,0) ) { mag <- 0;
	} else if ( identical(value,NaN) ) { mag <- 0;
	} else if ( identical(value,NA) ) { mag <- 0;
	} else if ( identical(value,Inf) ) { mag <- 0;
	} else if ( identical(value,-Inf) ) { mag <- 0;	
	} else if ( value < 0 ) { mag <- -1*log10(abs(value));
	} else { mag <- log10(value); }

        return(mag);
        }







# get the min, max, and range for each row of a matrix

get_ranges <- function(amatrix) {

	ranges <- matrix(NA,nrow(amatrix),3)
	ranges[,1] <- apply(amatrix,1,min)
	ranges[,2] <- apply(amatrix,1,max)
	ranges[,3] <- abs(ranges[,2] - ranges[,1])

        return(ranges)
        }











# compute summary statistics and tack them on to
# make the final matrices
#
# don't take the log of the surface volume in the original results!
#
#
# simplified to remove the transformed results
#
# Kazic, 10.4.2016


# oops, somehow colnames of proxies_matrix isn't getting passed
#
# stopped here

summary_statistics <- function(amatrix,wd) {


        cat("now computing summary statistics\n")

        nulls <- matrix(NA,nrow(amatrix),ncol=1)
        ranges <- get_ranges(amatrix)



# get the magnitudes of the numerical values for better visualization,
# preserving zeroes and signs
#
#	trans <- matrix(mapply(mags,amatrix,SIMPLIFY=TRUE),nrow(amatrix),ncol(amatrix))
#	tranges <- get_ranges(trans)





# write out tabular results as matrices with row names;
# swapped so that results gets the proxies' names from their matrices
#
# Kazic, 29.11.2015

        results <- cbind(amatrix,nulls,ranges)



# if amatrix doesn't have rownames, this won't create them
#
# Kazic, 20.2.2016
#
# now assumes colnames are correctly labelled
#
# Kazic, 10.4.2016
    
        rownames(results) <- rownames(amatrix)
#        cat(paste0(colnames(amatrix),"\n"))
        colnames(results) <- c(colnames(amatrix),"null","min","max","range")
	save(results,file=file.path(wd,"results.RData"))





#	tresults <- cbind(trans,nulls,tranges)
#	tr_rows <- paste("mag",rownames(amatrix))
#	rownames(tresults) <- tr_rows
#       colnames(tresults) <- colnames(results)
#
#	write.table(results,file=file.path(wd,"results.csv"),row.names=TRUE,col.names=TRUE)
#
#	save(tresults,file=file.path(wd,"tresults.RData"))
#	write.table(tresults,file=file.path(wd,"tresults.csv"),row.names=TRUE,col.names=TRUE)
#
#       return(list("results"=results,"tresults"=tresults))
#       return(list("results"=results))
        }












# need to recompute the proxies, since I messed up the surface vol, and
# want to label the columns better in preparation for mushing them all
# together.
#
#
# So, construct the directory names, open the results matrix, read in each
# surface and recompute its values, and at the end write out the results as
# both csv and RData.
#
# Kazic, 11.1.2016
#
# bother, labelling rows "s#:#" has R interpret it as a range when I try
# to use the label to access specific rows.  And, there is something screwy
# with the first entry:
#
# > m[1:10,1:2]
#     s10.1 X0.010234375
# 1   s10:2  0.010713668
# > m["s10:2",1]
# [1] <NA>
# 21366 Levels: s10:10 s10:100 
#
# just to be safe, eliminate the ":"
#
# Kazic, 12.1.2016


# removed production of mushables to speed up calculation, since by now I
# have sweeps of over a million surfaces.  Save the old version, just in case.
#
# Kazic, 10.4.2016





compute_proxies_from_surfaces <- function(sweep_num,swing) {

        wd <- build_wd(sweep_num)
#	file.rename(paste0(wd,"/results.RData"),paste0(wd,"/old.results.RData"))
	z_files <- list.files(wd,"z\\d+\\.RData")


        proxies_matrix = matrix(NA,nrow=length(proxies),ncol=length(z_files))
        surfaces <- regmatches(z_files,regexpr("z[0-9]+",z_files,perl=TRUE))
        colnames(proxies_matrix) <- surfaces


        for ( i in 1:length(z_files) ) {
                load(file.path(wd,z_files[i]))
                z_proxies <- compute_proxies(z,length(proxies),swing)
	        proxies_matrix[,surfaces[i]] <- z_proxies


                if ( i == 1 ) { rownames(proxies_matrix) <- rownames(z_proxies) }
		             # ugh, but I can't figure out a better way
                rm(z)
                }


         colnames(proxies_matrix) <- paste0("s",sweep_num,surfaces)




#        mushable <- t(proxies_matrix)
#
# save the transpose of the proxies_matrix without the added summary
# statistics, for mushing together for faster analysis by catting in the
# shell
#
#	save(mushable,file=file.path(wd,"mushable_proxies.RData"))
#	write.csv(mushable,file=file.path(wd,"mushable_proxies.csv"))





# this writes out the results matrices with the summary data
#
# have removed the return as useless in this context
#
# Kazic, 10.4.2016

        summary_statistics(proxies_matrix,wd)
        }






redo_proxies <- function(sweep_vec,swing) {
        for ( i in 1:length(sweep_vec) ) { compute_proxies_from_surfaces(sweep_vec[i],swing) }
        }










###################### plotting #########################



# now plot and save the surface for later examination
# have to clear3d() to leave just one plotting window, instead of
# thousands; note that *close* do not work!
#
# this is slow!  the problem is that all the graphical
# devices for the rgl plots are opened before closing any of them.  So as
# the number of plots rises, the graphics memory (I guess) becomes exhausted.
#
# ah, he anticipated this problem: see
#      http://www.inside-r.org/packages/cran/rgl/docs/rgl
# for writeWebGL caveats.
#
# just clear the device each time; closing it triggers rgl to open a new
# device at the next iteration
#
# Kazic, 28.11.2015


plot_em <- function(r,q,x,y,z,pal,wd) {

        plot_num <- paste0("z",r,collapse="")


# sweep_num a global defined in sweep.r, which calls this fcn

        plot_title <- paste0("sweep ",sweep_num,", ",plot_num,"  ",paste("(",paste(as.character(t(q)),sep=", ",collapse=", "),")",collapse=""))
	html_file <- paste0(plot_num,".html",collapse="")
	snapshot_file <- paste0(plot_num,"_snapshot.png",collapse="");


        col.ind <- cut(z,1000)


        clear3d()

# usual workaday plotting; have expanded z axis
#
# persp3d(x,y,z,col=pal[col.ind],xlim=c(-45,45),ylim=c(-8,8),zlim=c(-40,0),forceClipregion=TRUE,xlab="water",ylab="nitrogen",zlab="z")
persp3d(x,y,z,col=pal[col.ind],xlim=c(-45,45),ylim=c(-8,8),zlim=c(-60,60),forceClipregion=TRUE,xlab="water",ylab="nitrogen",zlab="z")


# for animation of sweep14: get all the z values on the same scale
#
#        persp3d(x,y,z,col=pal[col.ind],xlim=c(-45,45),ylim=c(-8,8),zlim=c(-60,35),forceClipregion=TRUE,xlab="water",ylab="nitrogen",zlab="z")


# shut off plot title for animations
#
        title3d(plot_title)


# the new version of rgl (0.95.1441) appears to separate writing the
# snapshot from the html files, so a separate call is necessary.  But now
# the file can be written directly.  It does now seem much slower, however.
# Kazic, 15.1.2016

        writeWebGL(dir=wd,filename=file.path(wd,html_file),
                   prefix="",font="Arial",width=700,height=700)
        snapshot3d(file.path(wd,snapshot_file),fmt="png",top=TRUE)
        }







# for separately plotting all surfaces in a directory, as a possible
# work-around for the rgl bug.  Called in plot_dir.r
#
# Kazic, 16.1.2016

# first_num is the first surface to plot, anticipating rgl will crash at
# surface 497; this way I can pick the plotting up from there.


plot_dirs_surfaces <- function(wd,first_num,p,z_files_vec) {

       for ( i in first_num:length(z_files_vec) ) {
                print(paste0("i: ",i))
                load(file.path(wd,z_files_vec[i]))
                q <- p[,i+1]
                plot_em(i,q,x,y,z,pal,wd)
                rm(z)
                }

       }











# heatmap via gplot:
# http://sebastianraschka.com/Articles/heatmaps_in_r.html
#
# heatmap.2(m,Rowv=NULL,Colv=NULL,col=brewer.pal(8,"Spectral"))
# heatmap.2(m,Rowv=NULL,Colv=NULL,col=colorRampPalette(brewer.pal(11,”Spectral”))(40))
# RColorBrewer and color interpolation:
# http://www.r-bloggers.com/r-using-rcolorbrewer-to-colour-your-figures-in-r/
#
# to do: color all of these consistently

plot_results <- function(dir,filename,amatrix,title) {
        png(file.path(paste0(dir,filename,collapse="")),width=10.5*400,height=7.5*400,res=400)
        par(mar=c(5, 5, 3, 2))
	m <- amatrix[nrow(amatrix):1,]
        color2D.matplot(m,c(1,0),c(0,1),c(0,0.85),show.legend=TRUE,nslices=20,axes=FALSE,xlab="parameter combinations",ylab="",main=title,border=NA)
        axis(side=1,at=1:ncol(amatrix),labels=seq(1:ncol(amatrix)))
        axis(side=2,at=1:nrow(amatrix),labels=rev(rownames(amatrix)),las=2)
        dev.off()
        }









# debug(summary_statistics)
