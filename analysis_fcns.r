# this is ../artistry/papers/current/ann_stapleton_ms/our_parts/simulations/analysis_fcns.r

# a subroutine library for analyze_results.r
#
# Kazic, 5.12.2015

rm(list=ls(all=TRUE))

require(SpatialTools)
require(MASS)
require(fields)
require(EnvStats)
require(akima)
require(graphics)
require(geometry)
require(rgl)
require(plotrix)
require(superheat)
require(viridis)



# for executing shape comparisons with multicore parallelism
# change the value of the cores option to increase the number of
# cores used on different machines. See
#
# https://cran.r-project.org/web/packages/doMC/vignettes/gettingstartedMC.pdf

# library(doMC)
# registerDoMC()
# options(cores=1)








############## plotting stuff #######################


# distributions of counts and densities are the same, but the y axes' ranges vary
# for colors, see ../notes/Rcolor.pdf
#
# modified to use with mushed matrices, so names ---> colnames, nrow ---> ncol, and of course
# m[i,1:data_res] ---> m[,i]
#
# Kazic, 26.1.2016

# plot_hists <- function(wd,some_matrix,data_res,names) {

plot_hists <- function(wd,some_matrix) {
	
#        for ( i in 1:nrow(some_matrix) ) {
        for ( i in 1:ncol(some_matrix) ) {
                png(paste0(wd,"/density_",colnames(some_matrix)[i],".png"))
                hist(as.numeric(some_matrix[,i]),freq=FALSE,main = paste("distribution of" ,colnames(some_matrix)[i]),xlab=colnames(some_matrix)[i],ylab="density",col="lightgreen")
                dev.off()
                png(paste0(wd,"/counts_",colnames(some_matrix)[i],".png"))
                hist(as.numeric(some_matrix[,i]),main = paste("counts of" ,colnames(some_matrix)[i]),xlab=colnames(some_matrix)[i],ylab="num",col="hotpink")
                dev.off()
                }
        }









# 2d scatterplot peak position
#
# remember, per [[file:paraboloid.org][point 3 in "more tweaks and checks <2015-11-25 Wed>"]], that
# to see the peaks oriented to match the surfaces, the x axis (nitrogen) is in the forward orientation,
# and the y axis (water) is in the reverse so that the plot matches the z matrix for the surfaces
#
# hints on labelling points are from
# https://chemicalstatistician.wordpress.com/2013/03/02/adding-labels-to-points-in-a-scatter-plot-in-r/
# 
# the optimization goal is to push the peak position towards the upper left quadrant of this plot,
# with as flat a curvature as possible


plot_2ds <- function(wd,results,data_res,combos) {

        png(paste0(wd,"/peak_locatn.png"))
        plot(as.numeric(results["col",1:data_res]),as.numeric(results["row",1:data_res]),main="peak positions",ylim=rev(range(results["row",1:data_res])),xlab="postn nitrogen pk in surface",ylab="postn water pk in surface",pch=19,col="red")
        with(as.data.frame(results), text(results["col",1:data_res],results["row",1:data_res],labels=combos,pos = 4))
        dev.off()
        }








###### dynamic 3d scatterplots
#
# these show distributions, but identifying desirable combinations of
# parameters is still hard; once one has an idea of the ranges, slicing the
# matrices using which is probably better (paraboloid.org, rmk "slicing results: ask for
# low curvature, low volume, don't care about z max").


plot_3ds <- function(wd,results,data_res,combos) {

# 3d plot initialization

        rgl.open(useNULL = rgl.useNULL()) 
        rgl.bg(color="white")
        pp <- dget("std_view.r")
        par3d(pp)


# automate generatn of saving of webgl versions
# writewebgl call fine here


# peak position and height

        plot3d(as.numeric(results["col",1:data_res]),as.numeric(results["row",1:data_res]),as.numeric(results["z_max",1:data_res]),main="peak postn and ht",ylim=rev(range(results["row",1:data_res])),xlab="postn nitrogen peak in surface",ylab="postn water peak in surface",zlab="z_max",pch=19,col="red",size=10)
        with(as.data.frame(results),text3d(results["col",1:data_res],results["row",1:data_res],results["z_max",1:data_res],combos,adj = c(1,1)))



# curvature, volume, and height

        open3d()
        plot3d(as.numeric(results["curvature",1:data_res]),as.numeric(results["surface_vol",1:data_res]),as.numeric(results["z_max",1:data_res]),main="curvature, volume, and ht",xlim=rev(range(results["curvature",1:data_res])),ylim=rev(range(results["surface_vol",1:data_res])),xlab="curvature",ylab="surface volume",zlab="z_max",pch=19,col="blue",size=10)
        with(as.data.frame(results),text3d(results["curvature",1:data_res],results["surface_vol",1:data_res],results["z_max",1:data_res],combos,adj = c(1,1)))



# level sets:  would be nice to flip the orientation of the level 0.50 axis so that all three
# zeroes are in the same corner, but can't get it to work . . .

        open3d()
        plot3d(as.numeric(results["level0.25",1:data_res]),as.numeric(results["level0.50",1:data_res]),as.numeric(results["level0.75",1:data_res]),main="level sets",xlab="level 0.25",ylab="level 0.5",zlab="level 0.75",col="forestgreen",size=10)
        with(as.data.frame(results),text3d(results["level0.25",1:data_res],results["level0.50",1:data_res],results["level0.75",1:data_res],combos,adj = c(1,1)))

        }







# show ranges for proxies passed in as a vector of names
#
# want lower end of "col" and upper end of "row"

show_ranges <- function(vec,results,data_res) {

        for ( i in 1:length(vec) ) {
	        header <- paste0(vec[i],": ")
		rng <- range(as.numeric(results[vec[i],1:data_res]),na.rm=TRUE)
	        print(paste0(header,rng))
		}
        }

















############ labelling and computing assorted summary matrices ######################



# preparation for considering all sweeps together: label each column in
# each summary matrix
#
# matrix_name must be string, that is in double quotes
#
# run on vec_sweep_nums <- c(5,6,8,9,10,11,12,13,14,15,16)
#
# call is:  label_cols_w_sweep_num(c(5,6,8,9,10,11,12,13,14,15,16,18),"","~/me")
#
# fixed problem with rowlabels
#
# Kazic, 12.1.2016
#
# amended to save out the matrix as an R object (.RData).
#
# Kazic, 9.4.2016



label_cols_w_sweep_num <- function(vec_sweep_nums,matrix_name,local) {

        for ( i in 1:length(vec_sweep_nums) ) {
	        header <- paste0("s",vec_sweep_nums[i],"z")
                wd <- build_wd_dbl(local,vec_sweep_nums[i])
                old <- read.csv(paste0(wd,"/",matrix_name,".csv"))

#		new_file <- paste0("new_",matrix_name,".csv")
		new_file <- paste0(matrix_name,".RData")



                new_matrix <- old[,2:ncol(old)]
		colnames(new_matrix) <- paste0(header,seq(1,ncol(old)-1))
		rownames(new_matrix) <- old[,1]
		
#		write.csv(new_matrix,file=file.path(wd,new_file))
                assign(matrix_name,new_matrix)
                save(list=matrix_name,file=file.path(wd,new_file))

                }
        }






# again, matrix_name must be a string: note output file is *mushable*!

invert_matrices <- function(vec_sweep_nums,matrix_name,local,output_dir) {
        root_dir <- "artistry/papers/current/ann_stapleton_ms/our_parts/results/parameter_sweeps"

        for ( i in 1:length(vec_sweep_nums) ) {
                wd <- build_wd_dbl(local,vec_sweep_nums[i])
                old <- read.csv(paste0(wd,"/mushable_",matrix_name,".csv"))
		new_file <- paste0("transposed_mushable_",matrix_name,vec_sweep_nums[i],".csv")
                new_matrix <- t(old)
                write.csv(new_matrix,file=file.path(paste0(local,"/",root_dir,"/",output_dir,"/",new_file)))
                }
        }




# call is: source("sweep_fcns.r");  mush_matrices(c(5,6,8,9,10,11,12,13,14,15,16,18),"proxies","~/me","analysis")

mush_matrices <- function(vec_sweep_nums,matrix_name,local,output_dir) {

	root_dir <- "artistry/papers/current/ann_stapleton_ms/our_parts/results/parameter_sweeps"
        for ( i in 1:length(vec_sweep_nums) ) {
                wd <- build_wd_dbl(local,vec_sweep_nums[i])
                load(file=file.path(paste0(wd,"/mushable_",matrix_name,".RData")))
		print(paste0("sweep num: ",vec_sweep_nums[i],"   dim(mushable):"))
		print(dim(mushable))
		if ( i == 1 ) { m <- mushable } else { m <- rbind(m,mushable) }
                }
	od <- paste0(local,"/",root_dir,"/",output_dir,"/")
	of <- paste0(od,"all_mushed_",matrix_name)
	print("dim(m):")
	print(dim(m))
        write.csv(m,file=paste0(of,".csv"))
        save(m,file=paste0(of,".RData"))
        }





# for the results of find_good_fits.r, print out the ranges of values for
# each equation parameter (a,b,c,d,e), for each allele, for each sweep
# tested.
#
# Kazic, 12.3.2016

# rm(list=ls(all=TRUE))
# setwd("~/me/artistry/papers/current/ann_stapleton_ms/our_parts/results/parameter_sweeps/analysis")
# setwd("/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts/results/parameter_sweeps/analysis")
# output_ranges(c(1457704772,1457705669,1457705750,1457705915,1457706093,1457706180,1457706229,1457706353,1457706529,1457706698,1457706827,1457708175),"sweep_ranges.org")
# output_ranges(c(),"sweep_ranges.org")


# file globbing trick is at
# https://stat.ethz.ch/pipermail/r-help/2008-December/183198.html



output_ranges <- function(timestamp_vector,file) {
        file_handle <- file(file,"a")

        date <- Sys.time()
	cat("* Results of Sweep Analyses for Timestamps",paste0(timestamp_vector),"\n\nRun on",as.character(date),"\n\n",file=file_handle)

        sumry <- matrix(nrow=9,ncol=length(timestamp_vector))
	rownames(sumry) <- c("mo17","b73","x1a","x2b","x3b","x1b","x2a","x3a","opt")
	colnames(sumry) <- as.character(timestamp_vector)

	for ( i in 1:length(timestamp_vector) ) {
	        opt <- ranges(timestamp_vector[i],file_handle)
                sumry[,i] <- opt
		rm(opt)
                }


        cat("* Summary Table for Timestamps\n\n#+tblname: timestamp_summary\n",file=file_handle)
        cat("\n\n|",colnames(sumry),"\n",sep=" | ",file=file_handle)		
        cat("|-----+----+------+--------+----------+----------+-----|\n",file=file_handle)		

        for ( p in 1:nrow(sumry) ) { cat("| ",rownames(sumry)[p]," | ",paste0(sumry[p,],sep=" | "),"\n",file=file_handle) }

        cat("\n\n",file=file_handle)
        close(file_handle)
	}





ranges <- function(timestamp,file_handle) {
        load(list.files(".")[grep(glob2rx(paste0("allele*",timestamp,".RData")),list.files("."))])
        load(list.files(".")[grep(glob2rx(paste0("ps*",timestamp,".RData")),list.files("."))])
        opt <- vector()

        cat("** timestamp ",timestamp,"\n\n",file=file_handle)


        for ( j in 1:nrow(alleles) ) { 
                cat("*** ",rownames(alleles)[j],"\n\n",file=file_handle)
                m <- matrix(nrow=length(rownames(ps))+1,ncol=length(colnames(alleles)))
                rownames(m) <- c("num",rownames(ps))
                colnames(m) <- c(colnames(alleles))
                cat("\n\n| ",colnames(alleles),"\n",sep=" | ",file=file_handle)		
                cat("|-----+----------+---------------------+-------------------+----------+----------+----------|\n",file=file_handle)		

                for ( k in 1:ncol(alleles) ) {
                        f <- unlist(strsplit(alleles[j,k], split=", "))

                        if ( length(f) > 0 ) { 
                                m[1,k] <- length(f)
                                for ( i in 1:nrow(ps) ) { 
		                        m[i+1,k] <- paste0("[ ",min(ps[i,f])," , ",max(ps[i,f]), " ]")
					}
                                } else { m[1,k] <- 0 }
                        }

                for ( n in 1:nrow(m) ) { cat("| ",rownames(m)[n]," | ",paste0(m[n,],sep=" | "),"\n",file=file_handle) }
                cat("\n\n",file=file_handle)

                if ( max(m[1,]) != 0 ) {
                        v <- m[1,]
			opt[j] <- head(sort(v[v!=0]))[1]
                        } else { opt[j] <- 0 }

                rm(m)
                }


        rm(alleles)
        rm(ps)
	return(opt)
        }



















############## shape matching #############################




# whole_matching_shebang("/Volumes/b","/results/parameter_sweeps",c("b73","x1a","x1b","x2a","x2b","x3a","x3b"),0.05,"matching_surfaces_ranges.org")



whole_matching_shebang <- function(local,wd,vector_exptl_surfaces,cutoff,ranges_file) {

        date <- Sys.time()
        timestamp <- as.integer(date)

        load(paste0("..",wd,"/analysis/match_surfs.rda"))


        output_file_handle <- file(paste0("/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts/results/parameter_sweeps/analysis/",ranges_file),"a")
        cat("\n\n\n* Ranges of Values for Shape Matching for Timestamp ",paste0(timestamp)," Run on",as.character(date),file=output_file_handle)


        load_standard_plotting_parameters(local)


        for ( i in 1:length(vector_exptl_surfaces) ) {

                quantles <- quantile(sort(match_surfs[which(match_surfs[,1] == vector_exptl_surfaces[i]),3]),c(0,cutoff,1),na.rm=TRUE,names=FALSE,type=7)
                q_cutoff <- quantles[2]

                mat_indices <- which(match_surfs[,1] == vector_exptl_surfaces[i])
                vector_surfaces <- as.vector(match_surfs[which(match_surfs[mat_indices,3] <= q_cutoff),2])

                cat(paste0("\n\n** Values' Ranges for Simulated Surfaces Matching Shape of ",vector_exptl_surfaces[i],"\n\n"),file=output_file_handle)
                grab_proxies(wd,vector_surfaces,vector_exptl_surfaces[i],timestamp,output_file_handle)

                quickie_plot_surfaces("/Volumes/b",vector_exptl_surfaces[i],vector_surfaces)
                }


        close(output_file_handle)
        }










# given the hit_list of surfaces from find_good_fits.r, grab the list of
# surfaces from all sweeps that match a particular experimental surface,
# put these together, and then compute combined proxy and parameter
# matrices for each target surface using grab_proxies/5.
#
# Kazic, 30.9.2016
  

# ugh!!!! The old version is about the dumbest possible way to do this computation!
# It repeatedly reloaded the same matrices to select out columns.
#
# Two better ways:
#    a. form the set of all surfaces, sort, and get their columns (not by a
#       for loop, but by slicing).
#    b. join all the matrices together into one big matrix, then subset or
#       just save that.
#
# Kazic, 3.10.2016





# combine_matrices_across_sweeps("../results/parameter_sweeps/analysis/surfaces_20_21_22_23_26_27_28_29_34_35_36_38_39_41_42_43_46_47_48_49_53_54_1475192395.RData","/results/parameter_sweeps","phenotype_space_log.org")




combine_matrices_across_sweeps <- function(surfaces_file,output_dir_suffix,output_log_file) {

        date <- Sys.time()
        setwd("/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts/simulations")
        load(surfaces_file)

        almost <- regmatches(surfaces_file,gregexpr("\\_[0-9]{10}\\.RData",surfaces_file,perl=TRUE))[[1]][1]
        timestamp <- regmatches(almost,gregexpr("[0-9]{10}",almost,perl=TRUE))[[1]][1]
        sweeps <- regmatches(surfaces_file,gregexpr("[0-9]{2,10}",surfaces_file,perl=TRUE))[[1]]
        sweeps <- sweeps[-length(sweeps)]

        results_file <- paste0("../results/parameter_sweeps/analysis/honker_",timestamp,"_results.rda")
        p_file <- paste0("../results/parameter_sweeps/analysis/honker_",timestamp,"_p.rda")


        if ( !file.exists(p_file) ) {
	        honker_ps <- make_honker(sweeps,"p",timestamp,p_file)
                honker_results <- make_honker(sweeps,"results",timestamp,results_file)
                }

        combo_surfaces <- combine_surfaces_across_sweeps(hit_list)


        load(results_file)
        subset_results <- honker_results[,combo_surfaces]
        save(subset_results,file=paste0("../results/parameter_sweeps/analysis/subset_",timestamp,"_results.rda"))


        load(p_file)
        subset_ps <- honker_ps[,combo_surfaces]
        save(subset_ps,file=paste0("../results/parameter_sweeps/analysis/subset_",timestamp,"_p.rda"))



        output_file_handle <- file(paste0("/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts/results/parameter_sweeps/analysis/",output_log_file),"a")
        cat("\n\n\n* Ranges of Values of Proxy Searching by Criteria for Timestamp ",paste0(timestamp)," Run on ",as.character(date),"\n\n",file=output_file_handle)

#        grab_proxies(output_dir_suffix,combo_surfaces,"combo_",timestamp,output_file_handle)

        close(output_file_handle)
        }









# load each matrix, concatenate it into a single big honking matrix and save that.
#
# for maximum size of a matrix in R, see
# http://stackoverflow.com/questions/9984283/maximum-size-of-a-matrix-in-r.
# For results matrix, this is 2147483647/34 = 63,161,284 columns; 
# for p matrix, this is 2147483647/5 = 429,496,729 columns.
#
# 2^30 -1 = 1,073,741,823
#
#
# nb: sweep_dims gives the number of columns of values for parameter
# combos, NOT including the summary data of results matrices.  We rewrote
# the p matrices for sweeps 20 and 21 to omit the extraneous first column
# and give the matrices the proper row names.  This simplifies lots of
# subsequent calculations.
#
# Stapleton and Kazic, 19.10.2016



# make_honker(c(20,21,22,23,26,27,28,29,34,35,36,38,39,41,42,43,46,47,48,49,53,54),"p",1475192395,output_file)

# not sweep 19 for p, as it had the extra flipper parameter

# make_honker(c(6,7,8,9,10,11,12,13,14,15,16,17,18,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,53,54,55,56,57,58,59,60,62,63,64),"p",,output_file)


make_honker <- function(sweeps,type,timestamp,output_file) {

        sweep_dims <- as.matrix(read.csv("sweep_dims.csv"))
	sweeps <- sort(as.integer(sweeps))



# the number of honker columns is equal to the sum of the number of steps
# for both results and p matrices.  We just need to read in the right
# sections of each matrix.  For results and all ps, that's just 1:num_steps.

        honker_cols <- sum(sweep_dims[sweeps,2])


        for ( i in 1:length(sweeps) ) {
                print(paste0("starting incorporation of ",type," matrix for sweep ",sweeps[i]," at ",date()))
                sweep_root <- toggle_partition(sweeps[i])
		file <- paste0(sweep_root,"/results/parameter_sweeps/sweep",sweeps[i],"/",type,".RData")
                load(file)
                if ( exists("p") ) { mat <- p; rm(p) } else { mat <- results; rm(results) }
                if ( i == 1 ) {
		        honker <- matrix(nrow=nrow(mat),ncol=honker_cols)
                        rownames(honker) <- rownames(mat)
			cnames <- c()
			cols_so_far <- 1
		        }




# now make this p and results specific

                data_cols <- sweep_dims[which(sweep_dims[,1]==sweeps[i]),2]
                end <- cols_so_far + data_cols - 1
                honker[,cols_so_far:end] <- as.matrix(mat[,1:data_cols])
                cnames <- c(cnames,colnames(mat)[1:data_cols])


		cols_so_far <- end + 1
		rm(mat)
                }


        colnames(honker) <- cnames



# now, if we are doing a results honker, get the corresponding p honker
# get p honker's column names, and use that to reorder the columns of 
# results honker.  Write the reordered results honker out.
#
# hmm, make honker name and substitute p for results???


# buggy!  don't use! couldn't work out references to matrices
#
# Kazic, 13.11.2016

        ## if ( type == "results" ) {
        ##         p_file <- paste0("../results/parameter_sweeps/analysis/honker_",timestamp,"_p.rda")
        ##         load(p_file)
        ##         vars <- ls()
        ##         for ( i in 1:length(vars) ) {
        ##                 test <- regmatches(vars[i],gregexpr("[0-9]{2,10}",vars[i],perl=TRUE))[[1]]
        ##                 if ( length(test) > 0 ) { p_honker <- as.symbol(vars[i]) }
        ##                 }


        ##         ordered_cols <- colnames(p_honker)
        ##         h2 <- honker[,ordered_cols]
        ##         honker <- h2
        ##         rm(h2)
        ##         rm(p_honker)

        ##         }


	
##	honker_name <- regmatches(output_file,gregexpr("honker\\_\\d+\\_[presults]+",output_file,perl=TRUE))[[1]][1]
##        assign(honker_name,honker)
        save(list=honker,file=output_file)
#	return(honker_name)
        }









# simplify make_honker, above, and call from make_honker_denovo.r
#
# Kazic, 3.12.2016

# not sweep 19 for p, as it had the extra flipper parameter

# make_honker_denovo(c(6,7,8,9,10,11,12,13,14,15,16,17,18,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,53,54,55,56,57,58,59,60,62,63,64),"p","phenotype_space_log.org")


# make_honker_denovo(c(9,10),"results","test.org")

make_honker_denovo <- function(sweeps,type,log_file) {

        analysis_root <- toggle_partition(11)
        analysis_dir <- paste0(analysis_root,"/results/parameter_sweeps/analysis/")


        date <- Sys.time()
        start <- as.integer(date)
     

        sweep_dims <- as.matrix(read.csv("sweep_dims.csv"))
	sweeps <- sort(as.integer(sweeps))



# the number of honker columns is equal to the sum of the number of steps
# for both results and p matrices.  We just need to read in the right
# sections of each matrix.  For results and all ps, that's just 1:num_steps.
#
# Kazic, 19.10.2016

        honker_cols <- sum(sweep_dims[sweeps,2])


        for ( i in 1:length(sweeps) ) {
                print(paste0("starting incorporation of ",type," matrix for sweep ",sweeps[i]," at ",date()))
                sweep_root <- toggle_partition(sweeps[i])
		file <- paste0(sweep_root,"/results/parameter_sweeps/sweep",sweeps[i],"/",type,".RData")
                load(file)
                if ( exists("p") ) { mat <- p; rm(p) } else { mat <- results; rm(results) }
                if ( i == 1 ) {
		        honker <- matrix(nrow=nrow(mat),ncol=honker_cols)
                        rownames(honker) <- rownames(mat)
			cnames <- c()
			cols_so_far <- 1
		        }


                data_cols <- sweep_dims[which(sweep_dims[,1]==sweeps[i]),2]
                end <- cols_so_far + data_cols - 1
                honker[,cols_so_far:end] <- as.matrix(mat[,1:data_cols])
                cnames <- c(cnames,colnames(mat)[1:data_cols])


		cols_so_far <- end + 1
		rm(mat)
                }


        colnames(honker) <- cnames



# get ready to save the honker, sorting the columns of the final matrices
# in the same way, to speed access

        sorted <- sort(cnames)
	tmp <- matrix(NA,nrow=nrow(honker),ncol=ncol(honker))
        tmp <- honker[,sorted]


        honker_name <- paste0("honker_",as.character(start),"_",type)
        honker_file <- paste0(analysis_dir,"/",honker_name,".rda")


# this works!

        assign(honker_name,tmp)
        save(list=honker_name,file=honker_file)


        stop <- as.integer(Sys.time())


# write to logfile essentially as per previous entries

        write_honking_to_log(analysis_dir,log_file,date,start,stop,sweeps,honker_file)

        }









write_honking_to_log <- function(analysis_dir,log_file,date,start,stop,sweeps,honker_file) {
        zz <- file(paste0(analysis_dir,log_file), "a")
        elapsed <- stop - start
        cat(paste0("\n\n\n\n* Matrix combination for timestamp ",start," Run on ",as.character(date), "\n\nexecution time (clock): ",elapsed," sec\nfor sweeps: "),file=zz)
        cat(as.character(sweeps),file=zz)
	cat("\ngenerated by analysis_fcns.r:make_honker_denovo/3\n\n",file=zz)
	cat(paste0("Results found in:\n",honker_file,"\n\n"),file=zz)
        close(zz)
        }
















# old version, sets us up for great inefficiencies!

combine_surfaces_across_sweeps_old <- function(hit_list) {
	
	num_entries <- length(hit_list[[1]])
        num_sweeps <- length(hit_list)

        combos <- vector("list",length(num_entries - 1))
	
        for ( i in 2:num_entries  ) {
	        acc <- c()
                for ( j in 1:num_sweeps ) {
                        acc <- c(acc,hit_list[[j]][[i]][[2]])
                        }
			
                exptl_surf <- regmatches(hit_list[[1]][[i]][[1]],gregexpr("[mbxopta0-9]+",hit_list[[1]][[i]][[1]],perl=TRUE))[[1]][1]
                combos[[i-1]] <- list(exptl_surf,acc)	
                }
		
        return(combos)
        }




# new version, just a flattened, sorted list of all simulated surfaces
# think this might be faster than unlisting and picking out the identifiers

combine_surfaces_across_sweeps <- function(hit_list) {
	
	num_entries <- length(hit_list[[1]])
        num_sweeps <- length(hit_list)

        combos <- c()

        for ( i in 2:num_entries  ) {
                for ( j in 1:num_sweeps ) {
                        combos <- c(combos,hit_list[[j]][[i]][[2]])
                        }
                }

        combos <- sort(unique(combos))
        return(combos)
        }








# given a vector of sweeps' surfaces' names, parse each name, grab the
# proxies and parameters, and write these out to summary matrices.
#
# WARNING: number of rows of output matrices hard-wired.
#
# Kazic, 22.3.2016

# wd <-  "~/me/artistry/papers/current/ann_stapleton_ms/our_parts/results/parameter_sweeps";  vector_surfaces <-  c("s30z95090","s30z95860","s31z110065","s31z116145","s32z187012","s32z187663","s33z71005","s33z71010");  stem <- "test"
#
# grab_proxies("/results/parameter_sweeps",c("s30z95090","s30z95860","s31z110065","s31z116145","s32z187012","s32z187663","s33z71005","s33z71010"),"test")




# modified to toggle between athe partitions; so here, wd now is
# essentially just the suffix of the output directory
#
# Stapleton and Kazic, 2.8.2016
#
#
# Further modified to print the summary table of the parameters' values'
# ranges to a file, and improved file naming.
#
# Kazic, 3.8.2016

# grab_proxies("/results/parameter_sweeps",c("s29z729085","s36z368450","s36z387764","s39z1153889","s39z1153901"),"test",111111111,"matching_surfaces_ranges.org")



grab_proxies <- function(wd,vector_surfaces,file_stem,timestamp,output_file_handle) {
        prior_sweep <- regmatches(vector_surfaces[1],gregexpr("[0-9]+",vector_surfaces[1],perl=TRUE))[[1]][1]


# here's the hard-wiring of the matrices' nrows

	proxy_values <- matrix(nrow=34,ncol=length(vector_surfaces))
	parameter_values <- matrix(nrow=5,ncol=length(vector_surfaces))	



# load in the results and parameter matrices for the first surface
# so that we can initialize the output matrices

        sweep_root <- toggle_partition(prior_sweep)
        sweep_dir <- paste0(sweep_root,wd,"/sweep",prior_sweep)

	load(file=file.path(sweep_dir,"results.RData"))
	load(file=file.path(sweep_dir,"p.RData"))			
        if ( !exists("r") ) { r <- results; rm(results) }
        print(paste0("loading files ",file.path(sweep_dir,"results.RData")," and ",file.path(sweep_dir,"p.RData")))



        rownames(proxy_values) <- rownames(r)
        rownames(parameter_values) <- rownames(p)
        colnames(proxy_values) <- vector_surfaces
        colnames(parameter_values) <- vector_surfaces



        for ( i in 1:length(vector_surfaces) ) {

                cur_sweep <- regmatches(vector_surfaces[i],gregexpr("[0-9]+",vector_surfaces[i],perl=TRUE))[[1]][1]

                if ( cur_sweep != prior_sweep ) {

                        rm(r)
                        rm(p)
                        prior_sweep <- cur_sweep

                        sweep_root <- toggle_partition(prior_sweep)
                        sweep_dir <- paste0(sweep_root,wd,"/sweep",prior_sweep)

	                load(file=file.path(sweep_dir,"results.RData"))
	                load(file=file.path(sweep_dir,"p.RData"))			
                        if ( !exists("r") ) { r <- results; rm(results) }
                        }


                proxy_values[,vector_surfaces[i]] <- r[,vector_surfaces[i]]
	        parameter_values[,vector_surfaces[i]] <- p[,vector_surfaces[i]]
                }


        cat(paste0("#+tblname: parameter_ranges_",file_stem,"_",timestamp,"\n"),file=output_file_handle)
        output_value_ranges(parameter_values,output_file_handle)
   
        cat(paste0("\n\n\n#+tblname: proxy_ranges_",file_stem,"_",timestamp,"\n"),file=output_file_handle)
        output_value_ranges(proxy_values,output_file_handle)


        save(proxy_values,file=file.path("/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts/results/parameter_sweeps/analysis",paste0("proxies_",file_stem,"_",timestamp,".rda")))
        save(parameter_values,file=file.path("/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts/results/parameter_sweeps/analysis",paste0("parameters_",file_stem,"_",timestamp,".rda")))

        }









# given a matrix of values for a group of surfaces, compute
# the range of each value
#
# this assumes each type of value is a row.  This is true for the
# parameter values and proxy values.
#
# This version is simplified from output_ranges and is called from grab_proxies.
#
# Kazic, 4.8.2016


output_value_ranges <- function(matrix_values,output_file_handle) {
        output_matrix <- matrix(NA,nrow(matrix_values),ncol=2)

        if ( identical(nrow(matrix_values),5) ) { 
                rownames(output_matrix) <- c("| a","| b","| c","| d","| e") 
                } else { rownames(output_matrix) <- paste0("| ",rownames(matrix_values)) }

        cat("| | min | max |\n|---+-+-|\n",file=output_file_handle)


        for ( i in 1:nrow(matrix_values) ) {
                output_matrix[i,] <- range(matrix_values[i,],finite=TRUE)
                }


        write.table(output_matrix,file=output_file_handle,sep=" | ",row.names=TRUE,col.names=FALSE,eol=" |\n",quote=FALSE)
        }












# load the standard plotting parameters into the global environment
# using the <<- operator.  That trick is in:
# https://www.r-bloggers.com/environments-in-r/
#
# Note that the name of m in plot_surfaces/2
# above has been changed to m_for_z here to preclude namespace collisions
# with the m in get_intersectns/3.  Scoping rules should prevent this, but
# really I shouldn't be so lazy about variable names anyway.
# That m in get_intersectns/3 is more appropriate anyway there, as it's the
# slope.
#
# Kazic, 4.8.2016


# set the palette range to be that of B73
#
# Kazic, 7.2.2017


# modified zlimits to accommodate the range of the surfaces generated by fitting the linear models
#
# Kazic, 21.7.2017

load_standard_plotting_parameters <- function(local) {

        if ( !exists("pp") ) {
                require(rgl)
#	        pp <<- dget(paste0(local,"/simulations/tilted_ann6.r"))


# for compare_surfaces_n_plot

                pp <- dget(paste0(local,"/simulations/std_view.r"))
                }


        x <<- seq(-42,42,0.5)
	y <<- seq(-7.5,7.5,0.5)

# experimental B73 range
#
#        zlimits <<- c(-5,45)
#
# simulated surfaces range

        zlimits <<- c(-250,150)
	
        m_for_z <<- matrix(c(-45,-8,0,45,-8,0,
              -45,-8,0,-45,8,0,
              -45,8,zlimits[1],-45,8,zlimits[2]),nrow=3,byrow=TRUE)
        yticks <<- c(-6,-3,0,3,6)
	
#	zticks <<- c(-5,5,15,25,35,45)
	zticks <<- c(-250,-200,-150,-100,-50,0,50,100,150)


#        jet.colors <<- colorRampPalette( c("violet","blue","green","yellow","red") )
#        load(paste0(local,"/results/parameter_sweeps/ann_replots/zb73.rda"))
#        pal <<- jet.colors(cut(z,1000))
#        rm(z)


# change to the viridis palette

        pal <<- viridis(1000)

        rgl.open(useNULL = rgl.useNULL()) 	
	par3d(pp)
        bg3d(color="white")


#        save_dir <<- "/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts/results/parameter_sweeps/matching_surfaces"
#        save_dir <<- paste0(local,"/artistry/papers/current/ann_stapleton_ms/our_parts/results/parameter_sweeps/matching_surfaces")
#        save_dir <<- paste0(local,"/artistry/papers/current/ann_stapleton_ms/our_parts/results/parameter_sweeps/analysis")

        return(pal)
	
        }








# made this much less dependent on the context of previously declared
# variables in the R process.
#
# Kazic, 12.2.2017
  

quickie_plot_surfaces <- function(local,exptl_surface,vector_surfaces) {

#        save_dir <<- "/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts/results/parameter_sweeps/matching_surfaces"
#        save_dir <<-        "~/me/artistry/papers/current/ann_stapleton_ms/our_parts/results/parameter_sweeps/matching_surfaces"
#
        save_dir <<- paste0(local,"/results/parameter_sweeps/analysis")


        x <<- seq(-42,42,0.5)
	y <<- seq(-7.5,7.5,0.5)
        zlimits <<- c(-5,45)
        m_for_z <<- matrix(c(-45,-8,0,45,-8,0,
              -45,-8,0,-45,8,0,
              -45,8,zlimits[1],-45,8,zlimits[2]),nrow=3,byrow=TRUE)
        yticks <<- c(-6,-3,0,3,6)
	zticks <<- c(-5,5,15,25,35,45)
	
        jet.colors <<- colorRampPalette(c("violet","blue","green","yellow","red") )
        pal <- jet.colors(1000)
	
        pp <- dget(paste0(local,"/simulations/std_view.r"))
        rgl.open(useNULL = rgl.useNULL()) 	
	par3d(pp)
        bg3d(color="white")


        prior_sweep <- regmatches(vector_surfaces[1],gregexpr("[0-9]+",vector_surfaces[1],perl=TRUE))[[1]][1]
        

        for ( i in 1:length(vector_surfaces) ) {
                cur_sweep <- regmatches(vector_surfaces[i],gregexpr("[0-9]+",vector_surfaces[i],perl=TRUE))[[1]][1]
                if ( cur_sweep != prior_sweep ) {
                        prior_sweep <- cur_sweep 
                        }

                sweep_root <- toggle_partition(prior_sweep)
                wd <- paste0(sweep_root,"/results/parameter_sweeps/")

                sweep_dir <- paste0(wd,"sweep",prior_sweep)
		surface <- regmatches(vector_surfaces[i],gregexpr("[0-9]+",vector_surfaces[i],perl=TRUE))[[1]][2]
                load(paste0(sweep_dir,"/z",surface,".RData"))


#                if ( i == 1 ) { col.ind <- cut(z,1000) } # pal <<- jet.colors(cut(z,1000)) }	# 

		
#                persp3d(x,y,z,col=pal[col.ind],xlim=c(-45,45),ylim=c(-8,8),zlim=zlimits,forceClipregion=TRUE,xlab="",ylab="",zlab="",axes=FALSE,box=FALSE,lit=TRUE,specular="black")

                persp3d(x,y,z,col=pal,xlim=c(-45,45),ylim=c(-8,8),zlim=zlimits,forceClipregion=TRUE,xlab="",ylab="",zlab="",axes=FALSE,box=FALSE,lit=TRUE,specular="black")

                segments3d(x=as.vector(t(m_for_z[,c(1,4)])),y=as.vector(t(m_for_z[,c(2,5)])),z=as.vector(t(m_for_z[,c(3,6)])))
                axis3d('x--',pos=c(-45,-8,0)) ; axis3d('y--',pos=c(-45,-8,0),at=yticks) ; axis3d('z-+',pos=c(-45,8,zlimits[1]),at=zticks)
                mtext3d("water",'x--', line = 0.75, at = NULL, pos = NA) ; mtext3d("nitrogen",'y--', line = 0.5, at = NULL, pos = NA) ; mtext3d("z",'z-+', line = 1.25, at = NULL, pos = NA)
		title3d(paste0(vector_surfaces[i]," matching exptl surface ",exptl_surface))



                writeWebGL(file=paste0(save_dir,"/",exptl_surface,"_",vector_surfaces[i],".html"),prefix="",font="Arial",width=700,height=700)
		snapshot3d(file=paste0(save_dir,"/",exptl_surface,"_",vector_surfaces[i],"_snapshot.png"),fmt="png",top=TRUE)
                rm(z)
		clear3d()
                }
		
        }































############ older 3d plotting stuff ###############


# given a vector of surfaces, plot these, being careful of
# rgl's quirks and using the standard view for :ann:'s surfaces
#
# no need for opening a null or scaled 3d window; all this stuff is in the pp
#
#
# need names(vector_surfaces*) if an adorned list is passed in;
# so make sure to strip adornments from the list that's used to call the function
#
# call is plot_surfaces("~/me",c("")) etc.
#
# Kazic, 9.4.2016


# modified to toggle between partitions and to test to be sure we're
# on athe ---- otherwise it sees if it can do it
#
# Stapleton and Kazic, 2.8.2016


plot_surfaces <- function(local,vector_surfaces) {

        try(if ( ( !identical(local,"/Volumes/b") ) & ( length(vector_surfaces) > 450 ) ) stop("too many surfaces!") )


        if ( !exists("pp") ) {
#                require(rgl)
	        pp <- dget(paste0(local,"/artistry/papers/current/ann_stapleton_ms/our_parts/simulations/tilted_ann6.r"))
#	        par3d(pp)
                }


        x=seq(-42,42,0.5) ;  y=seq(-7.5,7.5,0.5)
        zlimits <- c(-5,45)
        m <- matrix(c(-45,-8,0,45,-8,0,
              -45,-8,0,-45,8,0,
              -45,8,zlimits[1],-45,8,zlimits[2]),nrow=3,byrow=TRUE)
        yticks <- c(-6,-3,0,3,6) ; zticks <- c(-5,5,15,25,35,45)
	
        jet.colors <- colorRampPalette( c("violet","blue","green","yellow","red") )
        pal <- jet.colors(1000)


	par3d(pp)
        bg3d(color="white")


        save_dir <- "/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts/results/parameter_sweeps/matching_surfaces"
        prior_sweep <- regmatches(vector_surfaces[1],gregexpr("[0-9]+",vector_surfaces[1],perl=TRUE))[[1]][1]
        

        for ( i in 1:length(vector_surfaces) ) {
                cur_sweep <- regmatches(vector_surfaces[i],gregexpr("[0-9]+",vector_surfaces[i],perl=TRUE))[[1]][1]
                if ( cur_sweep != prior_sweep ) {
                        prior_sweep <- cur_sweep 
                        }

                sweep_root <- toggle_partition(prior_sweep)
                wd <- paste0(sweep_root,"/results/parameter_sweeps/")

                sweep_dir <- paste0(wd,"sweep",prior_sweep)
		surface <- regmatches(vector_surfaces[i],gregexpr("[0-9]+",vector_surfaces[i],perl=TRUE))[[1]][2]
                load(paste0(sweep_dir,"/z",surface,".RData"))


                if ( i == 1 ) { col.ind <- cut(z,1000) }

		
                persp3d(x,y,z,col=pal[col.ind],xlim=c(-45,45),ylim=c(-8,8),zlim=zlimits,forceClipregion=TRUE,xlab="",ylab="",zlab="",axes=FALSE,box=FALSE,lit=TRUE,specular="black")

                segments3d(x=as.vector(t(m[,c(1,4)])),y=as.vector(t(m[,c(2,5)])),z=as.vector(t(m[,c(3,6)])))
                axis3d('x--',pos=c(-45,-8,0)) ; axis3d('y--',pos=c(-45,-8,0),at=yticks) ; axis3d('z-+',pos=c(-45,8,zlimits[1]),at=zticks)
                mtext3d("water",'x--', line = 0.75, at = NULL, pos = NA) ; mtext3d("nitrogen",'y--', line = 0.5, at = NULL, pos = NA) ; mtext3d("z",'z-+', line = 1.25, at = NULL, pos = NA)
		title3d(names(vector_surfaces[i]))


                writeWebGL(file=paste0(save_dir,"/",names(vector_surfaces[i]),".html"),prefix="",font="Arial",width=700,height=700)
		snapshot3d(file=paste0(save_dir,"/",names(vector_surfaces[i]),"_snapshot.png"),fmt="png",top=TRUE)
                rm(z)
		clear3d()
                }
		
        }










##################### shape comparison functions #####################




# standard sets of contours for the experimental surfaces
#
# mo17 set revised 16.1.2017

domed <- c(seq(37,24,-1))
sim_domed <- c(seq(36,24,-1))

# 4 is too close to the peak; 
# 0.275 is too close to the local maximum in the trough
#
# these follow the coloring of the image pretty nicely

mo <-  c(3,2.35,2,1.5,1,0.8,0.6)     












# given a list of contours, find the point on each that intersects the
# "axis of symmetry" of the set of contours:  if they were ellipses, this
# would be the line joining their semi-major axes.
#
# Remove fragments of contours by thresholding the number of points.
# Note threshold for number of points in a contour is hardwired for now.  Instead, ask for all
# contours whose lengths are more than mean(length) - 2*std_dev(length)???
#
# Also, remove closed contours to simplify finding the axis.
#
#
#
# This definition of axis_points is awful and kludgey, but I am stuck.
# Because the ellipses are rotated relative to y, I can't just ask for the minimum
# of the function --- which would lie on the "semi-major axis" if this was a complete ellipse.
#
# Nor is the sign of the derivative helpful: starting from the right end and moving clockwise,
# because these will vary with the contour; and the sign will change in at most two places: once
# when x or y < 0, and once again when one or the other of them flips.  Those points can bound
# the region of interest; but computing the Euclidean distance within that region still doesn't help,
# because the maximum is not where I expect it visually.
#
# Nor can I guess how much to rotate the surface to align the contours parallel to y, and anyway
# I fear calculating a rotation matrix for each contour, for each candidate surface, would be far
# too slow.
#
# So instead, find the range of indices between the minimum values of x and y, take the midpoint, and
# add 1.  This seems to come out pretty close empirically, and damn if it doesn't look right.
# But oh boy, is this ad hockery.

# And indeed it is . . . see below!
#
# Kazic, 5.6.2017



# hey, return all the new_contours, but only use contours whose z >= 24 to find the axis_points:
# otherwise, the points wander leftward in the experimental surfaces.
#
# Kazic, 20.4.2016
#
#
# nah, just use contours >= 24: simpler and faster.
#
# Kazic, 20.4.2016



# have now added a condition for the first if test: if it fails,
# a list of length 0 is returned.  This amounts to saying that if
# a contour is too small, then we return a list of length 0.
#
# Kazic, 7.7.2016


# added another conditional:  if there are no new_contours, then don't
# bother computing the axis_points and levels and return a list of length
# 0. The test is a bit convoluted, but is correct.
#
# Stapleton and Kazic, 14.7.2016




# this is the old diffnt/1 with the kludgey axis finding.  The new one, new_diffnt/3, is below.
#
# Kazic, 5.6.2017

diffnt <- function(contours) {

        new_contours <- list()
	levels <- c()
	
        for ( i in 1:length(contours) ) {


# adjusted the threshold length of each contour downward to accommodate the simulated surfaces
#
# the contours must have enough points and not be a closed cycled
#
# Stapleton and Kazic, 9.11.2016
#
#                if ( ( length(contours[[i]]$x) >= 50 )
#                if ( ( length(contours[[i]]$x) >= 10 )
#
#
# hmmm, in the simulated surfaces, the contours can be cyclic!  so they will fail this
# next test.  So remove this test, at least for now.  It worked fine for the experimental
# surfaces, but those were much better behaved than these.
#
# Kazic and Stapleton, 26.3.2017
#
#
#		   & ( contours[[i]]$x[1] != contours[[i]]$x[length(contours[[i]]$x)] )
#		   & ( contours[[i]]$y[1] != contours[[i]]$x[length(contours[[i]]$y)] ) ) {

                if ( length(contours[[i]]$x) >= 10 ) {
                         if ( length(new_contours) == 0 ) {
			        new_contours <- c(contours[i])
                                levels <- c(contours[[i]]$level)
                                } else {
                                new_contours <- c(new_contours,contours[i])
				levels <- c(levels,contours[[i]]$level)
                                }
                        }
                }



        if ( !identical(as.integer(length(new_contours)),as.integer(0)) ) {


                axis_points <- matrix(nrow=length(new_contours),ncol=2)


                for ( i in 1:length(new_contours) ) {


# The old scheme was pretty kludgey, and didn't really get the central axis
# at the point of maximum curvature of the contours:  see
#
# ../results/parameter_sweeps/ann_replots/z*_intersectns_15.png for results.
#
# Kazic, 5.6.2017

                        t <- which.min(new_contours[[i]]$x)
                        b <- which.min(new_contours[[i]]$y)
                        index <- (t-b)/2 + b + 1

                        axis_points[i,] <- c(new_contours[[i]]$x[index],new_contours[[i]]$y[index])
                        }

                rownames(axis_points) <- levels


# seems like one can't define the list on the fly in the return statement!
# instead, build it separately, and then be sure to use the same variable name
# in the function calling this one.  Not sure about this, but this combination of
# voodoo worked.
#
# Kazic, 18.4.2016

                the_list <- list("axis_points"=axis_points,"new_contours"=new_contours)
                } else { the_list <- list() }


       return(the_list)
       }







# this is the new version with improved axis finding
#
# modified old scheme for choosing axis_points to make this more robust to
# different surfaces, and to handle the x3b case; mo17 is already dealt with
# separately elsewhere.
#
# The old scheme was pretty kludgey, and didn't really get the central axis
# at the point of maximum curvature of the contours:  see
#
# ../results/parameter_sweeps/ann_replots/z*_intersectns_15.png for results.
#
# So instead, try again to maximize the curvature.  It retains the shorter contour
# length and screen for cyclic contours of the original.
#
# Kazic, 5.6.2017

# abandoned for now
#
# Kazic, 12.6.2017



new_diffnt <- function(peak_indices,max_xy,contours) {

        new_contours <- list()
	levels <- c()
	
        for ( i in 1:length(contours) ) {

                if ( length(contours[[i]]$x) >= 10 ) {
                         if ( length(new_contours) == 0 ) {
			        new_contours <- c(contours[i])
                                levels <- c(contours[[i]]$level)
                                } else {
                                new_contours <- c(new_contours,contours[i])
				levels <- c(levels,contours[[i]]$level)
                                }
                        }
                }



        if ( !identical(as.integer(length(new_contours)),as.integer(0)) ) {

                axis_points <- matrix(nrow=length(new_contours),ncol=2)


                for ( i in 1:length(new_contours) ) {


# this MAY NOT BE the curvature and DOES NOT give the semi-major axis, but may be good enough for us

                        ncont <- matrix(c(new_contours[[i]]$x,new_contours[[i]]$y),length(new_contours[[i]]$x),2,byrow=FALSE)
                        curves <- diff(diff(ncont))
			cmax_index <- which.max(curves)
                        cmax <- ncont[cmax_index,]
                        if ( cmax <= max_xy ) {
                                axis_points[i,] <- c(new_contours[[i]]$x[cmax_index],new_contours[[i]]$y[cmax_index])




# not the right approach, but stumped for now:  need something better; maybe combine
# into one condition in selecting the axis point using a which
#
# stopped here; test resulting meshes (obsolete!)

                                } else {
                                        print(paste0("Warning! max curvature above peak for contour ",i))
					axis_points[i,] <- c(NA,NA)
                                        }
                        }

                rownames(axis_points) <- levels


# seems like one can't define the list on the fly in the return statement!
# instead, build it separately, and then be sure to use the same variable name
# in the function calling this one.  Not sure about this, but this combination of
# voodoo worked.
#
# Kazic, 18.4.2016

                the_list <- list("axis_points"=axis_points,"new_contours"=new_contours)
                } else { the_list <- list() }


       return(the_list)
       }
























# given the axis through the peak, find the rays that lie around it, two above and
# two below.
#
#
#
# convenient plotting:
#
# > quartz()
# > png("../results/parameter_sweeps/analysis/contours_b73.png") 
# > image(x,y,z,xlab="water",ylab="nitrogen",main="B73 w/ contours, axis, and original axis_points"); contour(x,y,z,levels=c(seq(40,15,-1)),add=TRUE)
# > points(the_list$axis_points,col="red")
# > 
# > points(max_xy[1],max_xy[2],col="black")
# > dev.off()






# given the peak, the axis, and theta (in degrees!), the offset angle
# between each ray (which gets adjusted below), find the rays relative to the
# peak and return them.
#
# Return a vector where each column is the slope of the ray,
# labelled by ray name.
#
#
#
# Here, we are interested in the slope of each ray:  they will all pass
# through the peak.  So, I don't have to use an affine transformation to
# rotate and translate an arbitrary point:  I just need to use a little
# old-fashioned geometry. 
#
# You might think all one needs to do is change the angles to get different rays.
# But tan is very nonlinear and is undefined in multiple places!
# So instead, just get the value for ray1 relative to
# alpha, and then divide that slope by several constants.
#
# Kazic, 24.4.2016


get_rays <- function(max_xy,slope,new_intercept,thetadeg) {
        theta <- as.numeric(thetadeg)/180 * pi
        x_intercept <- - new_intercept/slope
        alpha <- atan(max_xy[2]/(max_xy[1] - x_intercept))

        ccw <- - pi/2 + alpha
        beta1 <- tan(ccw - (alpha + 2*theta + pi))
        beta2 <- beta1/4
	beta3 <- beta1/12             # splits the arc a bit better than 16 for most surfaces
	beta4 <- beta1/32

        slopes <- c(beta1,beta2,slope,beta3,beta4)
	intercepts <- max_xy[2] - (slopes*max_xy[1])
        rays <- rbind(intercepts,slopes)
	colnames(rays) <- c("r1","r2","r0","r3","r4")
        rownames(rays) <- c("y_intercepts","slopes")

        return(rays)
        }




# from bottom to top, the order is r1, r2, r0, r3, r4
# notice the intercepts are now x-intercepts, not ys.

get_mo17_rays <- function(max_xy,slope,new_intercept,thetadeg) {
        theta <- as.numeric(thetadeg)/180 * pi
        x_intercept <- - new_intercept/slope
        alpha <- atan(max_xy[2]/(max_xy[1] - x_intercept))

        ccw <- - pi/2 + alpha
        beta1 <- tan(ccw - (alpha + 2*theta + pi))
        beta2 <- beta1/2.5
	beta3 <- beta1/8             # splits the arc a bit better than 16 for most surfaces
	beta4 <- beta1/20

        slopes <- c(beta1,beta2,slope,beta3,beta4)
	intercepts <- max_xy[2] - (slopes*max_xy[1])
        rays <- rbind(intercepts,slopes)
	colnames(rays) <- c("r1","r2","r0","r3","r4")
        rownames(rays) <- c("x_intercepts","slopes")

        return(rays)
        }







# returns the points on each ray as a 100 x 6 matrix which columns are
#
#  ly,r1x,r2x,r0x,r3x,r4x
#
# since all share the same ys.
#
# correct!
#
#
# err, not quite:  if we have a slope of 0, then lx will be Inf.  In that
# case, set it to ????? (never did, I guess, 16.1.2017)


make_ray_pts <- function(rays) {

        lx <- matrix(ncol=5,nrow=100)
        ly <- seq(-7.5,6.5,length.out=100)

        lx <- outer(ly,rays[1,],FUN="-")/matrix(rays[2,],nrow=100,ncol=5,byrow=T)
        ray_pts <- cbind(ly,lx)

        return(ray_pts)
        }







# Find the intersection of the axis and each contour:  first find a contour point
# closest to the axis, then interpolate between it and the contour neighbor that straddles
# the axis to find the z of the intersection point.  This is then the mesh point.
#
# Returns a 10 x 1 column vector of the form 
# c(r1.x, r1.y, r2.x, r2.y, . . . ,r4.x,r4.y). 
#
#
# for all pairwise Euclidean distances between the points, see
# http://stackoverflow.com/questions/24746892/how-to-calculate-euclidian-distance-between-two-points-defined-by-matrix-contain
#
#
# I don't think it matters very much if the two contour points straddle the
# ray:  the line between them isn't going to vary that much.  So just take
# the points corresponding to the minimum Euclidean
# distance and point before it, and get the line from them.  From the
# plots, this gets it dead on.
#
#
# Algebra refresh:
#
# y = m1x + b1
# y = m2x + b2
#
# m1x + b1 = m2x + b2
# m1x - m2x = b2 - b1
# x(m1 - m2) = b2 - b1
# x = (b2 - b1)/(m1 - m2)
# y = m1x + b1
#
# 1 ==> rays
# 2 ==> con
#
# x = (conb - raysb)/(raysm - conm)
# y = raysm*x + raysb
#
# correct!





# ha, we realized that the minimum distance between a ray
# and a contour can be the very first point of the contour.  In
# this case, the old calculation of subtracting the point above 
# will fail.  So instead, test to see if that minimal point is
# the first one in the contour and subtract the next contour 
# point.  Otherwise we proceed as before.
#
# old calculation:
#
# just in case the closest_index is the last row, always calculate
# the slope between it and the preceeding row
#
#                m <- (con[closest_index[1],2] - con[closest_index[1]-1,2])/(con[closest_index[1],1] - con[closest_index[1]-1,1])
#
#
#
#
# Stapleton and Kazic, 25.7.2016



get_intersectns <- function(ray_pts,new_contour,rays) {

        con <- cbind(new_contour[[1]]$x,new_contour[[1]]$y)
        interp <- matrix(nrow=nrow(rays),ncol=ncol(rays))
	
        for ( i in 2:ncol(ray_pts) ) {

                dm <- dist2(con,matrix(c(ray_pts[,i],ray_pts[,1]),ncol=2))
                closest_index <- which(dm == min(dm), arr.ind = TRUE)



# if the closest index is in the first row of the contour,
# then take the next point in the contour and interpolate 
# between those two;


                if ( identical(as.numeric(closest_index[1]),1) ) { 

                        m <- (con[closest_index[1],2] - con[closest_index[1]+1,2])/(con[closest_index[1],1] - con[closest_index[1]+1,1])


# otherwise, we are anywhere else in the contour (including the last
# row!), and so interpolate between that point in the contour and the one
# above it.
                        
                        } else {
                                m <- (con[closest_index[1],2] - con[closest_index[1]-1,2])/(con[closest_index[1],1] - con[closest_index[1]-1,1])

                                }





# don't bother dividing by 0!

                if ( identical(con[closest_index[1],1],0) ) {
		        b <- con[closest_index[1],2] } else {
	                b <- (con[closest_index[1],2]) - (m*con[closest_index[1],1])
			}
                interp[,i-1] <- c(b,m)
                }


        x <- (interp[1,] - rays[1,])/(rays[2,] - interp[2,])
        y <- rays[2,]*x + rays[1,]


        con_mesh <- matrix(c(x,y),ncol=2)
        rownames(con_mesh) <- c(paste0("c",new_contour[[1]]$level,colnames(rays)))
        
        return(con_mesh)
        }









# the mo17 version
#
# Algebra refresh again:
#
# (x1,y1) is the closest point
# (x2,y2) is the previous/next point on the contour
# _c -> contour
# _r -> ray
#
# m_c = (y1-y2)/(x1-x2)
# b_c = y1 - m_c * x1
# 
#
# intersect ray and locally linearized contour
#
# y = m_c * x + b_c
# y = m_r * x + b_r
#
# m_c x + b_c = m_r x + b_r
# x = (b_r - b_c)/(m_c - m_r)
# y = m_r * x + b_r
#
# x = (raysb - conb)/(conm - raysm)
# y = raysm*x + raysb


# need to straddle the ray a bit better and fit a tangent line, since
# sometimes the closest point is below the ray
#
# yep, this is much better!
#
# Kazic, 17.1.2017

get_mo17_intersectns <- function(ray_pts,new_contour,rays) {

        con <- cbind(new_contour[[1]]$x,new_contour[[1]]$y)
        interp <- matrix(nrow=nrow(rays),ncol=ncol(rays))


        for ( i in 2:ncol(ray_pts) ) {

                dm <- dist2(con,matrix(c(ray_pts[,1],ray_pts[,i]),ncol=2,byrow=FALSE))


# this closest_index is in the matrix dm! the row value ---
# closest_index[1] --- corresponds to the row of the contour.

                closest_index <- which(dm == min(dm), arr.ind = TRUE)



# if the closest index is in the first row of the contour,
# then take the next point in the contour and interpolate 
# between those two, pulling out the calculation of the adjacent
# point for insurance;


                row <- closest_index[1] + 1
                if ( identical(as.numeric(closest_index[1]),1) ) { 
                        adj <- row + 2
			far <- row + 4


# otherwise, we are anywhere else in the contour (including the last
# row!), and so interpolate between that point in the contour and two
# below it by fitting a tangent to the three contour points.

                        } else { adj <- row - 2; far <- row - 4 }


                tangnt <- lsfit(x=con[c(row,adj,far),1],y=con[c(row,adj,far),2])
                m <- tangnt$coefficients[2]
                b <- tangnt$coefficients[1]		 

                interp[,i-1] <- c(b,m)
                }

        
        x <- (rays[1,] - interp[1,])/(interp[2,] - rays[2,])
        y <- rays[2,]*x + rays[1,]


        con_mesh <- matrix(c(x,y),ncol=2)
        rownames(con_mesh) <- c(paste0("c",new_contour[[1]]$level,colnames(rays)))

        return(con_mesh)
        }


















# for each ray, find the points where it intersects the new_contours by linear
# interpolation; then linearly interpolate the value of the surface at that mesh
# point.  Return a matrix of those interpolated surface values, with columns labelled
# by contour names and rows labelled by concatenate(ray number,coordinate).
# So, ten rows by definition.
#
# hmm, but maybe mesh should be an n x 2 matrix for easier calculation of
# the Euclidean distances?  And then return a list, one element the list of
# contours used and the other the revised mesh?
#
# But checking for rownames is a fast, convenient screen.  So leave as is
# for now, and then reshape the mesh if we decide the Euclidean distance is
# the way to go in comparing meshes.
#
#
# This sequence for ly seems to cover all the experimental surfaces y-ranges,
# and gives pretty good resolution for the intersection.


# for use of lapply:
#
# http://www.r-bloggers.com/r-tutorial-on-the-apply-family-of-functions/
# http://www.r-bloggers.com/efficient-accumulation-in-r/
# https://rollingyours.wordpress.com/category/r-programming-apply-lapply-tapply/
#
# but I'm not making any headway on this!


get_mesh <- function(new_contours,rays) {
        ray_pts <- make_ray_pts(rays)



# the obvious way is slow, but faster than lapply and cbind!

         mesh <- matrix(nrow=0,ncol=2)

         for ( i in 1:length(new_contours) ) {
                 intersectns <- get_intersectns(ray_pts,new_contours[i],rays)
                 mesh <- rbind(mesh,intersectns)
                 }


        return(mesh)
        }











# given a surface and an angle theta in degrees,
# get its mesh for contours at fixed values of z
#
# Note that if a requested contour is above the peak, it won't be
# returned.  However, for B73, contour at z=36 gives two points, one in the
# upper right corner --- which is pretty far away from what we want to
# determine.  So filter out bogus contour fragments: for B73, the fragment
# has very few points.
#
# But then, won't this flexible system cause problems when comparing the
# computed to the experimental surfaces?  Have to check the z-levels and
# compare (x,y) values for contours at the same levels.  Have therefore
# labelled  the rows of axis_points with their levels.
#
# lower bound of contours set to 24:  many surfaces have lower contours,
# but not all,  and some surfaces twist away from the axis as the contours
# descend.
#
#
# Let theta be a variable for now, as I may have to tweak it once I see the
# rays on all the experimental surfaces.
#
# Kazic, 24.4.2016


# return an empty big list if we have any ray with zero slope
#
# Kazic, 6.11.2016

# this is the old version that doesn't pass the surface's peak to diffnt/1
#
# Kazic, 5.6.2017

define_mesh <- function(x,y,z,contour_set,thetadeg) {

        if ( is.vector(contour_set) ) { contours <- contourLines(x,y,z,levels=contour_set) } else { contours <- contourLines(x,y,z,levels=get(contour_set)) }
        the_list <- diffnt(contours)



# need at least three points to define an axis, a pretty loose line!

        if ( nrow(the_list[[1]]) >= 3 ) {

                lsl <- lsfit(x=the_list$axis_points[,1],the_list$axis_points[,2])


# adjust the axis so it lies on the peak; lsl$coefficients[2] is the slope of
# the least squares-fit line

                peak_indices <- arrayInd(which.max(z),c(nrow(z),ncol(z)))
                max_xy <- c(x[peak_indices[1]],y[peak_indices[2]])
                new_intercept <- y[peak_indices[2]] - lsl$coefficients[2]*x[peak_indices[1]]


#                if ( !identical(lsl$coefficients[2],0) ) {
                if ( !identical(lsl$coefficients[2],as.integer(0)) ) {


                        rays <- get_rays(max_xy,lsl$coefficients[2],new_intercept,thetadeg)
                        mesh <- get_mesh(the_list$new_contours,rays)

                        clevels <- c()
                        for ( k in 1:length(the_list$new_contours) ) { clevels <- c(clevels,as.numeric(the_list$new_contours[[k]][1])) }


                        big_list <- list("the_list"=the_list,"max_xy"=max_xy,"slope"=lsl$coefficients[2],"new_intercept"=new_intercept,"rays"=rays,"contours"=clevels,"mesh"=mesh)

                        } else { big_list <- list() }

                } else { big_list <- list() }


       return(big_list)

       }








# same as the old version of define_mesh/4, but pass in more information
# and return just the mesh.
#
#
#
#
# if lsl$coefficients[2] --- the slope --- = 0, then Inf will be returned
# in the subsequent steps as we divide by the slope.  So only call this for
# non-zero slopes.  
#
# Kazic, 1.8.2016


quick_mesh <- function(axis_points,new_contours,x,y,z,thetadeg) {

       lsl <- lsfit(x=axis_points[,1],axis_points[,2])
       peak_indices <- arrayInd(which.max(z),c(nrow(z),ncol(z)))
       max_xy <- c(x[peak_indices[1]],y[peak_indices[2]])
       new_intercept <- y[peak_indices[2]] - lsl$coefficients[2]*x[peak_indices[1]]

       if ( !identical(lsl$coefficients[2],0) ) {
               rays <- get_rays(max_xy,lsl$coefficients[2],new_intercept,thetadeg)
               } else { rays <- get_rays(max_xy,0.0001,new_intercept,thetadeg) }


       mesh <- get_mesh(new_contours,rays)
       return(mesh)

       }







# this is the new version that passes the location of the surface peak to new_diffnt/2
#
# Kazic, 5.6.2017

new_define_mesh <- function(x,y,z,contour_set,thetadeg) {

        if ( is.vector(contour_set) ) { contours <- contourLines(x,y,z,levels=contour_set) } else { contours <- contourLines(x,y,z,levels=get(contour_set)) }


        peak_indices <- arrayInd(which.max(z),c(nrow(z),ncol(z)))
        max_xy <- c(x[peak_indices[1]],y[peak_indices[2]])


        the_list <- new_diffnt(peak_indices,max_xy,contours)



# need at least three points to define an axis, a pretty loose line!

        if ( nrow(the_list[[1]]) >= 3 ) {

                lsl <- lsfit(x=the_list$axis_points[,1],the_list$axis_points[,2])


# adjust the axis so it lies on the peak; lsl$coefficients[2] is the slope of
# the least squares-fit line

                new_intercept <- y[peak_indices[2]] - lsl$coefficients[2]*x[peak_indices[1]]


#                if ( !identical(lsl$coefficients[2],0) ) {
                if ( !identical(lsl$coefficients[2],as.integer(0)) ) {


                        rays <- get_rays(max_xy,lsl$coefficients[2],new_intercept,thetadeg)
                        mesh <- get_mesh(the_list$new_contours,rays)

                        clevels <- c()
                        for ( k in 1:length(the_list$new_contours) ) { clevels <- c(clevels,as.numeric(the_list$new_contours[[k]][1])) }


                        big_list <- list("the_list"=the_list,"max_xy"=max_xy,"slope"=lsl$coefficients[2],"new_intercept"=new_intercept,"rays"=rays,"contours"=clevels,"mesh"=mesh)

                        } else { big_list <- list() }

                } else { big_list <- list() }


       return(big_list)

       }









# stop being cute.  For all but mo17 and x3b, find the peak, take the
# line to the leftmost y axis, and go a fixed angle down from that line to
# produce five more rays
#
# then, find a set of either second derivatives or curvatures along each
# ray, at fixed relative intervals --- e.g., every 2% of the length of the ray



# see https://en.wikipedia.org/wiki/Curvature,
# "local expressions" and also "curvature of a graph"
#
## For a plane curve given parametrically in Cartesian coordinates as γ(t) = (x(t),y(t)), the curvature is
#
##     κ = | x ′ y ″ − y ′ x ″ | ( x ′ 2 + y ′ 2 ) 3 2 , {\displaystyle \kappa ={\frac {|x'y''-y'x''|}{\left(x'^{2}+y'^{2}\right)^{\frac {3}{2}}}},} {\displaystyle \kappa ={\frac {|x'y''-y'x''|}{\left(x'^{2}+y'^{2}\right)^{\frac {3}{2}}}},}
#
## where primes refer to derivatives d⁄dt with respect to the parameter t. The signed curvature k is
#
##     k = x ′ y ″ − y ′ x ″ ( x ′ 2 + y ′ 2 ) 3 2 . {\displaystyle k={\frac {x'y''-y'x''}{\left(x'^{2}+y'^{2}\right)^{\frac {3}{2}}}}.} {\displaystyle k={\frac {x'y''-y'x''}{\left(x'^{2}+y'^{2}\right)^{\frac {3}{2}}}}.}






############# better idea ###############
#
#
# 1.  First, find m fixed ray segments that lie between the peak and the lower
# right quadrant of the evaluation interval, that is where
#
#        -42 <= x <= max_xy[1]
#
# and where
#
#        -7.5 <= y <= max_xy[2]
#
#
# 2.  For those ray segments, find n (50?) "equally spaced" points along the rays.
# By equally spaced, we mean the (x,y) coordinates of n points that have
#
#    1/n * total_Euclidean_length_of_the_ray
#
# Don't think we have to interpolate for these, just advance along the
# line.  Well, it's kinda an interpolation.  These are our "2D mesh points"
# for each ray. 
#
#
# 3.  For each set of 2D mesh points, interpolate the value of z at that
# point. So now we have an n x (m*3) mesh matrix of the form
#
#     rows = n mesh points
#     3 cols = (x,y,z) coordinates of each mesh point
#     m sets of 3 columns are the m fixed rays.
#
#
# 4.  For each interpolated z of the mesh points, find diff(z)/total z to give
# a relative change in z for the Euclidean-equally-spaced mesh points.
# Save that to another matrix to be safe.
#
#
# 5.  Now compute those matrices for all surfaces --- we don't have to scootch
# the simulated surfaces around.  And then, compute the Euclidean distance
# between each experimental and simulated surfaces' diff(z)/total z at each
# point. 
#
#
#
# This might work for Mo17 (choose the upper right quadrant), but what to do
# for x3b??????  Maybe we just say they're incomparable and color
# everything black for that part of the plot.
#
# Stapleton and Kazic, 13.6.2017
  





# hmmm, to compare shapes we need something that in some sense is
# scale-independent
#
## Browse[1]> diff(diff(hi)/diff(x[1:10]))
## [1] 0 0 0 0 0 0 0 0
## Browse[1]> diff(diff(low)/diff(x[1:10]))
## [1] 0 0 0 0 0 0 0 0
## Browse[1]> low; hi
##  [1]  1  2  3  4  5  6  7  8  9 10
##  [1]  100  200  300  400  500  600  700  800  900 1000
##
##
## Browse[1]> low <- c(1,2,3,10,20,30,100,400)
## Browse[1]> hi <- low * 100
## Browse[1]> low
## [1]   1   2   3  10  20  30 100 400
## Browse[1]> hi
## [1]   100   200   300  1000  2000  3000 10000 40000
## Browse[1]> length(diff(low))
## [1] 7
## Browse[1]> diff(diff(low)/diff(x[1:8]))
## [1]   0  12   6   0 120 460
## Browse[1]> diff(diff(hi)/diff(x[1:8]))
## [1]     0  1200   600     0 12000 46000
##
##
## Browse[1]> diff(diff(low)/diff(x[1:8]))
## [1]   0  12   6   0 120 460
## Browse[1]> diff(diff(hi)/diff(x[1:8]))
## [1]     0  1200   600     0 12000 46000
## Browse[1]> 
## [1]   200   200  1400  2000  2000 14000 60000
## Browse[1]> diff(low)/diff(x[1:8])
## [1]   2   2  14  20  20 140 600
## Browse[1]> diff(low)/diff(x[1:8])/diff(hi)/diff(x[1:8])
## [1] 0.04 0.04 0.04 0.04 0.04 0.04 0.04



# https://www.math.upenn.edu/~shiydong/Math501X-6-Geodesics.pdf
# https://mathematica.stackexchange.com/questions/129207/how-to-estimate-geodesics-on-discrete-surfaces
# https://www.ck12.org/book/ck-12-trigonometry---second-edition/section/2.5/






# for lines through peak defined by get_fixed_rays, get
# the values of x, y, and z along those rays.
#
# these have to be by interpolation over the surface z
# and, we want the diff relative to total z along that line segment
#
# get differentials of relative zs along those rays and return a matrix with those
# and their xs and ys, since the latter will vary among the surfaces based on the position
# of the peaks.
#
# return BOTH the matrix of (x,y,z) relative mesh points AND the matrix
# of relative differentials along the rays between successive mesh points.




# pick num_points points along each ray segment lying between the peak and the
# lower left quadrant

# https://math.stackexchange.com/questions/175896/finding-a-point-along-a-line-a-certain-distance-away-from-another-point




# this doesn't work for surfaces simulating mo17 with or without peaks;
# the problem is that the peaks are in the lower right region of the plane,
# rather than in the upper left.  So this would need a special case, and
# then we would still have to figure out HOW to compare two such different surfaces.
#
# Kazic, 19.7.2017


# find_relative_diffs_along_rays <- function(x,y,xs,ys,z,num_points) {

find_relative_diffs_along_rays <- function(type,x,y,z,num_points) {
        peak_indices <- arrayInd(which.max(z),c(nrow(z),ncol(z)))
        max_xy <- c(x[peak_indices[1]],y[peak_indices[2]])
        rays <- get_fixed_rays(type,max_xy)
        ray_edge_termini <- segment_rays(type,rays) 


        ray_lengths <- sqrt((ray_edge_termini[,1] - max_xy[1])^2 + (ray_edge_termini[,2] - max_xy[2])^2)
	ray_length_incr <- matrix(NA,num_points,ncol(rays))
	frac <- 1/num_points
        for ( i in 1:ncol(rays) ) { ray_length_incr[,i] <- seq(ray_lengths[i]*frac,ray_lengths[i],ray_lengths[i]*frac) }


        mesh <- get_relative_xy_ray_mesh(ray_lengths,ray_length_incr,ray_edge_termini,max_xy,num_points,frac)


# now interpolate the value of z at each of the points

#        xyz_mesh <- get_interpolated_zs(xs,ys,z,mesh,num_points)
        xyz_mesh <- get_interpolated_zs(x,y,z,mesh,num_points)


	ztot <- max(z) - min(z)
        zcols <- ncol(xyz_mesh)/3
        z_rels_diff_mat <- matrix(NA,num_points-1,zcols)
	for ( j in 1:zcols) { z_rels_diff_mat[,j] <- diff(xyz_mesh[,j*3]/ztot) }


        return(list(xyz_mesh,z_rels_diff_mat))
        }














# the returned mesh is in the form num_points rows x 2 * num rays, where each pair of
# columns is the (x,y) coordinates of the mesh point:  this is their projection of the
# surface mesh points on the evaluation plane.

get_relative_xy_ray_mesh <- function(ray_lengths,ray_length_incr,ray_edge_termini,max_xy,num_points,frac) {

        mesh <- matrix(NA,num_points,length(ray_lengths)*2)
	counter <- 0

        for ( i in 1:length(ray_lengths) ) {

                x_index <- i + counter
                y_index <- x_index + 1
		
                mesh[,x_index] <-  ray_length_incr[,i]/ray_lengths[i] * ray_edge_termini[i,1] + (frac * max_xy[1])
                if ( identical(i,as.integer(1)) ) { mesh[,y_index] <-  max_xy[2] 
		        } else {
                        mesh[,y_index] <-  ray_length_incr[,i]/ray_lengths[i] * ray_edge_termini[i,2] + (frac * max_xy[2])
			        }
		
                counter <- counter + 1
                }


# an ugly kludge to handle the round-off errors!

        over_pts_indices <- which(mesh[num_points,] < -42)
        if ( length(over_pts_indices) > 0 ) {
	        over_pts <- mesh[num_points,over_pts_indices]
                cat("kludging round-off errors for points ",over_pts,"of the xy mesh\n")
                mesh[num_points,which(mesh[num_points,] < -42)] <- -42
	        }


        return(mesh)
        }








# interpp oscillates z values, which is wierd.  So use bicubic from package akima
# instead.
#
# Stapleton and Kazic, 25.6.2017

get_interpolated_zs <- function(x,y,z,mesh,num_points) {
        cols <- as.integer((ncol(mesh)/2)*3)
	upper <- as.integer(ncol(mesh) - 1)
        all_pts <- matrix(NA,num_points,cols)
	counter <- 0
	
	for ( i in seq(1,upper,2) ) {
#                bar <- interpp(xs,ys,c(z),mesh[,i],mesh[,i+1],linear=TRUE,extrap=FALSE,duplicate="error")
                bar <- bicubic(x,y,z,mesh[,i],mesh[,i+1])
                lower_bd <- i + counter
		upper_bd <- lower_bd + 2
                all_pts[,lower_bd:upper_bd] <- matrix(c(bar$x,bar$y,bar$z),num_points,3,byrow=FALSE)
		counter <- counter + 1
                }
		
        return(all_pts)
        }












# find the matrix of rays' slopes and y intercepts, given the
# max_xy (the peak) and the fixed vector of arbitrary slopes
#
# slopes eyeballed to split the quadrant nicely
#
# https://stackoverflow.com/questions/1571294/line-equation-with-angle
# https://math.stackexchange.com/questions/105770/find-the-slope-of-a-line-given-a-point-and-an-angle
# http://www.mathopenref.com/coordslope.html
# http://www.mathopenref.com/coordintersection.html
#
# > atan(slopes)
# [1] 0.00000000 0.04995840 0.09867846 0.17324567 0.30970294 0.64350111
# > diff(atan(slopes))
# [1] 0.04995840 0.04872006 0.07456721 0.13645728 0.33379816

# get_fixed_rays <- function(max_xy,thetadeg,x,y,z) {
#	image(x,y,z,xlab="water",ylab="nitrogen",main="slopes <- c(0,0.05,0.099,0.175,0.32,0.75)"); points(max_xy[1],max_xy[2]) ; plot_rays(rays,max_xy)



# nb: "trough" type slopes only work for experimental Mo17 surface;
# the simulated surfaces have the peak on the right edge, like the domed
# surfaces
#
# Stapleton and Kazic, 16.7.2017

get_fixed_rays <- function(type,max_xy) {
        if ( type == "domed") { slopes <- c(0,0.05,0.099,0.175,0.32,0.75) }
	else { slopes <- c(0,-0.05,-0.099,-0.175,-0.32,-0.75) }
	
        y_intercepts <- max_xy[2] - slopes * max_xy[1]
        rays <- matrix(c(y_intercepts,slopes),2,6,byrow=TRUE)
	
	colnames(rays) <- c("r0","r1","r2","r3","r4","r5")
        rownames(rays) <- c("y_intercepts","slopes")

        return(rays)
        }

  
















# get the right ray segments: rays intersecting the left and bottom
# edges in the evaluation interval
#
# http://www.mathopenref.com/coordsegment.html


segment_rays <- function(type,rays) {

        if ( type == "domed" ) {

                left_y <- rays[2,] * -42 + rays[1,]                                # y = mx + b
		bottom_x <- ( -7.5 - rays[1,])/rays[2,]                            # (y - b)/m = x
                y_ray_indices <- which(left_y >= -7.5 & left_y <= 7.5)
                x_ray_indices <- which(bottom_x >= -42 & bottom_x <= 42)
                ray_termini <- matrix(NA,6,2)
                ray_termini[y_ray_indices,1] <- -42
                ray_termini[y_ray_indices,2] <- left_y[y_ray_indices]
                ray_termini[x_ray_indices,1] <- bottom_x[x_ray_indices]
		ray_termini[x_ray_indices,2] <- -7.5
                } else {


# again, this only works for the experimental Mo17 surface, where the peak is
# in the upper left corner
#
# Stapleton and Kazic, 16.7.2017

# use the left_y for right_y, but now the meaning is different
#
# we want the value of the segment at the right edge of the plane, so that's where
# x = 42

                        left_y <- rays[2,] * 42 # + rays[1,] 
			bottom_x <- ( -7.5 - rays[1,])/rays[2,]                            # (y - b)/m = x
                        y_ray_indices <- which(left_y >= -7.5 & left_y <= 7.5)
                        x_ray_indices <- which(bottom_x >= -42 & bottom_x <= 42)
                        ray_termini <- matrix(NA,6,2)
                        ray_termini[y_ray_indices,1] <- 42
                        ray_termini[y_ray_indices,2] <- left_y[y_ray_indices]
                        ray_termini[x_ray_indices,1] <- bottom_x[x_ray_indices]
			ray_termini[x_ray_indices,2] <- -7.5

                        }


        return(ray_termini)
        }







# given a matrix of rays' intercepts and slopes and the peak, plot these on the current
# graphics device

plot_rays <- function(rays,max_xy) {

        x_texts <- 3 - (rays[1,]/rays[2,])

        for ( i in 1:ncol(rays) ) {
                abline(rays[1,i],rays[2,i],col="blue")
                text(x_texts[i],0,colnames(rays)[i],col="lightblue")
                }
        }








# plot all the surfaces to see the results; thetadeg is in degrees!
#
# nb:  doesn't work for Mo17!
#
# see
# http://stackoverflow.com/questions/1169456/in-r-what-is-the-difference-between-the-and-notations-for-accessing-the
# for help accessing different levels of the list; but still had to kludge
# for now
#
#
# correct!
#
# Kazic, 30.4.2016


# surfaces <- c("B73","X1A","X1B","X2A","X2B","X3A","X3B"); thetadeg <- 15; i <- 1
#
# make_images(c("B73","X1A","X1B","X2A","X2B","X3A","X3B"),"domed",15)

# this is for exptl surfaces


# this assumes the surfaces are a vector either of experimental surfaces or simulated
# surfaces, but not a mixture.
#
# Stapleton and Kazic, 9.11.2016

make_images <- function(surfaces,contour_set,thetadeg) {

        if ( identical(grep("s[0-9]{2,3}",surfaces[1]),as.integer(1)) ) {
	        sweep <- as.logical(TRUE)
#                wd <- "~/me/artistry/papers/current/ann_stapleton_ms/our_parts/results/parameter_sweeps/"
                wd <- "/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts/results/parameter_sweeps/"
	        analysis_dir <- "analysis/"
                } else {
			sweep <- as.logical(FALSE)
#                        wd <- "~/me/artistry/papers/current/ann_stapleton_ms/our_parts/results/parameter_sweeps/ann_replots/"
                        wd <- "/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts/results/parameter_sweeps/ann_replots/"
	                analysis_dir <- "../ann_replots/"
                        }



        setwd(wd)
	x=seq(-42,42,0.5) ;  y=seq(-7.5,7.5,0.5)

        for ( i in 1:length(surfaces) ) {
	        if ( identical(sweep,TRUE) ) {
		        sweep_num <- regmatches(surfaces[i],gregexpr("[0-9]{2,3}",surfaces[i],perl=TRUE))[[1]][1]
			surf <- regmatches(surfaces[i],gregexpr("[0-9]+",surfaces[i],perl=TRUE))[[1]][2]
		        load(paste0(wd,"sweep",sweep_num,"/z",surf,".RData"))
			surf <- tolower(surf)
                        } else {
                                load(paste0(wd,"z",surfaces[i],".RData"))
				surf <- tolower(surfaces[i]) }



# to create a variable name dynamically and save a data structure to it,
# see
# http://stackoverflow.com/questions/4675755/how-to-save-with-a-particular-variable-name,
# the first answer
#
# Kazic, 17.11.2016

                big_list <- define_mesh(x,y,z,contour_set,thetadeg)
	        if ( identical(sweep,FALSE) ) {
		        newlist <- paste0(surf,"_big_list")
			assign(newlist,big_list)
		        save(list=newlist,file=file.path(analysis_dir,paste0("z",surf,"_mesh.rda")))
			}

                png(paste0(analysis_dir,"z",surf,"_intersectns_",thetadeg,".png"))
                image(x,y,z,xlab="water",ylab="nitrogen",main=paste0(surfaces[i]," with contours, axis, and rays, thetadeg = ",thetadeg))




# arrgh, can't figure out a better way to get the levels, see remark above
# this doesn't work:
#
# big_list$the_list$new_contours[1:length(big_list$the_list$new_contours)][[1]]


#                clevels <- c()
#                for ( j in 1:length(big_list$the_list$new_contours) ) { clevels <- c(clevels,as.numeric(big_list$the_list$new_contours[[j]][1])) }
#		contour(x,y,z,levels=clevels,add=TRUE,col="gray")

		contour(x,y,z,levels=big_list$contours,add=TRUE,col="gray")
		
#                points(big_list$the_list$axis_points,col="red")
                points(big_list$max_xy[1],big_list$max_xy[2],col="black")
#                abline(big_list$new_intercept,big_list$slope,col="black")

                plot_rays(big_list$rays,big_list$max_xy)
                points(big_list$mesh,col="blue")

		dev.off()
                }
         }



quickie_plot_image <- function(x,y,z,mesh_list,out_dir,name_string) {
        png(paste0(out_dir,"z",name_string,".png"))
        image(x,y,z,xlab="water",ylab="nitrogen",main=name_string)
	contour(x,y,z,levels=mesh_list$contours,add=TRUE,col="gray")
        points(mesh_list$max_xy[1],mesh_list$max_xy[2],col="black")
        plot_rays(mesh_list$rays,mesh_list$max_xy)
        points(mesh_list$mesh,col="blue")
	dev.off()
        }




# save each experimental surface's mesh as a list with its vector of
# contour names and matrix of mesh points, named so that all can be loaded
# together (lower-cased line names).

# get_meshes(c("b73","x1a","x1b","x2a","x2b","x3a","x3b"),domed,15)


get_meshes <- function(surfaces,contour_set,thetadeg) {
#        wd <- "~/me/artistry/papers/current/ann_stapleton_ms/our_parts/results/parameter_sweeps/ann_replots/"
        wd <- "/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts/results/parameter_sweeps/ann_replots/"
	setwd(wd)
	analysis_dir <- "../analysis/"
	x=seq(-42,42,0.5) ;  y=seq(-7.5,7.5,0.5)

	contours <- vector("list",length(surfaces))
	meshes <- vector("list",length(surfaces))	


        for ( i in 1:length(surfaces) ) {
                load(paste0(wd,"z",surfaces[i],".rda"))
                big_list <- define_mesh(x,y,z,contour_set,thetadeg)

                for ( j in 1:length(big_list$the_list$new_contours) ) { contours[[i]][j] <- big_list$the_list$new_contours[[j]][[1]] }

                meshes[[i]] <- big_list$mesh
                }


        the_meshes <- vector("list",3)
	the_meshes[[1]] <- surfaces
	the_meshes[[2]] <- contours
	the_meshes[[3]] <- meshes


        save(the_meshes,file=file.path(analysis_dir,"meshes.rda"))
        }








# for each mesh in the meshes, return a list of boxes bounding each mesh
# point, calculated using the box_bound.
#
# So for a point (x,y):
#
#
# 0.9x, 0.9y             0.9x, 1.1y
#
#
# 1.1x, 0.9y             1.1x, 1.1y
#
# the bounds are 0.9 <= x <= 1.1, and 0.9 <= y <= 1.1, if box_bound = 0.1

# outer product
# http://stattrek.com/matrix-algebra/vector-multiplication.aspx?Tutorial=matrix


# > a
#      [,1] [,2]
# [1,]    1    2
# [2,]    3    4
# [3,]    5    6
# > s <- c(0.9,1.1)
# > outer(a,s)
# , , 1

#      [,1] [,2]
# [1,]  0.9  1.8
# [2,]  2.7  3.6
# [3,]  4.5  5.4

# , , 2

#      [,1] [,2]
# [1,]  1.1  2.2
# [2,]  3.3  4.4
# [3,]  5.5  6.6


# > mins <- outer(a,s)[,,1]
# > maxs <- outer(a,s)[,,2]
# > mins
#      [,1] [,2]
# [1,]  0.9  1.8
# [2,]  2.7  3.6
# [3,]  4.5  5.4
# > maxs
#      [,1] [,2]
# [1,]  1.1  2.2
# [2,]  3.3  4.4
# [3,]  5.5  6.6


# > mins <= a && a <= maxs
# [1] TRUE


# # nope

# > boxes <- matrix(nrow=0,ncol=4)
# > boxes[,1:2] <- mins
# > boxes[,3:4] <- maxs
# > boxes
#      [,1] [,2] [,3] [,4]


# # yes

# > boxes <- matrix(nrow=3,ncol=4)
# > boxes[,1:2] <- mins
# > boxes[,3:4] <- maxs
# > boxes
#      [,1] [,2] [,3] [,4]
# [1,]  0.9  1.8  1.1  2.2
# [2,]  2.7  3.6  3.3  4.4
# [3,]  4.5  5.4  5.5  6.6





bound_boxes <- function(box_bound,meshes) {
        bounds <- c(1 - box_bound, 1 + box_bound)
        boxels <- list()


        for ( i in 1:length(meshes) ) {
        	boxes <- matrix(nrow=nrow(meshes[[i]]),ncol=4)
                box <- meshes[[i]] %o% bounds
                boxes[,1:2] <- box[,,1]
                boxes[,3:4] <- box[,,2]		

                boxels[[i]] <- boxes
                }


        return(boxels)
        }















box_match <- function(boxes,mesh) {
        if ( identical(TRUE,boxes[,1:2] <= mesh && mesh <= boxes[,3:4]) ) { result <- TRUE } else { result <- NULL }

        return(result)
        }









# eventually, return NULL if test fails, once we know the bounds for the
# various criteria
#
# for boxes, failure is easy (FALSE); but frobenius and eucd both return a number,
# so what is failure for these?  Will know once I compare all the experimental
# surfaces to each other:  surfaces outside each category help set an upper bound.


compare_shapes_aux_aux <- function(criterion,mesh,exptl_mesh,boxes) {


# Frobenius norm, the square of RMSD

        if ( criterion == "fro" ) {
                net <- exptl_mesh - mesh
                result <- norm(net,type="F")


# Euclidean distance:  only interested in the distances along the diagonal,
# as these compare corresponding mesh points
#
# see compute_surfaces/10 for a surely faster approach than dist2/diag

                } else if ( criterion == "eucd" ) {
                        dm <- dist2(exptl_mesh,mesh)
                        result <- max(diag(dm))


# boxes

                } else if ( criterion == "boxes" ) {
		        result <- box_match(boxes,mesh)


                } else { result <- NULL }

        return(result)
        }






# added some conditionals to make the output more informative.  The key
# thing is to notice if some surfaces that were submitted failed the
# initial screens and didn't have their shapes examined.  num_missing is
# the count of these.
#
# Stapleton and Kazic, 14.7.2016





# the structure of the input list of returned_matches is confusing!
#
# returned_matches is a list of lists.  Each element at the top
# level is indexed by i:  this returns a list, not a singular value.
#
# The second, nested list, indexed by j, is a list of two elements.
#
#        The first, [[i]][[1]][[1]], is the name of the simulated surface that
#        matched something.  Here, j = 1.
#
#        The second and subsequent elements are the vectors of match information
#        for each experimental surface.  So [[i]][[2]] is the vector, and its elements
#        are:
#
#                [[i]][[j >= 2]][1], the experimental surface matched;
#                [[i]][[j >= 2]][2], the numerical value of the match; and
#                [[i]][[j >= 2]][3], the criterion used.
#
# Kazic, 29.7.2016


  
write_matches <- function(local,file,start,stop,surfaces_timestamp,criterion,box_bound,surfaces,returned_matches,cutoff) {

       shape_timestamp <- as.integer(stop)
       elapsed <- shape_timestamp  - as.integer(start)

       wd <- "/results/parameter_sweeps/"
       zz <- file(file.path(paste0(local,wd,"analysis"),file), "a")
       cat(paste0("\n\n** shape comparison timestamp ",shape_timestamp,"\n\nstarted on date: ",start, "\nfinished on date: ",stop,"\nexecution time (clock): ",elapsed,"\nfor simulated surfaces found at timestamp: ",surfaces_timestamp,"\nsurfaces: "),file=zz)
       cat(surfaces,file=zz)
       cat(paste0("\ncriterion: ",criterion,"\nbox bounds: " ,box_bound,"\n\n"),file=zz)

       num_missing <- returned_matches[[1]]
       matches <- returned_matches[[2]]




       if ( !identical(as.integer(num_missing),as.integer(0)) ) { cat(paste0("Noto bene: ",num_missing,"surfaces were not checked for good reasons!\n\n"),file=zz) }



       if ( !identical(as.integer(length(matches)),as.integer(0)) ) {


                cat("\n\n*** matched surfaces\n\n",file=zz)
             
                cat(paste0("#+tblname: shape_matching_",surfaces_timestamp,"\n"),file=zz)
                cat("| exptl surface | sim surface | distance | criterion |\n|-+-+-+-|\n",file=zz)


	
                for ( i in 1:length(matches) ) {
	 	
                        sim_surface <- matches[[i]][[1]][1]
	 	        for ( j in 2:length(matches[[i]]) ) {
                                cat(paste0(" | ",matches[[i]][[j]][1]," | ",sim_surface," | ",matches[[i]][[j]][2]," | ",matches[[i]][[j]][3]," |\n"),file=zz)
                                }
                        }



# write out the ranges for the parameters for the matching surfaces
#
# the trick is to flatten the list into a list of three elements, each
# of which is a vector.  Messing with dataframes introduces too many
# unreliabilities:  in particular, it is hard to keep numeric vectors
# numeric! 
#
# Kazic, 27.8.2016


                cat("\n\n\n*** Ranges of Values for Matched Shapes\n\n",file=zz)
                remats <- restructure_match_list(matches)


                for ( i in 1:length(surfaces) ) {
	

# get distances for a particular exptl surface

                        good_dists <- sort(as.numeric(remats[[2]][which(remats[[1]] == surfaces[i])]))


                        if ( length(good_dists) > 0 ) {
                                cat(paste0("\n\n**** Values' Ranges for Simulated Surfaces Matching Shape of ",surfaces[i],"\n\n"),file=zz)


# can't construct a quantile on only one value, and can't get a range for that,
# either


                                if ( length(good_dists) > 1 ) {
				
                                        quantles <- quantile(good_dists,c(0,as.numeric(cutoff),1),na.rm=TRUE,names=FALSE,type=7)
                                        q_cutoff <- quantles[2]

                                        vector_surfaces <- remats[[3]][which(as.numeric(remats[[2]]) <= q_cutoff)]
                                        grab_proxies(wd,vector_surfaces,as.character(surfaces[i]),surfaces_timestamp,zz)
				
                                        } else {
                                                 vector_surfaces <- remats[[3]][which(remats[[1]] == surfaces[i])]
                                                 grab_proxies(wd,vector_surfaces,as.character(surfaces[i]),surfaces_timestamp,zz)
                                                 }
				}
                        }	        

               } else { cat("NO MATCHES FOUND\n\n\n",file=zz) }

       close(zz)
       }
























# well, why not just look at all the surfaces in a directory?  So let's construct the simulated vector,
# and then pass that in to compare_surfaces/8
#
# local <- directory where meshes live, i.e., /Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts
# sweep_root <- partition where surfaces live, so either:
#
#               /Volumes/a/toni/ann_sweeps   OR
#               /Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts
#
# sweep_num <- the particular sweep number we want to search


# may be incomplete; we're thinking we don't want to do this, anyway
#
# Stapleton and Kazic, 14.7.2016



# compare_surfaces_in_dir("/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts","/Volumes/a/toni/ann_sweeps",36,domed,15,"fro",0.1,"shapes.org")
# compare_surfaces_in_dir("/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts","/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts",36,domed,15,"fro",0.1,"shapes.org")

# compare_surfaces_in_dir <- function(local,sweep_root,sweep_num,contour_set,theta,criterion,box_bound,file) {
#
#        setwd(paste0(sweep_root,"/results/parameter_sweeps/sweep",sweep_num))
#        simulated <- list.files(pattern='z.*\\.RData',recursive=TRUE)
#
#
#
#        start <- as.integer(as.POSIXct(Sys.time()))
#        now <- date()
#       
#        load(paste0(local,"/results/parameter_sweeps/analysis/meshes.rda"))
#
#        surfaces <- the_meshes[[1]]
#        contours <- the_meshes[[2]]
#        meshes   <- the_meshes[[3]]
#        if ( identical(criterion,"boxes") ) { boxes <- bound_boxes(box_bound,meshes) } else { boxes <- list() }
#
#        wd <- paste0(sweep_root,"/results/parameter_sweeps/")
#        setwd(wd)
#        x=seq(-42,42,0.5) ;  y=seq(-7.5,7.5,0.5)
#
#        prior_sweep <- sweep_num
#
#
#        returned_matches <- match_surfaces(simulated,wd,prior_sweep,surfaces,contours,meshes,x,y,contour_set,theta,criterion,boxes)
#
#        stop <- as.integer(as.POSIXct(Sys.time()))
#
#        write_matches(wd,file,start,now,stop,criterion,box_bound,surfaces,returned_matches)
#        }























# for a list of surfaces (previously identified through find_good_fits.r
# and its criteria), compare the surfaces' intrinsic shapes to those of the
# experimental surfaces.  The first comparison is just the list of contours
# produced by the automatic criteria; the second is matching the two
# meshes. 
#
# For the second comparison, there are two possibilities.  One is to
# compare the meshes using a pairwise Euclidean distance matrix (eucd), and the
# form of the mesh data has been set to easily allow this.  The second is
# to compare the locations of the points within some arbitrary margin, such
# as 10% (boxes).
#
# The data for the experimental meshes are stored as a list.  The first
# element is the vector of the names of the surfaces; the second element is
# a list of the contours used for each surface; and the third is the list
# of mesh points and their coordinates for each surface.


# precompute and save the boxes for the target surfaces' meshes 
#
# Kazic, 7.7.2016


# compare_shapes("~/me/artistry/papers/current/ann_stapleton_ms/our_parts",domed,15,c("b73","x1a","x1b","x2a","x2b","x3a","x3b"),"fro",0.1,"shapes.org")

# compare_shapes("/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts",domed,15,c(),"fro",0.1,"shapes.org")

# compare_shapes("/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts","/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts",domed,15,c("s36z239508.RData", "s36z258813.RData", "s36z276480.RData", "s36z278118.RData", "s36z294282.RData", "s36z295667.RData", "s36z314972.RData", "s36z315099.RData", "s36z332521.RData", "s36z334277.RData", "s36z334278.RData", "s36z334286.RData", "s36z334404.RData", "s36z351953.RData", "s36z353347.RData", "s36z353465.RData", "s36z353591.RData", "s36z353709.RData", "s36z353718.RData", "s36z371140.RData", "s36z371385.RData", "s36z372534.RData", "s36z372652.RData", "s36z372770.RData", "s36z372896.RData", "s36z372897.RData", "s36z372905.RData", "s36z373023.RData", "s36z373032.RData", "s36z390318.RData", "s36z390319.RData", "s36z390327.RData", "s36z390445.RData", "s36z390572.RData", "s36z390699.RData", "s36z390825.RData", "s36z390834.RData", "s36z392201.RData", "s36z392328.RData"),"fro",0.1,"shapes.org")


# compare_shapes("/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts","/Volumes/a/toni/ann_sweeps",domed,15,c(),"fro",0.1,"shapes.org")

# compare_shapes("/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts","/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts",domed,15,c("s36z239508.RData","s22z2.RData"),"fro",0.1,"shapes.org")

# compare_shapes("/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts","/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts",domed,15,c("s22z2.RData"),"fro",0.1,"shapes.org")



# modified to be called from match_sim_surfaces.r, which loads the surfaces
# stored in the file for a given timestamp from find_good_fits.r.  Here,
# compare_shapes/7 finds the correct partition for each surface in
# match_surfaces/11 and returns a list of the number of surfaces that
# failed the initial screening and the matching simulated surfaces.
# write_matches/9 writes these results to the file.
#
# Kazic, 15.7.2016


# we again have the apples to oranges problem:  see
# compare_shapes_by_z_dists/3 below.
#
# Kazic, 20.11.2016



# newish

compare_shapes <- function(local,contour_set,thetadeg,cutoff,surfaces_timestamp,hit_list,criterion,box_bound,file) {

       start <- as.POSIXct(Sys.time())

       
       load(paste0(local,"/results/parameter_sweeps/analysis/meshes.rda"))

       surfaces <- the_meshes[[1]]
       contours <- the_meshes[[2]]
       meshes   <- the_meshes[[3]]
       if ( identical(criterion,"boxes") ) { boxes <- bound_boxes(box_bound,meshes) } else { boxes <- list() }

       x=seq(-42,42,0.5) ;  y=seq(-7.5,7.5,0.5)



# need to convert hit_list to simulated: want to concatenate the vectors of
# actual surfaces together


       simulated <- c()


       for ( i in 1:length(hit_list) ) {
               for ( j in 2:length(hit_list[[i]]) ) { simulated <- c(simulated,hit_list[[i]][[j]][[2]]) }
              }








# toggled off as we are now doing this for real
#
# Kazic, 24.6.2016
#
#       prior_sweep <- "ann_replots"
#
#
# toggled on

       prior_sweep <- regmatches(simulated[1],gregexpr("[0-9]+",simulated[1],perl=TRUE))[[1]][1]
       



       returned_matches <- match_surfaces(simulated,prior_sweep,surfaces,contours,meshes,x,y,contour_set,thetadeg,criterion,boxes)

       stop <- as.POSIXct(Sys.time())

       write_matches(local,file,start,stop,surfaces_timestamp,criterion,box_bound,surfaces,returned_matches,cutoff)
       }










# this is simpler than a look-up table, but presumes all new sweeps will be
# in /Volumes/a/toni/ann_sweeps

# note that the two conditions are not disjoint; the first permits
# switching between laielohelohe and athe.
#
# Kazic, 20.11.2016

# modified threshold to reflect copying of sweeps from /athe/b/.... to
# /athe/a/...
#
# Kazic, 2.12.2016

toggle_partition <- function(prior_sweep) {

        machine <- Sys.info()["nodename"]

        if ( machine == "athe.rnet.missouri.edu" ) {
                if ( prior_sweep <= 10 ) { sweep_root <- "/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts"
	                } else { sweep_root <- "/Volumes/a/toni/ann_sweeps" }



                } else { sweep_root <- "/Users/toni/me/artistry/papers/current/ann_stapleton_ms/our_parts" }


        return(sweep_root)
        }









# idea is to write a vector of indices for filtered
# candidate exptl surfaces (based on matching contour vectors), and use
# that vector to access the meshes for comparison.  criteria =
# {boxes,eucd,fro} and these are toggled.


# see
# http://stackoverflow.com/questions/2436688/append-an-object-to-a-list-in-r-in-amortized-constant-time-o1
# for efficiently growing flat lists.
#
# My approach is to check the length of the prior
# version and then add the next element.


# output list deeply nested, otherwise correct!

# for help with exception handling, see:
# http://www.r-bloggers.com/error-handling-in-r/


match_surfaces <- function(simulated,prior_sweep,surfaces,contours,meshes,x,y,contour_set,thetadeg,criterion,boxes) {


        matches <- list()
        num_tested <- 0


        for ( i in 1:length(simulated) ) {


# toggled off as we are now doing this for real
#
# Kazic, 24.6.2016
#
#                sweep_dir <- "ann_replots"
#                load(paste0(wd,sweep_dir,"/z",simulated[i],".rda"))
#                test <- simulated[i]





                cur_sweep <- regmatches(simulated[i],gregexpr("[0-9]+",simulated[i],perl=TRUE))[[1]][1]
                if ( cur_sweep != prior_sweep ) { prior_sweep <- cur_sweep }



                sweep_root <- toggle_partition(prior_sweep)


                test <- regmatches(simulated[i],gregexpr("[0-9]+",simulated[i],perl=TRUE))[[1]][2]
                load(paste0(sweep_root,"/results/parameter_sweeps/sweep",prior_sweep,"/z",test,".RData"))


# and print the name of the simulated surface to the screen
#

                surface_name <- paste0("s",prior_sweep,"z",test)
                cat(surface_name,"\n")


# for the simulated surface, get its tidied up contours, and construct the list
# of the names of those tidied contours (names == z-values of the contours)
# 
# test_contours may be empty!  in that case, just increment the count of
# num_tested and go on.
#
# Use get() to get the vector of values stored as a variable mapped to
# contour_set:
#
# http://stackoverflow.com/questions/10429829/how-to-get-value-when-a-variable-name-is-passed-as-a-string

                test_contours <- contourLines(x,y,z,levels=get(contour_set))


                num_tested <- num_tested + 1
                if ( length(test_contours) >= 1 ) {


                        the_list <- diffnt(test_contours)




# if the_list returned from diffnt/1 is empty, then skip all this; in
# either case, clear all the intermediate variables for a particular surface


                        if ( !identical(as.integer(length(the_list)),as.integer(0)) ) {

                                pruned_contours <- c(length=0)

                                for ( j in 1:length(the_list$new_contours) ) { pruned_contours[j] <- the_list$new_contours[[j]][[1]] }



# test the pruned_contours against the contours for each experimental
# surface (1:length(contours[[k]])), recording in a vector match_vec the
# index number of all experimental surfaces that match.  Then pass
# match_vec, z of the simulated surface, and the meshes to the
# mesh-matching algorithm, compare_shapes_aux/7.  Calling compare_shapes_aux here
# is better because *.RData doesn't have to be loaded twice.
#
# Otherwise, if match vector remains NULL (empty), then go to the next simulated surface. 


                                match_vec <- c()

                                for ( k in 1:length(surfaces) ) {
                                        if ( identical(as.numeric(pruned_contours),contours[[k]]) ) {
                                                len <- length(match_vec)
				                match_vec[len+1] <- k
				                }
                                        }



# declarative clarity trumps efficiency, with grins


                                if ( !identical(match_vec,NULL) ) {
                                        if ( length(levels(factor(the_list$axis_points[,2]))) > 1 ) {
                                               mesh <- quick_mesh(the_list$axis_points,the_list$new_contours,x,y,z,thetadeg)
                                               comprsns <- compare_shapes_aux(match_vec,mesh,meshes,surfaces,criterion,boxes)
                                               } else {
                                                       mesh <- c()
                                                       comprsns <- c()
                                                       print(cat("mesh generation fails for surface ",surface_name,"\n"))
                                                       }


# some more exception handling has been embedded in quick_mesh/6,
# rather than using the tryCatch here.
#
# Kazic, 1.8.2016
#
#
#                                        tryCatch(mesh <- quick_mesh(the_list$axis_points,the_list$new_contours,x,y,z,thetadeg),error = function(e) { print(paste0("mesh fails for ",surface_name,"\n"))})

#                                        comprsns <- compare_shapes_aux(match_vec,mesh,meshes,surfaces,criterion,boxes)
                                        } else {
                                        comprsns <- c()
                                        mesh <- c()
                                        
                                        }





# more efficient list expsn?

                                if ( !identical(comprsns,NULL) ) {
                                        len <- length(matches)
                                        matchee <- regmatches(simulated[i],gregexpr(".RData",simulated[i],perl=TRUE),invert=TRUE)[[1]][1]
                                        matches[[len+1]] <- c(matchee,comprsns)
                                        rm("matchee")
			                }

                                rm(list=c("pruned_contours","match_vec","mesh","comprsns"))

                                } 
                        
                                rm("the_list")
                        }

                rm(list=c("z","test_contours"))
                }



        num_missing <- length(simulated) - num_tested
        returned_matches <- list(num_missing,matches)

        return(returned_matches)
        }














# for a global measure, compute rmsd over the euclidean distance matrix
# check casp and dietlind for newer statistics that might be better than
# rmsd.  Like rmsd a lot better than just silly summing over absolute
# values.   Read zemla2003 and perez2012 to be sure, but we chickens prefer
# rmsd. 

# see rustamov2013 for map-based transformation fcns for shape comparison
# http://dl.acm.org/citation.cfm?doid=2461912.2461959
#
# http://math.stackexchange.com/questions/507742/distance-similarity-between-two-matrices
#
# http://www.cis.upenn.edu/~cis515/
# http://www.cis.upenn.edu/~jean/home.html
# http://www.cis.upenn.edu/~cis515/cis515-notes-15.html
# https://inst.eecs.berkeley.edu/~ee127a/book/login/l_mats_norms.html
#
# rmsd:
# https://en.wikipedia.org/wiki/Root-mean-square_deviation_of_atomic_positions
#
# rmsd is the sqrt of frobenius norm
#
# R: https://stat.ethz.ch/R-manual/R-devel/library/base/html/norm.html
# R: http://finzi.psych.upenn.edu/library/bio3d/html/rmsd.html
#
# can a more efficient list append operation be chosen?


compare_shapes_aux <- function(match_vec,mesh,meshes,surfaces,criterion,boxes) {


# initialize a list of length 0

        results <- list()


        for ( i in 1:length(match_vec) ) {
	        matchee <- match_vec[i]
                result <- compare_shapes_aux_aux(criterion,mesh,meshes[[matchee]],boxes[[matchee]])


                if ( !identical(result,NULL) ) {
		        len <- length(results)
                        results[[len+1]] <- c(surfaces[matchee],result,criterion)
                        }
                }


# if the length of results remains equal to 0, then we have no
# results, so set results to NULL.  Otherwise, return what we got.

        if ( identical(as.integer(length(results)),as.integer(0)) ) { results <- NULL }

        return(results)
        }











# we want the ranges of parameter values when we're getting the matched
# shapes directly, but it's easiest to restructure the list of returned
# matches into something simpler, I hope.  Right now that looks like a data 
# frame.
#
# The first column is the names of the matched experimental surfaces;
# the second column is the distance, by whatever criterion, between the
# simulated surface and the matching surface --- BUT AS STRINGS!!!; 
# and the third column is the name of the simulated surface.
#
# Stapleton and Kazic, 26.8.2016



restructure_match_list <- function(matches) {


        exptls <- vector(length=0)
        distances <- vector(length=0)
        sim_surfs <- vector(length=0)



        for ( i in 1:length(matches) ) {
                for ( j in 2:length(matches[[i]]) ) {
                        exptls[length(exptls)+1] <- matches[[i]][[j]][1]
                        distances[length(distances)+1] <- matches[[i]][[j]][2]
                        sim_surfs[length(sim_surfs)+1] <- matches[[i]][[1]]

                        }
                }

        restructured_list <- list(exptls,distances,sim_surfs)
        return(restructured_list)

        }




######################### slope analysis ##################
#
# We want to test the hypothesis that the explored parameter space missed
# the experimental surfaces because we over-emphasized water compared to
# nitrogen, producing the set of right hunch-backed surfaces.  To do this,
# we want to explore a large parameter range, but we really don't want to
# look at bazillions of contour plots.  So instead, for each simulated
# surface, we will compute the slope and intercept of the least-squared fit
# of the central defining axis.  If we're right, then we should see the
# slope rotate around the peak.  We're taking the peak and axis for each
# surface individually, because that's what we want to see, rather than
# using the peak and slope of each of the experimental surfaces.  We can
# post-filter the correctly rotated surfaces using the (x,y,z) position of
# the peaks (so these data go in a separate matrix so we don't end up
# rewriting code).
#
# We are going to filter out contours that are too short or don't have
# at least three points to define a central axis.  We are going to shift
# the line up so that it intersects the peak, just as before.
# 
# Kazic and Stapleton, 21.12.2016


grab_surface_rotation_data <- function(x,y,z,contour_set) {
        contours <- contourLines(x,y,z,levels=get(contour_set))
        the_list <- diffnt(contours)
        if ( nrow(the_list[[1]]) >= 3 ) {

                lsl <- lsfit(x=the_list$axis_points[,1],the_list$axis_points[,2])


# adjust the axis so it lies on the peak; lsl$coefficients[2] is the slope of
# the least squares-fit line

                peak_indices <- arrayInd(which.max(z),c(nrow(z),ncol(z)))
                max_xy <- c(x[peak_indices[1]],y[peak_indices[2]])
                new_intercept <- y[peak_indices[2]] - lsl$coefficients[2]*x[peak_indices[1]]

                slope <- lsl$coefficients[2]

                } else {
                        new_intercept = NA
			slope = NA
                        }
	rot_list <- list("new_intercept"=new_intercept,"slope"=slope)		
        return(rot_list)
        }







################### shape analysis ##################


# for a matrix of strings, where the first column is the name of an
# experimental surface and the second a matching simulated surface, compute:
#
#      a difference surface between the experimental and simulated
#      surfaces, and plot that;
#
#      the vector of differences between the mesh points of the
#      experimental and simulated surfaces, accumulating these in a matrix; and
#
#      a parallel plot of the polygons of the matrix, coloring by
#      experimental surfaces; may be one plot/experimental surface.
#
#
# To prepare the matrix, save the first two columns of the desired output
# from shapes.org as a separate csv file, sort by sweep for easier reading
# of final output, M-x org-table-export to csv and save the csv file =
# match_file. 


# errr, what we really want to do is read in the csv file, and for a vector
# (list?) of the form c(exptl_surface,lower_bd_quant,upper_bd_quant), grab
# the matching surfaces for that exptl_surface between those quantile
# bounds, and pass those surfaces and their Frobenius distances for
# computation.  At least it would be nice to have the Frobenius distances
# on the plots, so rearrange the csv file above to include those, and pass
# to image.plot.
#
# Kazic, 14.9.2016





# diff_matched_shapes("x2as_sim_s53.csv","domed",15,"shapes.org")


diff_matched_shapes <- function(match_file,contour_set,thetadeg,log_file) {

        date <- Sys.time()
        start <- as.integer(date)

        analysis_dir <- "/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts/results/parameter_sweeps/analysis/"
#        analysis_dir <- "~/me/artistry/papers/current/ann_stapleton_ms/our_parts/results/parameter_sweeps/analysis/"
        load(paste0(analysis_dir,"meshes.rda"))

        surfaces <- the_meshes[[1]]
        contours <- the_meshes[[2]]
        meshes   <- the_meshes[[3]]

        x=seq(-42,42,0.5) ;  y=seq(-7.5,7.5,0.5)

        match_mat <- as.matrix(read.csv(file.path(analysis_dir,match_file),header=FALSE,comment.char="#"))

        compare_surfaces(surfaces,match_mat,contours,meshes,contour_set,thetadeg,x,y,start,analysis_dir)
        stop <- as.integer(Sys.time())



        write_to_log_file(analysis_dir,log_file,match_file,date,contour_set,thetadeg,start,stop,nrow(match_mat))

        }











# refactor this

compare_surfaces <- function(surfaces,match_mat,contours,meshes,contour_set,thetadeg,x,y,start,analysis_dir) {

        pal <- colorRampPalette( c("red","yellow","green","blue","violet") )
        colors <- pal(100)
        breaks <- seq(-100,100,2)
        nlevel <- 7



# just load the ones we need!
# > is.vector(levels(match_mat[,1]))

        for ( i in 1:length(surfaces) ) {



# load z matrix of exptl surface

                load(paste0("/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts/results/parameter_sweeps/ann_replots/z",surfaces[i],".rda"))
#                load(paste0("~/me/artistry/papers/current/ann_stapleton_ms/our_parts/results/parameter_sweeps/ann_replots/z",surfaces[i],".rda"))
                es <- as.symbol(paste0("z",surfaces[i]))
                es <- z


# oops, get the (x,y) values, not just the indices of their position in the vectors

                exptl_max_pos_uncentered <- as.vector(which(es == max(es), arr.ind = TRUE))
                exptl_max_pos_centered <- c(x[exptl_max_pos_uncentered[1]],y[exptl_max_pos_uncentered[2]])

                
                peaks <- matrix(exptl_max_pos_centered,nrow=2,ncol=2,byrow=TRUE)
                rm(z)





# construct the matrix of (x,y) values for the mesh points of the
# experimental surface

                exptl_mat <- matrix(c(meshes[[i]][,1],meshes[[i]][,2]),nrow=length(meshes[[i]][,1]),ncol=2)
		

# select all elements of match_mat whose first element matches the
# experimental surface to get the list of simulated surfaces


#                sims <- as.vector(match_mat[which(match_mat[,1]==surfaces[i]),2])
#                if ( length(sims) > 0 ) {



                if ( length(match_mat[which(match_mat[,1]==surfaces[i])]) > 0 ) {

                        sims_mat <- matrix(nrow=length(which(match_mat[,1]==surfaces[i])),ncol=2)
                        sims_mat[,1:2] <- c(match_mat[which(match_mat[,1]==surfaces[i]),2],match_mat[which(match_mat[,1]==surfaces[i]),3])
                        ordered <- sims_mat[order(sims_mat[,2]),]

                        diffmat <- c()
			diffcols <- c()
			
#                        for ( j in 1:length(sims) ) {
                        for ( j in 1:nrow(ordered) ) {


# we don't want to calculate the difference surface and its entailments
# if we've already done so.  So test to see if the file exists; and
# construct the diffmat by appending the vectors of differences between
# pairs of mesh points to a vector, then reshaping that vector into a matrix
# at the end.  This saves slicing out columns from a matrix or worse, cbind.
# And construct the vector of colnames as we go along, so we get the names and the
# final number of columns for reshaping.
#
# Kazic, 14.9.2016


                                file_stem <- paste0(surfaces[i],"_",ordered[j,1])
                                diff_file <- paste0(analysis_dir,file_stem,".rda")



#                                if ( !file.exists(diff_file) ) {

                                        # print(paste0(surfaces[i]," i: ",i," j: ",j," s: ",sims[j]))
                                        # sweep <- regmatches(sims[j],gregexpr("[0-9]+",sims[j],perl=TRUE))[[1]][1]
                                        # sweep_root <- toggle_partition(sweep)
		                        # test <- regmatches(sims[j],gregexpr("[0-9]+",sims[j],perl=TRUE))[[1]][2]
                                        # frobenius <- match_mat[which(match_mat[,2]==sims[j] & match_mat[,1]==surfaces[i]),3]



                                        print(paste0(surfaces[i]," i: ",i," j: ",j," s: ",ordered[j,1]))
                                        sweep <- regmatches(ordered[j,1],gregexpr("[0-9]+",ordered[j,1],perl=TRUE))[[1]][1]
                                        sweep_root <- toggle_partition(sweep)
		                        test <- regmatches(ordered[j,1],gregexpr("[0-9]+",ordered[j,1],perl=TRUE))[[1]][2]
                                        frobenius <- ordered[j,2]



                                        load(paste0(sweep_root,"/results/parameter_sweeps/sweep",sweep,"/z",test,".RData"))
                                        sim_max_pos_uncentered <- as.vector(which(z == max(z), arr.ind = TRUE))
                                        sim_max_pos_centered <- c(x[sim_max_pos_uncentered[1]],y[sim_max_pos_uncentered[2]])


# difference surface

                                        diff_surface <- es - z

                                        peaks[2,] <- sim_max_pos_centered
        	  		        plot_diff_surface(analysis_dir,file_stem,x,y,diff_surface,colors,peaks,breaks,nlevel,frobenius)
                                        save(diff_surface,file=diff_file)
	
	
                                        



# compute euclidean distance between points in the exptl and sim meshes, and save vector
# to diffmat
#
# to get the vector of distances between the two sets of mesh points,
# direct calculation is faster than dist2 and taking the diagonal; see
# http://stackoverflow.com/questions/24746892/how-to-calculate-euclidian-distance-between-two-points-defined-by-matrix-contain
# nb: rowSums assumes nrow >= 2 for each matrix

                                        sim_mesh_list <- define_mesh(x,y,z,contour_set,thetadeg)
                                        diffmat[length(diffmat)+1:nrow(exptl_mat)] <- sqrt(rowSums((exptl_mat - sim_mesh_list[[7]])^2))
                                        diffcols[length(diffcols)+1] <- ordered[j,1]
                                        rm(z)
#	                                }
                                }


# reshape diffmat into the actual matrix
#
# http://stackoverflow.com/questions/4357866/construct-dynamic-sized-array-in-r
# http://www.burns-stat.com/pages/Tutor/R_inferno.pdf


                                dim(diffmat) <- c(nrow(exptl_mat),length(diffcols))
                                rownames(diffmat) <- rownames(meshes[[i]])
			        colnames(diffmat) <- diffcols
                                


# save and plot diffmat

                                diff_stem <- paste0(analysis_dir,"diffmat_",surfaces[i],"_",start)
                                save(diffmat,file=paste0(diff_stem,".rda"))

#
#                               parallel_plot_fit_quality(diff_stem,diffmat,surfaces[i])
#
# rather than a parallel plot, get a stripchart of the diffs with nice
# color coding

                                stripchart_fit_quality(diff_stem,diffmat,surfaces[i])

                                }
	                }
	        }
        
        






# refactored, I couldn't stand it anymore
#
# Stapleton and Kazic, 2.11.2016


# new_compare_surfaces(c("x2a","x3a"),vector,"domed",15)

new_compare_surfaces <- function(target_surfs,vector_sim_surfs,contour_set,thetadeg) {


        analysis_dir <- "/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts/results/parameter_sweeps/analysis/"
#        analysis_dir <- "~/me/artistry/papers/current/ann_stapleton_ms/our_parts/results/parameter_sweeps/analysis/"
        load(paste0(analysis_dir,"meshes.rda"))

        surfaces <- the_meshes[[1]]
        contours <- the_meshes[[2]]
        meshes   <- the_meshes[[3]]

        x=seq(-42,42,0.5) ;  y=seq(-7.5,7.5,0.5)


        pal <- colorRampPalette( c("red","yellow","green","blue","violet") )
        colors <- pal(100)
        breaks <- seq(-100,100,2)
        nlevel <- 7



        date <- Sys.time()
        start <- as.integer(date)


        for ( i in 1:length(target_surfs) ) {


# load z matrix of exptl surface

                load(paste0("/Volumes/b/artistry/papers/current/ann_stapleton_ms/our_parts/results/parameter_sweeps/ann_replots/z",target_surfs[i],".rda"))
#                load(paste0("~/me/artistry/papers/current/ann_stapleton_ms/our_parts/results/parameter_sweeps/ann_replots/z",target_surfs[i],".rda"))
                es <- as.symbol(paste0("z",target_surfs[i]))
                es <- z


# oops, get the (x,y) values, not just the indices of their position in the vectors

                exptl_max_pos_uncentered <- as.vector(which(es == max(es), arr.ind = TRUE))
                exptl_max_pos_centered <- c(x[exptl_max_pos_uncentered[1]],y[exptl_max_pos_uncentered[2]])

                
                peaks <- matrix(exptl_max_pos_centered,nrow=2,ncol=2,byrow=TRUE)
                rm(z)





# construct the matrix of (x,y) values for the mesh points of the
# experimental surface
# rownames a bit of overhead, but might come in handy
#
# damn, they did!
#
# Kazic, 20.11.2016

                k <- which(surfaces==target_surfs[i])
                exptl_mesh <- matrix(c(meshes[[k]][,1],meshes[[k]][,2]),nrow=length(meshes[[k]][,1]),ncol=2)
		rownames(exptl_mesh) <- rownames(meshes[[k]])
		

                if ( length(vector_sim_surfs) > 0 ) {

                        sims_mat <- matrix(nrow=length(vector_sim_surfs),ncol=2)

                        sims_mat[,1] <- vector_sim_surfs
			sims_mat[,2] <- NA

                        diffmat <- c()
			diffcols <- c()


                        for ( j in 1:length(vector_sim_surfs) ) {

                                file_stem <- paste0(target_surfs[i],"_",vector_sim_surfs[j])
                                diff_file <- paste0(analysis_dir,file_stem,".rda")




# we don't want to calculate the difference surface and its entailments
# if we've already done so.  So test to see if the file exists; and
# construct the diffmat by appending the vectors of differences between
# pairs of mesh points to a vector, then reshaping that vector into a matrix
# at the end.  This saves slicing out columns from a matrix or worse, cbind.
# And construct the vector of colnames as we go along, so we get the names and the
# final number of columns for reshaping.
#
# Kazic, 14.9.2016
#
#
# well, most of that's fine, but if the simulated surface doesn't have the same mesh
# points as the experimental one, then its mesh point construction will fail and compromise
# the reshaped matrix.  And similarly, if the diff_file exists, we'd have to hunt out the values
# of the mesh points from somewhere . . . so just recalculate, it'll be fast enough.
#
# Kazic, 5.11.2016


                                print(paste0(target_surfs[i]," i: ",i," j: ",j," s: ",vector_sim_surfs[j]))
                                sweep <- regmatches(vector_sim_surfs[j],gregexpr("[0-9]+",vector_sim_surfs[j],perl=TRUE))[[1]][1]
                                sweep_root <- toggle_partition(sweep)
	                        test <- regmatches(vector_sim_surfs[j],gregexpr("[0-9]+",vector_sim_surfs[j],perl=TRUE))[[1]][2]


                                load(paste0(sweep_root,"/results/parameter_sweeps/sweep",sweep,"/z",test,".RData"))
                                sim_max_pos_uncentered <- as.vector(which(z == max(z), arr.ind = TRUE))
                                sim_max_pos_centered <- c(x[sim_max_pos_uncentered[1]],y[sim_max_pos_uncentered[2]])


# difference surface

                                diff_surface <- es - z

                                peaks[2,] <- sim_max_pos_centered
	
	
                                        



# compute euclidean distance between points in the exptl and sim meshes, and save vector
# to diffmat.
#
# sim_mesh_list[[7]] is the mesh point matrix for the simulated surface
#
# to get the vector of distances between the two sets of mesh points,
# direct calculation is faster than dist2 and taking the diagonal; see
# http://stackoverflow.com/questions/24746892/how-to-calculate-euclidian-distance-between-two-points-defined-by-matrix-contain
# nb: rowSums assumes nrow >= 2 for each matrix
#
# directly calculate the Frobenius distance and overwrite the previous NA in sims_mat with
# the new value
#
# do the plotting here, instead of above, because we only here calculate the Frobenius distance


                                sim_mesh_list <- define_mesh(x,y,z,contour_set,thetadeg)

                                if ( !identical(as.integer(length(sim_mesh_list)),as.integer(0))
                                   & ( identical(sim_mesh_list$contours,contours[which(surfaces==target_surfs[i])][[1]]) ) ) { 
                                        diffmat[length(diffmat)+1:nrow(exptl_mesh)] <- sqrt(rowSums((exptl_mesh - sim_mesh_list[[7]])^2))
                                        frobenius <- norm(exptl_mesh-sim_mesh_list[[7]],type="F")

                                        } else {
					        diffmat[length(diffmat)+1:nrow(exptl_mesh)] <- NA
                                                frobenius <- Inf
                                                }

                                diffcols[length(diffcols)+1] <- vector_sim_surfs[j]
				sims_mat[j,2] <- frobenius

	  		        plot_diff_surface(analysis_dir,file_stem,x,y,diff_surface,colors,peaks,breaks,nlevel,frobenius)
                                save(diff_surface,file=diff_file)

                                rm(z)

                                }



# reshape diffmat into the actual matrix
#
# http://stackoverflow.com/questions/4357866/construct-dynamic-sized-array-in-r
# http://www.burns-stat.com/pages/Tutor/R_inferno.pdf


                                dim(diffmat) <- c(nrow(exptl_mesh),length(diffcols))
                                rownames(diffmat) <- rownames(exptl_mesh)
			        colnames(diffmat) <- diffcols
                                


# save and plot diffmat

                                diff_stem <- paste0(analysis_dir,"diffmat_",target_surfs[i],"_",start)
                                save(diffmat,file=paste0(diff_stem,".rda"))

#
#                               parallel_plot_fit_quality(diff_stem,diffmat,target_surfs[i])
#
# rather than a parallel plot, get a stripchart of the diffs with nice
# color coding


# filter out columns with just NAs, trick from
# http://stackoverflow.com/questions/15968494/how-to-delete-columns-that-contain-only-nas

                                filtered <- diffmat[,colSums(is.na(diffmat)) != nrow(diffmat)]

                                if ( !identical(as.integer(ncol(filtered)),as.integer(0)) ) {
                                        stripchart_fit_quality(diff_stem,filtered,target_surfs[i])
					}
                                }
	                }
	        }
        
        





### refactor it yet again to not bother with mesh points
### and to plot the the simulated surfaces.
#
# descended from new_compare_surfaces/4
#
# Kazic, 2.2.2017


# target_surfs <- c("mo17","b73","x1a","x1b","x2a","x2b","x3a","x3b")
# vector_sim_surfs <- c("s76z3281","s78z1","s78z11","s78z99","s78z109","s78z1373","s78z1383","s78z1471","s78z1481","s78z19209","s78z19219","s78z19307","s78z19317","s78z20581","s78z20591","s78z20679","s78z20689","s78z2956","s78z2966","s78z3054","s78z3064","s78z4328","s78z4338","s78z4426","s78z4436","s78z22164","s78z22174","s78z22262","s78z22272","s78z23536","s78z23546","s78z23634","s78z23644","s78z5911","s78z5921","s78z6009","s78z6019","s78z7283","s78z7293","s78z7381","s78z7391","s78z25119","s78z25129","s78z25217","s78z25227","s78z26491","s78z26501","s78z26589","s78z26599","s78z8866","s78z8876","s78z8964","s78z8974","s78z10238","s78z10248","s78z10336","s78z10346","s78z28074","s78z28084","s78z28172","s78z28182","s78z29446","s78z29456","s78z29544","s78z29554","s78z11821","s78z11824","s78z11919","s78z11922","s78z13193","s78z13196","s78z13291","s78z13294","s78z31029","s78z31032","s78z31127","s78z31130","s78z32401","s78z32404","s78z32499","s78z32502","s78z17731","s78z17734","s78z17829","s78z17832","s78z19103","s78z19106","s78z19201","s78z19204","s78z36939","s78z36942","s78z37037","s78z37040","s78z38311","s78z38314","s78z38409","s78z38412","s78z14776","s78z14779","s78z14874","s78z14877","s78z16148","s78z16151","s78z16246","s78z16249","s78z33984","s78z33987","s78z34082","s78z34085","s78z35356","s78z35359","s78z35454","s78z35457")



# compare_surfaces_n_plot(c("mo17","b73"),c("s78z1","s78z11"))

compare_surfaces_n_plot <- function(target_surfs,vector_sim_surfs) {

        local <- toggle_partition(0)
        pal3d <- load_standard_plotting_parameters(local)


# for 2d difference surfaces

        pal2d <- colorRampPalette( c("red","yellow","green","blue","violet") )
        colors <- pal2d(100)
        breaks <- seq(-100,100,2)
        nlevel <- 7



        analysis_dir <- paste0(local,"/results/parameter_sweeps/analysis/")
        ann_dir <- paste0(local,"/results/parameter_sweeps/ann_replots/")

        load(paste0(ann_dir,"exptl_meshes.rda"))
        load(paste0(ann_dir,"mo17_mesh.rda"))
        load(paste0(ann_dir,"exptl_meshes_tops.rda"))



#        x=seq(-42,42,0.5) ;  y=seq(-7.5,7.5,0.5)
        xs <- rep(x,31) ; ys <- rep(y,169)



        date <- Sys.time()
        start <- as.integer(date)


        for ( i in 1:length(target_surfs) ) {



# grab or form the correct matrix of mesh points

                if ( target_surfs[i] == "mo17" ) {

			mesh_xs <- colnames(mo17_mesh)[1]
			mesh_ys <- colnames(mo17_mesh)[2]

                        tmpz <- rownames(mo17_mesh)
                        gz <- regmatches(tmpz,gregexpr("[0-9\\.]+",tmpz,perl=TRUE))
                        z <- as.numeric(unlist(gz)[seq(1,length(gz)*2,2)])
                        mesh <- cbind(z,mo17_mesh)

		        } else {
			        mesh_xs <- paste0(target_surfs[i],":x")
			        mesh_ys <- paste0(target_surfs[i],":y")
				tops <- nrow(exptl_meshes) - exptl_meshes_tops[mesh_xs] + 1
				first <- exptl_meshes_tops[mesh_xs]
				last <- nrow(exptl_meshes)
                                mesh <- matrix(c(exptl_meshes[first:last,c("z",mesh_xs,mesh_ys)]),nrow=tops,ncol=3,byrow=FALSE)
                                colnames(mesh) <- c("z",mesh_xs,mesh_ys)
                                }




# load z matrix of exptl surface

                load(paste0(ann_dir,"z",target_surfs[i],".rda"))
                es <- as.symbol(paste0("z",target_surfs[i]))
                es <- z


# oops, get the (x,y) values, not just the indices of their position in the vectors

                exptl_max_pos_uncentered <- as.vector(which(es == max(es), arr.ind = TRUE))
                exptl_max_pos_centered <- c(x[exptl_max_pos_uncentered[1]],y[exptl_max_pos_uncentered[2]])

                
                peaks <- matrix(exptl_max_pos_centered,nrow=2,ncol=2,byrow=TRUE)
                rm(z)




# now iterate over all simulated surfaces and compare these to the chosen
# experimental surface

                if ( length(vector_sim_surfs) > 0 ) {
		
                        diffmat <- c()
			diffnames <- c()


                        for ( j in 1:length(vector_sim_surfs) ) {

                                file_stem <- paste0(target_surfs[i],"_",vector_sim_surfs[j])
                                diff_file <- paste0(analysis_dir,file_stem,".rda")



                                print(paste0(target_surfs[i]," i: ",i," j: ",j," s: ",vector_sim_surfs[j]))
                                sweep <- regmatches(vector_sim_surfs[j],gregexpr("[0-9]+",vector_sim_surfs[j],perl=TRUE))[[1]][1]
                                sweep_root <- toggle_partition(sweep)
	                        test <- regmatches(vector_sim_surfs[j],gregexpr("[0-9]+",vector_sim_surfs[j],perl=TRUE))[[1]][2]


                                load(paste0(sweep_root,"/results/parameter_sweeps/sweep",sweep,"/z",test,".RData"))
                                sim_max_pos_uncentered <- as.vector(which(z == max(z), arr.ind = TRUE))
                                sim_max_pos_centered <- c(x[sim_max_pos_uncentered[1]],y[sim_max_pos_uncentered[2]])


# difference surface

                                diff_surface <- es - z

                                peaks[2,] <- sim_max_pos_centered
	
	


# Euclidean distance between experimental mesh points and the z positions of
# the (x,y) coordinates for each of the simulated surfaces

                                interpolated <- interpp(xs,ys,c(z),mesh[,mesh_xs],mesh[,mesh_ys],linear=TRUE,extrap=FALSE,duplicate="error")

                                distances <- sqrt((interpolated$z)^2 + (mesh[,"z"])^2)


				next_index <- length(diffmat) + 1
                                so_far <- length(diffmat) + length(distances)
                                diffmat[next_index:so_far] <- distances
                                diffnames <- c(diffnames,vector_sim_surfs[j])


                                if ( i == 1 ) {
                                        plot_sim_surface(analysis_dir,x,y,z,pal3d,vector_sim_surfs[j])
                                        }
					
	  		        plot_diff_surface(analysis_dir,file_stem,x,y,diff_surface,colors,peaks,breaks,nlevel,0)
                                save(diff_surface,file=diff_file)

                                rm(z)

                                }



# reshape diffmat into the actual matrix
#
# http://stackoverflow.com/questions/4357866/construct-dynamic-sized-array-in-r
# http://www.burns-stat.com/pages/Tutor/R_inferno.pdf

                        diffcols <- as.integer(length(diffmat) / nrow(mesh))
                        dim(diffmat) <- c(nrow(mesh),diffcols)
                        rownames(diffmat) <- rownames(mesh)
		        colnames(diffmat) <- diffnames
                                

  
# save and plot diffmat

                        diff_stem <- paste0(analysis_dir,"diffmat_",target_surfs[i],"_",start)
                        save(diffmat,file=paste0(diff_stem,".rda"))


# filter out columns with just NAs, trick from
# http://stackoverflow.com/questions/15968494/how-to-delete-columns-that-contain-only-nas

                        filtered <- diffmat[,colSums(is.na(diffmat)) != nrow(diffmat)]

                        if ( !identical(as.integer(ncol(filtered)),as.integer(0)) ) {
                                stripchart_fit_quality(diff_stem,filtered,target_surfs[i])
				}
                        }
	        }
	}
        
        








        
        
	
# http://www.endmemo.com/program/R/points.php        


# image.plot from package fields puts in a perfectly adequate color bar for
# our purposes.  hoooray!

# plot_diff_surface("../results/parameter_sweeps/analysis/","foo",seq(-42,42,0.5),seq(-7.5,7.5,0.5),z,colors,peaks,breaks,nlevel,frobenius)


plot_diff_surface <- function(analysis_dir,file_stem,x,y,diff_surface,colors,peaks,breaks,nlevel,frobenius) {
        png(paste0(analysis_dir,file_stem,"_diffsurf.png"))

#         image.plot(x,y,diff_surface,add=FALSE,col=colors,graphics.reset=TRUE,horizontal=TRUE,breaks=breaks,nlevel=nlevel,xlab="water",ylab="nitrogen",main=paste0(file_stem," difference surface, Fro = "))


        image.plot(x,y,diff_surface,add=FALSE,col=colors,breaks=breaks,nlevel=nlevel,graphics.reset=FALSE,horizontal=TRUE,xlab="water",ylab="nitrogen",main=paste0(file_stem," diff surf, Fro = ",frobenius," circ = exptl, sq = sim"))
        points(peaks[,1],peaks[,2],pch=c(16,15),cex=2)

        dev.off()
        }












# must transpose matrix for parallel plot as variables are now the mesh points

parallel_plot_fit_quality <- function(diff_stem,diffmat,surface_name) {
        png(paste0(diff_stem,"_pplot.png"))
	title <- paste0("distances from mesh points of ",surface_name)
        parcoord(t(diffmat),col=c("violet","blue","green","yellow","red"),lty=1,var.label=TRUE,main=title)
        dev.off()
        }




# nicer, easier to read, I hope
#
# the mesh points for each surface are listed from the peak
# going outward.  The lowest contour for each surface is closest
# to the minima (for the domed surfaces), and each contour has five
# mesh points.  So just color the last five points in the diffmat
# differently.
#
# Kazic, 29.9.2016


# use stripChart and get rid of confidence intervals, silly here; but
# better than iterating over stripchart:
# https://www.rdocumentation.org/packages/EnvStats/versions/2.1.1/topics/stripChart
#
# label rotation trick from:
# https://www.r-bloggers.com/rotated-axis-labels-in-r-plots/


# hmmm, rotate 90???

stripchart_fit_quality <- function(diff_stem,diffmat,surface_name) {

        rows <- nrow(diffmat)
        green_pts <- rows - 10
	green_ptsp <- green_pts + 1
        blue_pts <- rows - 5
	blue_ptsp <- blue_pts + 1
        red_pts <- rows
	xlocs <- seq(1,ncol(diffmat),1)
	xlabels <- colnames(diffmat)

        png(paste0(diff_stem,"_scplt.png"))
	title <- paste0("distances from mesh points of ",surface_name)

        stripChart(diffmat[1:green_pts,],col="green",main=title,vertical=TRUE,show.ci=FALSE,n.text="none",location.scale.text="none",xaxt="n",xlab="")
        axis(1,at=xlocs,labels=FALSE)
        text(x=xlocs,y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),labels=xlabels,srt=90,adj=0.9,xpd=TRUE,cex=0.5)

        stripChart(diffmat[green_ptsp:blue_pts,],col="blue",pch=16,add=TRUE,vertical=TRUE,show.ci=FALSE,n.text="none",location.scale.text="none")
        stripChart(diffmat[blue_ptsp:red_pts,],col="red",pch=16,add=TRUE,vertical=TRUE,show.ci=FALSE,n.text="none",location.scale.text="none")


        dev.off()
        }















# write basic statistics out to shapes.org

write_to_log_file <- function(analysis_dir,log_file,match_file,date,contour_set,thetadeg,start,stop,num_compared) {

        zz <- file(paste0(analysis_dir,log_file), "a")
        elapsed <- start - stop
        cat(paste0("\n\n** shape differencing timestamp ",start,"\n\nstarted on date: ",date, "\nexecution time (clock): ",elapsed,"\nfor simulated surfaces found in file: ",match_file,"\ncontour_set: ",contour_set,"\nthetadeg: ",thetadeg,"\ntotal number simulated surfaces compared: ",num_compared,"\n\n"),file=zz)

        cat(paste0("all output files are found in ",analysis_dir,":\n"),file=zz)
        cat("   + difference surfaces between experimental and simulated are in files of the form: EXPTLSURF_SIMSURF.rda\n",file=zz)
	cat("   + their plots are in files of the form: EXPTLSURF_SIMSURF_diffsurf.png\n",file=zz)
	cat("   + matrices of Euclidean distances between mesh points for experimental and matching simulated surfaces are in files of the form diffmat_EXPTLSURF_STARTINGTIMESTAMP.rda\n",file=zz)	
	cat("   + parallel plots of Euclidean distances are in files of the form EXPTLSURF_STARTINGTIMESTAMP_pplot.png\n\n\n",file=zz)	
        close(zz)
        }







# given a list of experimental surfaces, load their individually named
# big_list, and add their meshes and rays to two combined matrices.
#
#
# notice the use of get to grab the object associated with the constructed
# name:  this is the ``inverse'' of assign.
#
# Kazic, 18.11.2016

# x2b on no data

bind_mesh_rays <- function() {
        surfaces <- c("b73","x1a","x1b","x2a","x2b","x3a","x3b")
	sweep_root <- toggle_partition(0)
	setwd(paste0(sweep_root,"/simulations"))
        ann_replots <- "../results/parameter_sweeps/ann_replots/"

        for ( i in 1:length(surfaces) ) {
                load(paste0(ann_replots,"z",surfaces[i],"_mesh.rda"))
                surf_list <- get(paste0(surfaces[i],"_big_list"))

                if ( identical(as.numeric(i),1) ) {
                        exptl_meshes <- matrix(NA,nrow=nrow(surf_list$mesh),ncol=(length(surfaces)*2 + 1))
                        rownames(exptl_meshes) <- rownames(surf_list$mesh)
                        colnames(exptl_meshes) <- c("z",interleave(paste0(surfaces,":x"),paste0(surfaces,":y")))
                        exptl_meshes[,1] <- as.numeric(unlist(regmatches(rownames(exptl_meshes),gregexpr("[0-9]{2}",rownames(exptl_meshes),perl=TRUE))))

                        b <- length(rownames(exptl_meshes))

                        exptl_rays <- matrix(NA,nrow=5,ncol=length(surfaces)*2)
                        rownames(exptl_rays) <- colnames(b73_big_list$rays)
                        colnames(exptl_rays) <- c(interleave(paste0(surfaces,":m"),paste0(surfaces,":b")))
                        }

                t <- length(rownames(exptl_meshes))-length(rownames(surf_list$mesh)) + 1
                col1 <- i * 2
		col2 <- col1 + 1
		col3 <- col1 - 1
		col4 <- col3 + 1
                exptl_meshes[t:b,col1:col2] <- surf_list$mesh[,1:2]
                exptl_rays[,col3:col4] <- c(surf_list$rays[2,],surf_list$rays[1,])
		print(paste0("i: ",i," surface: ",surfaces[i]," t: ",t," b: ",b," col3: ",col3," col4: ",col4))
                }


        save(exptl_meshes,file=paste0(ann_replots,"exptl_meshes.rda"))
        save(exptl_rays,file=paste0(ann_replots,"exptl_rays.rda"))	
        }









# cribbed from an old mailing list post

interleave <- function(v1,v2) {
        ord1 <- 2*(1:length(v1))-1
        ord2 <- 2*(1:length(v2))
        c(v1,v2)[order(c(ord1,ord2))]
        }






# this is a revision of compare_shapes/10
#
#
# we realized that we are once again comparing apples to oranges.  The
# mesh points for each simulated surfaces shift as the rays shift.  This is
# deliberate, but it means that the distances between mesh points are not
# the normals from the experimental mesh points to the simulated mesh
# points. So the shape matching looks good numerically, even if we look at
# the distributions of the distances between an experimental and a
# simulated surface; but when we look at the difference plots, or the
# contour plot for a simulated surface, they're way off.  In effect, by
# permitting the simulated surface's mesh points to rotate relative to
# those of the experimental surface, we have sheared the vectors of
# distances.
#
# We see two options.
#
#     1.  use the (x,y) values from the experimental mesh points and find z
#         for the simulated surface, then calculate z_exptl - z_sim.  This
#         is the "fixed mesh point" option.  It should be appropriately signed. 
#
#     2.  use the rays of the experimental surface on the simulated
#         surface, compute the intersections of the contours and fixed
#         rays, use those as the mesh points, and then compute the
#         distances.  This is the "fixed ray" option.
#
#
#
# I've constructed two matrices with data from the experimental surfaces,
# one for the mesh points and the other for the rays, using
# bind_mesh_rays/0 above.
#
#
# I considered computing the normals from the experimental mesh points to the
# simulated surface, and then the magnitude of those normals.  I thought
# this would be equivalent to the "fixed mesh point" option, but this is
# untrue.  As the surface bends, the direction of the normal will
# change. If we're interested in the bending, we can calculate the normals
# or a grid of derivatives and do arrow plots.  But otherwise, no.
#
#
# The "fixed ray" option will still permit shear of the contours relative to
# those of the experimental surfaces.  I think the place where we want to
# permit this rotation is among the experimental surfaces, not between
# experimental and simulated surfaces, where we're looking for the best
# match.
#
#
# Once we have the distances, then plotting them in a stripchart and
# computing the variance for each comparison as a summary statistic makes
# sense.  I re-checked the mesh points, and they nicely avoid the nasty
# corners and edges as I remembered.
#
# Kazic, 20.11.2016



# compare_shapes_by_z_dists("1478728934","domed","shapes.org")


compare_shapes_by_z_dists <- function(surfaces_timestamp,contour_set,log_file) {

       date <- Sys.time()
       start <- as.integer(date)
     

       sweep_root <- toggle_partition(0)
       load(paste0(sweep_root,"/results/parameter_sweeps/ann_replots/exptl_meshes.rda"))
       load(paste0(sweep_root,"/results/parameter_sweeps/ann_replots/exptl_meshes_tops.rda"))
       analysis_dir <- paste0(sweep_root,"/results/parameter_sweeps/analysis/")
       surfaces_file <- list.files(analysis_dir)[grep(glob2rx(paste0("surfaces*",surfaces_timestamp,".RData")),list.files(analysis_dir))]
       load(paste0(analysis_dir,surfaces_file))
       load(paste0(analysis_dir,"meshes.rda"))

       surfaces <- the_meshes[[1]]
       contours <- the_meshes[[2]]



# hit_list is stored in ../results/parameter_sweeps/analysis/surfaces*TIMESTAMP.RData
#
# may be faster than looping through hit_list, but is certainly more transparent
#
# neither of these sorts the surfaces, but hit_list is built by sweep, so each sweep's 
# surfaces are adjacent in the vector simulated.  Further sorting might result in a speed
# improvement in disk access.
#
# Kazic, 28.11.2016

       flattened <- unlist(hit_list)
       sim <- grep("s[0-9]+z[0-9]+",flattened,perl=TRUE,value=TRUE)
       simulated <- sort(sim)



       summary_matrix <- compare_surfaces_by_z_dists(simulated,surfaces,contours,contour_set,exptl_meshes,exptl_meshes_tops,surfaces_timestamp,analysis_dir)

       stop <- as.integer(Sys.time())

       write_summary_matches(analysis_dir,log_file,date,start,stop,surfaces_timestamp,surfaces_file,surfaces,contour_set,summary_matrix)
       }













# this is a revision of new_compare_surfaces/4
#
# this version shunts surfaces that yield internal NAs on interpolation to
# a separate matrix
#
# Kazic, 3.12.2016


compare_surfaces_by_z_dists <- function(simulated,surfaces,contours,contour_set,exptl_meshes,exptl_meshes_tops,surfaces_timestamp,analysis_dir) {

        x=seq(-42,42,0.5) ;  y=seq(-7.5,7.5,0.5)
        xs <- rep(x,31) ; ys <- rep(y,169)
        b <- nrow(exptl_meshes)

	

# set up results vectors here for final reshaping at the end
#
# a version for the good ones, and another to hold the surfaces that have NAs
# on interpolation

        dists_vec <- vector(mode="numeric", length=0)
        vars_vec <-  vector(mode="numeric", length=0)
        match_names_vec <- vector(mode="character",length=0)
        dists_matrix <- c()


	na_dists_vec <- vector(mode="numeric", length=0)
        na_match_names_vec <- vector(mode="character",length=0)
	na_dists_matrix <- c()
	

        for ( i in 1:length(simulated) ) {

                sweep <- regmatches(simulated[i],gregexpr("[0-9]+",simulated[i],perl=TRUE))[[1]][1]
	        test <- regmatches(simulated[i],gregexpr("[0-9]+",simulated[i],perl=TRUE))[[1]][2]
                sweep_root <- toggle_partition(sweep)
                load(paste0(sweep_root,"/results/parameter_sweeps/sweep",sweep,"/z",test,".RData"))



# produce the list of contours for the simulated surface
# check that against the contours for the exptl surfaces

                test_contours <- contourLines(x,y,z,levels=get(contour_set))
                if ( length(test_contours) >= 1 ) {
                        the_list <- diffnt(test_contours)
                        if ( !identical(as.integer(length(the_list)),as.integer(0)) ) {

                                pruned_contours <- vector(mode="numeric", length=0)

                                for ( j in 1:length(the_list$new_contours) ) { pruned_contours[j] <- the_list$new_contours[[j]][[1]] }



# returns the vector of indices of experimental surfaces, not the
# names of the surfaces; just overwrite match_vec with the names

                                match_vec <- c()

                                for ( k in 1:length(surfaces) ) {
                                        if ( identical(as.numeric(pruned_contours),contours[[k]]) ) {
                                                len <- length(match_vec)
				                match_vec[len+1] <- k
				                }
                                        }
                       
                                match_vec <- surfaces[match_vec]




                                if ( !identical(as.numeric(length(match_vec)),0) ) {

                                       num_matching_surfs <- length(match_vec)


# for those sim surfaces that have matching contours, grab the
# corresponding mesh points and compute the distance along the z axis
#
# use match_vec to construct the list of columns to subset exptl_meshes,
# and interpolate over the subsetted matrix
#
# nifty trick for generating ordered combined keys is at
# http://stackoverflow.com/questions/16143700/pasting-two-vectors-with-combinations-of-all-vectors-elements

                                       exptl_cols <- as.vector(t(outer(match_vec,c(":x",":y"),paste0)))


                                       for ( qq in seq(1,length(exptl_cols),2) ) {


# iterate over the matching exptl surfaces, subtract
# interpp can't tolerate NAs, phooey

                                               t <- exptl_meshes_tops[exptl_cols[qq]]

#		                               print(paste0("sim: s",sweep,"z",test," exptl: ",exptl_cols[qq]," t: ",t))


# if the nominal value of z is off the simulated surface, interpp will return NA.
# asking it to extrapolate produces nonsense.  So decided to preserve the NAs -- they're
# informative after all -- and then color-code surfaces containing at least one NA.
#
# Kazic, 25.11.2016

                                               interpolated <- interpp(xs,ys,c(z),exptl_meshes[t:b,exptl_cols[qq]],exptl_meshes[t:b,exptl_cols[qq+1]],linear=TRUE,extrap=FALSE,duplicate="error")





# write out these magnitudes to a summary matrix using the reshaping trick;
# be careful to preserve upper NAs and surfaces' names
#
# eucl = sqrt((exp_z - sim_z)^2)
#      = abs(exp_z - sim_z)
#
# if exp > sim, then ind = -eucl
# if exp < sim, then ind = +eucl
#
# take the sqrt of the absolute value of the difference of squares
# to avoid getting NaNs (R's response to the sqrt of a negative number)
#
# the mask will tell us what we want to know, anyway




# do we have interior NAs in the interpolated?  if so, write
# those values to a separate matrix for later consideration.
#
# The test is based on the fact that sort by default removes NAs, so
# if any are present then the lengths of the two vectors will be unequal:
# see
# http://www.ats.ucla.edu/stat/r/faq/missing.htm


                                                na_test <- sort(interpolated$z)
                                                if ( identical(length(na_test),length(c(interpolated$z))) ) {

#                                                        eucl_dists <- sqrt((exptl_meshes[t:b,"z"] - interpolated$z)^2)
                                                        eucl_dists <- abs(exptl_meshes[t:b,"z"] - interpolated$z)
                                                        mask <- ifelse(as.numeric(exptl_meshes[t:b,"z"] < interpolated$z),1,-1)
                                                        signed_dists <- eucl_dists * mask


                                                        dists_so_far <- as.numeric(length(dists_vec))
						        vars_so_far <- as.numeric(length(vars_vec))


						        na_start <- dists_so_far + 1
						        na_end <- dists_so_far + t - 1
						        signed_start <- na_end + 1
						        signed_end <- signed_start + length(signed_dists) - 1
						        dists_vec[na_start:na_end] <- NA
						        dists_vec[signed_start:signed_end] <- signed_dists



# compute variance of each distribution of distances, write summary matrix
# of those results with surfaces' names and exptl names, with key as above

                                                        vars_vec[vars_so_far+1] <- var(signed_dists,na.rm=TRUE)

                                                        matches_so_far <- length(match_names_vec)
						        names_index <- (qq %/% 2) + 1
                                                        match_names_vec[matches_so_far+1] <- paste0(match_vec[names_index],":",simulated[i])
						        } else {

                                                                na_dists_so_far <- as.numeric(length(na_dists_vec))
						                na_na_start <- na_dists_so_far + 1
						                na_na_end <- na_dists_so_far + t - 1
						                na_dists_vec[na_na_start:na_na_end] <- NA
								na_na_next <- na_na_end + 1
								na_na_penult <- na_na_next + length(interpolated$z) -1
						                na_dists_vec[na_na_next:na_na_penult] <- interpolated$z


                                                                na_matches_so_far <- length(na_match_names_vec)
						                na_names_index <- (qq %/% 2) + 1
                                                                na_match_names_vec[na_matches_so_far+1] <- paste0(match_vec[na_names_index],":",simulated[i])


#                                                                print(paste0("len inter: ",length(interpolated$z),"\n"))
#								 print(interpolated$z)
#                                                                print(paste0("nns: ",na_na_start," nne: ",na_na_end," nnn: ",na_na_next," nnp: ",na_na_penult,"\n"))
#
#                                                                cols_so_far <- length(na_dists_vec)/b
#                                                                dum <- matrix(na_dists_vec,nrow=b,ncol=cols_so_far)
#
#                                                                colnames(dum) <- na_match_names_vec
#								 print(dum)
#								 cat("\n\n\n")
                                                                }
                                                }
                                        }
                                }
                        }
                rm(z)	
                }




# reshape and label the dists_matrix, add in the vector of variances

        dists_rows <- nrow(exptl_meshes)
        tmp1 <- matrix(dists_vec,nrow=dists_rows,ncol=length(match_names_vec),byrow=FALSE)
	tmp2 <- matrix(NA,nrow=dists_rows+1,ncol=length(match_names_vec))
	tmp2[1:dists_rows,] <- tmp1
	tmp2[dists_rows+1,] <- vars_vec
        rownames(tmp2) <- c(rownames(exptl_meshes),"variances")
	colnames(tmp2) <- match_names_vec


# at the very end, sort the reshaped matrices by the column names to group by
# experimental surface.

        sorted <- sort(match_names_vec)
	dists_matrix <- tmp2[,sorted]


        labelled_mat <- paste0("shape_dists_matrix_",surfaces_timestamp)
	assign(labelled_mat,dists_matrix)
        mat_stem <- paste0(analysis_dir,"shape_dists_",surfaces_timestamp)		
        save(list=labelled_mat,file=paste0(mat_stem,".rda"))





# and save those surfaces with internal NAs

        tmp3 <- matrix(na_dists_vec,nrow=dists_rows,ncol=length(na_match_names_vec),byrow=FALSE)
	rownames(tmp3) <- rownames(exptl_meshes)
	colnames(tmp3) <- na_match_names_vec
        na_sorted <- sort(na_match_names_vec)
        na_dists_matrix <- tmp3[,na_sorted]


        na_labelled_mat <- paste0("shape_dists_matrix_na_",surfaces_timestamp)
	assign(na_labelled_mat,na_dists_matrix)
        na_mat_stem <- paste0(analysis_dir,"shape_dists_na_",surfaces_timestamp)		
        save(list=na_labelled_mat,file=paste0(na_mat_stem,".rda"))




# plot the summary matrix as a strip chart, partitioning by exptl surface

        summary_matrix <- partition_n_stripchart(mat_stem,dists_matrix,surfaces,b)

        return(summary_matrix)
        }







# partition the labelled_mat by experimental surface, and for each experimental
# surface, plot each simulated surface's distances from 0.  Mark the variance with a different
# glyph.
#
# Black glyphs for the mesh points indicate the first surface has no interior NAs, which are
# produced by interpp when the interpolated z value is off the simulated surface.  This is a very crude
# test because it assumes the the simulated surfaces in a partition are relatively uniform, and they may
# well not be.  But I can't get interpp to extrapolate for me, and I don't want to stripchart matrices
# column by column.
#
#                        nas <- which(is.na(labelled_mat[1:rows-1,partition_indices[1]]))
#                        if ( all(nas==seq(1,length(nas),1)) ) {
#
# Revised to test for NAs interior to the interpolated points over the entire partition.
# As before, black glyphs indicate no interior NAs, blue indicate at least one.
#
# Kazic, 26.11.2016
#
#
# Revised again because we are now writing surfaces with interior NAs due to interpolation
# to a separate matrix.  But changed all the colors from black to blue as it's a bit nicer.
#
# Kazic, 3.12.2016
#
#
# added in a basic boxplot to see if that's an improvement
#
# Kazic, 4.12.2016


partition_n_stripchart <- function(mat_stem,dists_matrix,surfaces,b) {

        summary_matrix <- matrix(NA,nrow=4,ncol=length(surfaces))
	rownames(summary_matrix) <- c("num_surfs","min_dist","max_dist","mean_var")
        colnames(summary_matrix) <- surfaces


        rows <- nrow(dists_matrix)
        labels <- colnames(dists_matrix)


        for ( i in 1:length(surfaces) ) {

                partition_indices <- grep(surfaces[i],labels,perl=TRUE,value=FALSE)
                partition_names <- labels[partition_indices]
		partition_names <- sub(paste0(surfaces[i],":"),"",partition_names, perl=TRUE)
                num_surfaces <- length(partition_indices)
                if ( !identical(as.numeric(num_surfaces),0) ) {
                        min <- min(dists_matrix[1:rows-1,partition_indices],na.rm=TRUE)
                        max <- max(dists_matrix[1:rows-1,partition_indices],na.rm=TRUE)
                        mean_var <- mean(dists_matrix[rows,partition_indices],na.rm=TRUE)
                        } else {
                                min <- 0
                                max <- 0
                                mean_var <- 0
                                }


		summary_matrix[,i] <- c(num_surfaces,min,max,mean_var)


                if ( !identical(as.numeric(num_surfaces),0) ) {

                        xlocs <- seq(1,num_surfaces,1)
			
                	png(paste0(mat_stem,"_",surfaces[i],"_scplt.png"),width=11,height=8,units="in",res=72)
			title <- paste0("Euclidean distances from mesh points for ",num_surfaces," simulated surfaces for experimental surface ",surfaces[i])


                	stripChart(dists_matrix[1:rows-1,partition_indices],col="blue",main=title,ylim=c(-50,50),vertical=TRUE,show.ci=FALSE,n.text="none",location.scale.text="none",xaxt="n",xlab="",na.action = NULL)


                	stripChart(t(dists_matrix[rows,partition_indices]),col="red",pch=16,add=TRUE,vertical=TRUE,show.ci=FALSE,n.text="none",location.scale.text="none",xaxt="n",xlab="",na.action = NULL)
                	axis(1,at=xlocs,labels=FALSE)
                        text(x=xlocs,y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),labels=partition_names,srt=90,adj=0.9,xpd=TRUE,cex=0.35)
		        
                	dev.off()


# boxplot

                	png(paste0(mat_stem,"_",surfaces[i],"_boxplt.png"),width=11,height=8,units="in",res=72)
			boxplot.matrix(dists_matrix,main=title,notch=TRUE,ylim=c(-50,50))
			dev.off()
			}
	        }
		
	return(summary_matrix)	
        }		


















# write a summary table to the shapes.org log file of exptl surface, simulated
# surfaces and their variances, with key as above

write_summary_matches <- function(analysis_dir,log_file,date,start,stop,surfaces_timestamp,surfaces_file,surfaces,contour_set,summary_matrix) {

        zz <- file(paste0(analysis_dir,log_file), "a")
        elapsed <- stop - start
	days <- floor(((elapsed/60)/60)/24)
        hms <- format(.POSIXct(elapsed,tz="GMT"), "%H:%M:%S")	
	
        cat(paste0("\n\n\n* Shape Matching by Distance along z at Mesh Points, comparison timestamp ",start,"\n\nstarted on date: ",date, "\nexecution time (clock): ",elapsed," sec or ",days," days ",hms," HH:MM:SS\n\n"),file=zz)

        cat(paste0("\nfor simulated surfaces found at timestamp: ",surfaces_timestamp,"\ncontour_set: ",contour_set,"\n\nsurfaces: "),file=zz)
        cat(surfaces,file=zz)
	cat(paste0("\nsurfaces file: ",surfaces_file,"\n\n"),file=zz)
        cat(paste0("\n\n\nall output files are found in ",analysis_dir,":\n"),file=zz)
        cat("   + matrix of Euclidean z-distances and their variances is in\n",file=zz)
        cat("     a file of the form shape_dists_SURFACES_TIMESTAMP.rda\n",file=zz)
	cat("   + stripchart plots for each allele are in files of the form\n",file=zz)
	cat("     shape_dists_SURFACES_TIMESTAMP_EXPTLSURF_scplt.png\n\n\n",file=zz)

        cat(paste0("#+tblname: shape_matching_by_z_dists_for_surfaces_timestamp_",surfaces_timestamp,"\n"),file=zz)
        write.table(summary_matrix,file=zz,sep=" | ",row.names=TRUE,col.names=TRUE)

        close(zz)
        }










# given a hit_list from a surfaces*.rda file, sort the simulated surfaces into those
# that match at least one experimental surface without (sheep, gooduns) and with (goats, baduns) internal NAs.
#
# Kazic, 3.12.2016

sort_sheep_from_goats <- function(simulated) {

        results_root <- toggle_partition(0)
        load(paste0(results_root,"/results/parameter_sweeps/ann_replots/exptl_meshes.rda"))
        load(paste0(results_root,"/results/parameter_sweeps/ann_replots/exptl_meshes_tops.rda"))	
        analysis_dir <- paste0(results_root,"/results/parameter_sweeps/analysis/")
        b <- nrow(exptl_meshes)

        load(paste0(analysis_dir,"meshes.rda"))
        surfaces <- the_meshes[[1]]
        contours <- the_meshes[[2]]	
        contour_set <- "domed"


        x=seq(-42,42,0.5) ;  y=seq(-7.5,7.5,0.5)
        xs <- rep(x,31) ; ys <- rep(y,169)


        gooduns <- vector(mode="character", length=0)
        baduns <- vector(mode="character", length=0)
	
        for ( i in 1:length(simulated) ) {	
                sweep <- regmatches(simulated[i],gregexpr("[0-9]+",simulated[i],perl=TRUE))[[1]][1]
	        test <- regmatches(simulated[i],gregexpr("[0-9]+",simulated[i],perl=TRUE))[[1]][2]
                sweep_root <- toggle_partition(sweep)
                load(paste0(sweep_root,"/results/parameter_sweeps/sweep",sweep,"/z",test,".RData"))



# produce the list of contours for the simulated surface
# check that against the contours for the exptl surfaces

                test_contours <- contourLines(x,y,z,levels=get(contour_set))
                if ( length(test_contours) >= 1 ) {
                        the_list <- diffnt(test_contours)
                        if ( !identical(as.integer(length(the_list)),as.integer(0)) ) {

                                pruned_contours <- vector(mode="numeric", length=0)

                                for ( j in 1:length(the_list$new_contours) ) { pruned_contours[j] <- the_list$new_contours[[j]][[1]] }

                                match_vec <- c()

                                for ( k in 1:length(surfaces) ) {
                                        if ( identical(as.numeric(pruned_contours),contours[[k]]) ) {
                                                len <- length(match_vec)
				                match_vec[len+1] <- k
				                }
                                        }
                       
                                match_vec <- surfaces[match_vec]



                                if ( !identical(as.integer(length(match_vec)),as.integer(0)) ) {

                                       num_matching_surfs <- length(match_vec)
                                       exptl_cols <- as.vector(t(outer(match_vec,c(":x",":y"),paste0)))


                                       for ( qq in seq(1,length(exptl_cols),2) ) {

                                               t <- exptl_meshes_tops[exptl_cols[qq]]
					       
                                               interpolated <- interpp(xs,ys,c(z),exptl_meshes[t:b,exptl_cols[qq]],exptl_meshes[t:b,exptl_cols[qq+1]],linear=TRUE,extrap=FALSE,duplicate="error")


                                                na_test <- sort(interpolated$z)
						names_index <- (qq %/% 2) + 1
						surf_name <- match_vec[names_index]
						
						if ( identical(length(na_test),length(c(interpolated$z))) ) {
                                                        next_goodun <- length(gooduns) + 1
							gooduns[next_goodun] <- paste0(surf_name,":",simulated[i])
                                                        } else {
                                                                next_badun <- length(baduns) + 1
                                                                baduns[next_badun] <- paste0(surf_name,":",simulated[i])
                                                                }
						}
                                        }
				}
			}
		        rm(z)	
		}
		save(gooduns,file=paste0(analysis_dir,"/gooduns.rda"))
		save(baduns,file=paste0(analysis_dir,"/baduns.rda"))
        }













# called from compare_surfaces_n_plot
#
# note a lot of the parameters are in the global environment, placed
# there by load_standard_plotting_parameters


plot_sim_surface <- function(save_dir,x,y,z,pal3d,file_stem) {

		
        persp3d(x,y,z,col=pal3d,xlim=c(-45,45),ylim=c(-8,8),zlim=zlimits,forceClipregion=TRUE,xlab="",ylab="",zlab="",axes=FALSE,box=FALSE,lit=TRUE,specular="black")

        segments3d(x=as.vector(t(m_for_z[,c(1,4)])),y=as.vector(t(m_for_z[,c(2,5)])),z=as.vector(t(m_for_z[,c(3,6)])))
        axis3d('x--',pos=c(-45,-8,0)) ; axis3d('y--',pos=c(-45,-8,0),at=yticks) ; axis3d('z-+',pos=c(-45,8,zlimits[1]),at=zticks)
        mtext3d("water",'x--', line = 0.75, at = NULL, pos = NA) ; mtext3d("nitrogen",'y--', line = 0.5, at = NULL, pos = NA) ; mtext3d("z",'z-+', line = 1.25, at = NULL, pos = NA)
	title3d(paste0(file_stem," simulated surface"))


        writeWebGL(file=paste0(save_dir,"/",file_stem,".html"),prefix="",font="Arial",width=700,height=700)
	snapshot3d(file=paste0(save_dir,"/",file_stem,"_snapshot.png"),fmt="png",top=TRUE)
	clear3d()
        }


        
        


# for each surface, experimental or simulated, compute a set of standardized rays
# relative to the peak of the surface; find regularly spaced mesh points along those rays;
# compute the discrete differential of the relative change in z at each point.
#
# compute and save the new relative mesh points for each of the simulated and experimental surfaces
# this uses the analysis_fcns.r:find_relative_diffs_along_rays/4 function
#
# surfaces includes both experimental and simulated surfaces
#
# Kazic and Stapleton, 2.7.2017


# todo: write out to log later
# todo: make nicer loop to cover the three cases and remove cut-and-paste



# make_n_save_relative_meshes(c("b73","x1a","x1b","x2a","x2b","x3a","x3b"),"test.org",10)


# stopped here:  need to revise to use new relative meshes and their models 19.7.2017
#
# huh?
#
# but clean up the cut-and-paste
#
# Kazic, 23.7.2017


make_n_save_relative_meshes <- function(surfaces,out_file_stem,num_points) {

        timestamp <- as.integer(as.POSIXct(Sys.time()))
        results_root <- toggle_partition(0)
	out_dir  <- paste0(results_root,"/results/parameter_sweeps/ann_replots/")
        out_file <- paste0(out_dir,out_file_stem)


        x=seq(-42,42,0.5) ;  y=seq(-7.5,7.5,0.5)



        for ( s in 1:length(surfaces) ) {


                if ( "mo17" %in% surfaces) { type <- "trough" } else { type <- "domed" }


# load files in this order as the experimental surfaces were saved with the variable name "z"

                load(paste0(out_dir,"z",surfaces[s],".rda"))		
                z_exp <- z

                load(paste0(out_dir,"z_from_lm_",surfaces[s],".rda"))
                z <- get(paste0("z_from_lm_",surfaces[s]))
		
                load(paste0(out_dir,"z_np_from_lm_",surfaces[s],".rda"))
                z_np <- get(paste0("z_np_from_lm_",surfaces[s]))



                foo <- find_relative_diffs_along_rays(type,x,y,z,num_points)
                foo_np <- find_relative_diffs_along_rays(type,x,y,z_np,num_points)
                foo_exp <- find_relative_diffs_along_rays(type,x,y,z_exp,num_points) 				


		name_stem <- paste0("s_",surfaces[s],"_",num_points,"pts")



# foo*[[1]] is the matrix of (x,y,z) positions of the relative mesh points
#
# foo*[[2]] is the matrix of discrete differentials in z relative to total z
# between successive mesh points along each ray


# linear model with peak included

                xyz_mesh <- foo[[1]]
         	z_rels_diff_mat <- foo[[2]]

		xyz_name <- paste0(name_stem,"_mesh")
		assign(xyz_name,xyz_mesh)
                save(list=xyz_name,file=paste0(out_dir,xyz_name,".rda"))
		
		z_rels_name <- paste0(name_stem,"_rels")
		assign(z_rels_name,z_rels_diff_mat)
                save(list=z_rels_name,file=paste0(out_dir,name_stem,"_rels.rda"))



# linear model with peak excluded

                xyz_mesh_np <- foo_np[[1]]
		z_rels_diff_mat_np <- foo_np[[2]]

                xyz_name_np <- paste0(name_stem,"_np_mesh")
		assign(xyz_name_np,xyz_mesh_np)
                save(list=xyz_name_np,file=paste0(out_dir,name_stem,"_np_mesh.rda"))
		
		z_rels_name_np <- paste0(name_stem,"_np_rels")
		assign(z_rels_name_np,z_rels_diff_mat_np)
                save(list=z_rels_name_np,file=paste0(out_dir,name_stem,"_np_rels.rda"))



# experimental surface

                xyz_mesh_exp <- foo_exp[[1]]
		z_rels_diff_mat_exp <- foo_exp[[2]]
		
		xyz_name_exp <- paste0(name_stem,"_exp_mesh")
		assign(xyz_name_exp,xyz_mesh_exp)
                save(list=xyz_name_exp,file=paste0(out_dir,name_stem,"_exp_mesh.rda"))
		
		z_rels_name_exp <- paste0(name_stem,"_exp_rels")
		assign(z_rels_name_exp,z_rels_diff_mat_exp)
                save(list=z_rels_name_exp,file=paste0(out_dir,name_stem,"_exp_rels.rda"))
                }
        }












# ok, test how well the computed linear models fit the experimental
# surfaces.  Get the parameters from the linear fit model and use these to
# compute a simulated surface, then compare corresponding mesh points by
# computing the magnitudes and directions of the vectors joining the two
# sets of mesh points.
#
#
# uses results of simulate_surfaces_using_lm/2 below.
#
#
# Then scootch the simulated surface so its peak coincides in z with the
# target experimental surface.  The idea is to adjust so we get the same
# mesh points, since we are interested in shapes, not so much quantitative
# values.
#
# Compute the mesh points of the scootched simulated surface.
#
# Because we scootch the simulated surfaces to the z_max of each
# experimental surface, and the next contour down from each of those
# z_maxes are different, many pairs of simulated and experimental surfaces
# will not have mesh points at all of the upper contours:  36 -- 31.  So
# for those comparisons, we return NAs.
#
#
# For those mesh points, compute the magnitudes and directions of the
# vectors connecting the mesh points of the simulated surface to
# corresponding mesh points of each of the experimental surfaces.
#
# magnitude: http://www.wikihow.com/Find-the-Magnitude-of-a-Vector
#
#
# The magnitudes are the total displacement, and convey how smoothly and
# parallelly the curvatures of the simulated and experimental surfaces change
# relative to each other.  A constant magnitude --- rho --- over the entire
# set of paired mesh points indicates the two surfaces are parallel
# everywhere but the peak.  Changes in rhos over the pairs of mesh points
# indicate changes in relative curvatures in those neighborhoods.  We could
# have computed the curvatures directly, but our experience was that these
# are relatively values numerically, so a parameter with a bigger dynamic
# range seemed more informative.
#
#
# The rotation of the surfaces about the scootched, coincident peak
# relative to each other is given by theta.  The sign of theta is correct
# because we compute this by atan2/2.  Different surfaces have different
# relative rotations at different corresponding mesh points, reflecting
# changes in position of the mesh points computed using the same algorithm.
#
#
# Lastly, the sharpness of the surfaces relative to each other is given by
# phi.  If the two surfaces have the same local curvature at the
# corresponding points, then phi = 0.  Nonzero phi indicates one surface is
# sharper or shallower than the other at that point.  So this is a proxy
# for the differences in the local gradients near those points.
#
#


# we'll have two matrices to store the results, one for the lm w/peak and
# the other for the lm w/o peak
#
#
# Kazic, 2.4.2017



# hmmm, maybe just do something separate for mo17??
# yes, just separate into two lists and write out separate matrices for
# mo17 comparisons
#
# Kazic, 9.4.2017



# now, if we plan to compute rhos, phis, and thetas for MILLIONS of
# simulated surfaces, we should at least:
#
#     compute the normals for the mesh points of the experimental surfaces
#     and store these in a matrix;
#
#     for each simulated surface and its mesh points, compute the normals
#     and store those in a matrix; when we shift the surface up or down to
#     align the peaks with those of the experimental surface, we just add or
#     subtract the shift in z to the z component of the normals.
#
# BUT, we are not doing this for millions; so we'll do the second thing for
# each simulated surface, since that's how the loop goes, but we won't
# bother with the rest or an extensive restructuring of the code.
#
# Stapleton and Kazic, 25.5.2017





# ok, now for the phis, NOT on the shifted vectors (which are constant for all points,
# since z = z axis throughout, and so completely useless) BUT on the normal to the
# simulated surface at each mesh point.
#
# two ideas:
#
# 1.  compute the normal to the simulated surface at its mesh point, then
# compute the cartographic azimuth of the sim -> exp vector relative to that normal
#
# 2.  do we really want the tensor instead?
# https://en.wikipedia.org/wiki/Tensor
#
# well, not sure.  For each pair of mesh points, the set of vectors connecting them could
# form a tensor.
#
# https://www.grc.nasa.gov/www/k-12/Numbers/Math/documents/Tensors_TM2002211716.pdf
# tensors for us
#
#
# Let's find the normal at each mesh point and compute the
# cartographic azimuth relative to that:
#
# https://en.wikipedia.org/wiki/Azimuth
#
# So the values of the azimuth will go from 0 to 2pi, moving clockwise
# from the local reference vector that is the normal.  Values between pi/2
# and 3pi/2 are in the 'southern hemisphere', that is, the experimental is
# below the simulated surface.
#
#
# If 0 =< phi =< pi/2, then vector is in upper right quadrant and
# simulated is below experimental.
#
# If pi/2 < phi =< pi, then vector is in lower right quadrant and
# experimental is below simulated.
#
# If pi < phi =< 3pi/2, then vector is in lower left quadrant and
# experimental is below simulated.
#
# If 3pi/2 < phi =< 0, then vector is in upper left quadrant and simulated
# is below experimental.
#
#
# http://www.mathworks.com/matlabcentral/newsreader/view_thread/151925#849830
#
# matlab norm === R Norm/2 in package pracma, second argument is 2


# we decide to calculate phis by forming the normal to the simulated surface at the
# sim surfaces' mesh points, and compute phi relative to those normals.
#
# http://tutorial.math.lamar.edu/Classes/CalcIII/GradientVectorTangentPlane.aspx
# http://mathworld.wolfram.com/NormalVector.html
#
# partial differentials/derivatives:
# http://r.789695.n4.nabble.com/Partial-Derivatives-in-R-td886992.html
#
# http://stats.stackexchange.com/questions/70967/mathematics-behind-how-the-normal-to-the-plane-is-derived-from-a-linear-model-wh
#
# https://www.plymouth.ac.uk/uploads/production/document/path/3/3742/PlymouthUniversity_MathsandStats_gradients_and_directional_deriavatives.pdf
#
# Stapleton and Kazic, 10.5.2017


# atan2 with vectors:
#
# http://stackoverflow.com/questions/21483999/using-atan2-to-find-angle-between-two-vectors
#
# https://www.mathworks.com/matlabcentral/answers/16243-angle-between-two-vectors-in-3d
# translating to R,
#
# atan2(Norm(extprod3d(v1,v2)),crossprod(v1,v2))




# Note we do NOT want the Hausdorff distance (see https://en.wikipedia.org/wiki/Hausdorff_distance),
# because this is a summary function finding the maximum of the infinum distance between two sets.
# So the infinum can be between any two points, one from each set, and not necessarily the
# corresponding mesh points.  To do that, we would need a loop or a matrix to compute all
# the Hausdorff distances, which we're already doing.





# http://mathworld.wolfram.com/NormalVector.html
#
# N = [ f_x(x_0,y_0) , f_y(x_0,y_0), -1 ]^T
#
# gradients: https://www.plymouth.ac.uk/uploads/production/document/path/3/3742/PlymouthUniversity_MathsandStats_gradients_and_directional_deriavatives.pdf
#
# https://www.cs.purdue.edu/homes/cmh/distribution/books/chap6.pdf
# https://www.cs.purdue.edu/homes/cmh/distribution/books/
#
# https://www.mathworks.com/help/matlab/math/calculate-tangent-plane-to-surface.html?s_tid=gn_loc_drop
#
#
# what we want to do is find two vectors tangent to the simulated surface
# at each mesh point; take the cross product (physicist's sense) between them
# to find the normal; and then compute the atan2 of the angle between that normal
# and the vector between the simulated mesh point and its corresponding exptl_mesh point,
# for each simulated mesh point.  This is our vector of phis.
#
#
# To ensure we're getting the correct normal (pointing from the simulated surface to
# the experimental if sim is below exptl; if sim above exptl, then the angle of
# inclination will be negative and that's what we want), we must ensure that the angle
# between our two noncollinear tangent vectors is between 0 and pi/2.
#
# But actually, if we just take the normal having the maximum k component, that SHOULD
# do it.  And we then don't need to test for the angles between the possible pair of
# adjacent tangent vectors.
#
#
# Anyway, we have a grid of four points around each mesh point, e.g.:
#
# > z_mesh
#      [,1]         [,2]       [,3]
# [1,]   36  -0.53918307  1.2981583
#
# > which(x >= z_mesh[1,2] - 0.5 & x <= z_mesh[1,2] + 0.5)
# [1] 83 84
# Browse[2]> x[83:84]
# [1] -1.0 -0.5
# Browse[2]> which(y >= z_mesh[1,3] - 0.5 & y <= z_mesh[1,3] + 0.5)
# [1] 18 19
# Browse[2]> y[18:19]
# [1] 1.0 1.5
# Browse[2]> adj_z[83:84,18:19]
#          [,1]     [,2]
# [1,] 36.94163 35.36802
# [2,] 36.93810 35.36448
# Browse[2]> gradient(adj_z[83:84,18:19])
# $X
#           [,1]      [,2]
# [1,] -1.573616 -1.573616
# [2,] -1.573616 -1.573616
#
# $Y
#              [,1]         [,2]
# [1,] -0.003536301 -0.003536301
# [2,] -0.003536301 -0.003536301
#
#
#
# Note we don't want a bigger patch of surface around the mesh point, as the
# gradients start changing rapidly, meaning our tangent vectors might
# not be tangent and instead go through the surface.



# Now one of the most important things to do BEFORE we compute the normal is to
# recognize that we're in a local coordinate system with z=constant (the contour value).  So
# when we make the grid points around the mesh point, their k components are all 0, NOT the
# value of the contour line.  Otherwise, the i and j components of the cross product will be
# nonzero and be really strange --- because in that case the vectors are all defined relative
# to (0,0,0) instead of (mesh_x,mesh_y,0).
#
# Remember the tangent vectors to the mesh point all lie in the same z plane, so no
# projections are needed.
#
#
#
# If the mesh point lies exactly in a row or column, the angle between the tangent
# vectors will be pi, and the normal will vanish.  If the angle between the tangents is
# between 0 and pi, the normal points up; if between pi and 2pi, it points down.  See the
# nice animation in
# https://en.wikipedia.org/wiki/Cross_product.
#
# The maximum magnitude of the normal will be at pi/2 and 3pi/4.  For numerical reasons,
# it's probably a good idea to have the angle between the tangents be as close to pi/2
# as possible.
#
#
# For example, to find this:
#
# mp <- z_mesh[1,2:3]; mp
# [1] -0.5391831  1.2981583
#
# For this case, we want the points {(-1,1),(-0.5,1),(-1,1.5),(-0.5,1.5)}, which we can
# get with expand.grid:
#
# grid_pts <- expand.grid(x[83:84],y[18:19])
#   Var1 Var2
# 1 -1.0  1.0
# 2 -0.5  1.0
# 3 -1.0  1.5
# 4 -0.5  1.5
#
# Var1 is x and Var2 is y here.
#
# Translating each vector to an origin at mp:
# > grid_pts - mp
#         Var1       Var2
# 1 -0.4608169  1.5391831
# 2 -1.7981583 -0.2981583
# 3 -0.4608169  2.0391831
# 4 -1.7981583  0.2018417

# > trans <- grid_pts - mp

# Note the R operator for the dot product of two vectors is called crossprod!
# crossprod(as.vector(as.numeric(trans[1,])),as.vector(as.numeric(trans[2,])))
#           [,1]
# [1,] 0.3697015



# need to write stuff out to log file . . . 






# from the expt mesh to the sim mesh; phi is the angle between that test_vec
# and the normal to the sim surface at its mesh point

#                test_vec <- c((exp_mesh[i,1:2] - sim_mesh[i,2:3]),sim_mesh[i,1])

# this version is the difference vector between the experimental mesh point and the
# simulated mesh point; z = 0 as they are both on the same contour

#                test_vec <- exp_mesh[i,] - c(sim_mesh[i,2:3],sim_mesh[i,1])



# when surfaces are severely rotated --- e.g., x3b in either model --- what
# to do?
#
# Kazic, 2.6.2017





# hmmm do we need to rescale the simulated surfaces to match those
# of the experimental surfaces?
#
# for example,
#

# max exp - min exp
# --------------------------    x  value sim z
# adj max sim - adj min sim

# 40
#----  x 50 = 10
# 200


# 40
# -- x 10 = 10
# 40 

# but do we want to rescale?  Maybe just compare the gradients?

# since we're trying to compare shapes, then scale should be the same . . .
# standardize sim to experimental range --- if just shape, then numerical difference
# doesn't matter and we can use our theta, phi after adjusting peak and rescaling
#
# and if not just shape, then the numerical difference matters; and so maybe this is rho, and
# maybe we don't even adjust the peak; or maybe the gradients


# and we also need to pick the mesh points in the lower left quadrant when there are
# symmetric points in the upper right quadrant.  Since we know where the peak is, we
# discard mesh points that are above and to the right of the peak, for B73-likeish things.


# and we still need to figure out what to do with x3b




# so the mesh and discrete differential matrices are organized by ray:
# mesh is ray:x, ray:y, ray:z; and differential is ray0, ray1, ...
#
# For any pair of surfaces, the position of the peak determines how far each ray
# extends relative to the experimental surface.
#
# So we want to compute, for all rays:
#  the Euclidean distances between each pair of points on a given ray;
#  the angle, theta, between pair of points on a given ray (maybe just theta by rays?);
#  and the difference matrix of the two differential matrices.
#
#
# output is vectorized over all rays, one row/pair of compared surfaces, one matrix for each type of datum
# then plot a heatmap for each matrix
#
# Stapleton and Kazic, 2.7.2017


# modified to use rel_mesh variable and file names
#
# Kazic, 19.7.2017


check_lm_fit_quality <- function(exptl_surfaces,sim_surfaces,out_file_stem,num_points) {

        timestamp <- as.integer(as.POSIXct(Sys.time()))
        results_root <- toggle_partition(0)
	out_dir  <- paste0(results_root,"/results/parameter_sweeps/ann_replots/")
        out_file <- paste0(out_dir,out_file_stem)


        rhos <- vector(length=0)
        thetas <- vector(length=0)
        diffntls <- vector(length=0)


        rhos_np <- vector(length=0)
        thetas_np <- vector(length=0)
        diffntls_np <- vector(length=0)

        row_names <- vector(length=0)
        row_names_np <- vector(length=0)	


        x=seq(-42,42,0.5) ;  y=seq(-7.5,7.5,0.5)

        num_rays <- 6




# no scootching, which would give a distance relative to the experimental
# surface


        for ( e in 1:length(exptl_surfaces) ) {

                load(paste0(out_dir,"s_",exptl_surfaces[e],"_",num_points,"pts_exp_mesh.rda"))
                load(paste0(out_dir,"s_",exptl_surfaces[e],"_",num_points,"pts_exp_rels.rda"))

                exp_mesh <- get(paste0("s_",exptl_surfaces[e],"_",num_points,"pts_exp_mesh"))
                exp_rels <- get(paste0("s_",exptl_surfaces[e],"_",num_points,"pts_exp_rels"))


                exp_mesh_stack <- stack_matrix(exp_mesh)



                for ( s in 1:length(sim_surfaces) ) {


# load the surfaces simulated using the two versions of the linear model

                        load(paste0(out_dir,"s_",sim_surfaces[s],"_",num_points,"pts_mesh.rda"))
                        load(paste0(out_dir,"s_",sim_surfaces[s],"_",num_points,"pts_rels.rda"))		
                        z_mesh <- get(paste0("s_",sim_surfaces[s],"_",num_points,"pts_mesh"))
                        z_rels <- get(paste0("s_",sim_surfaces[s],"_",num_points,"pts_rels"))
                        z_mesh_stack <- stack_matrix(z_mesh)


                        load(paste0(out_dir,"s_",sim_surfaces[s],"_",num_points,"pts_np_mesh.rda"))
                        load(paste0(out_dir,"s_",sim_surfaces[s],"_",num_points,"pts_np_rels.rda"))		
                        z_mesh_np <- get(paste0("s_",sim_surfaces[s],"_",num_points,"pts_np_mesh"))
                        z_rels_np <- get(paste0("s_",sim_surfaces[s],"_",num_points,"pts_np_rels"))
                        z_mesh_np_stack <- stack_matrix(z_mesh_np)




# compute stuff and write out results

                        local_rhos <- sqrt(rowSums(exp_mesh_stack - z_mesh_stack)^2)
                        local_rhos_np <- sqrt(rowSums(exp_mesh_stack - z_mesh_np_stack)^2)




# now define theta as the angle between the projections of the vectors for each pair of points
# for experimental and simulated.  These projections will be in the xy plane orthogonal to the z-axis.
#
# The projections into xy are used because we are interested only in the relative rotation
# of the simulated surface compared to the experimental.
#
# All we need to consider is the first and second components of the vectors.
# projection of 3d vector onto xy plane, just drop the last coordinate:
# https://stackoverflow.com/questions/23472048/projecting-3d-points-to-2d-plane
# confirmed with wbw.
#
#
# We DO NOT translate these vectors away from the origin of the plane.  Then (as before), 
#
# theta = atan(y/x), but better, atan2(y/x), atan2(y,x) (angle in xy plane; azimuthal angle)
#
# https://stackoverflow.com/questions/1897704/angle-between-two-vectors-in-r, atan2 answer
# http://www.euclideanspace.com/maths/algebra/vectors/angleBetween/



                        local_thetas <- atan2(exp_mesh_stack[,2],exp_mesh_stack[,1]) - atan2(z_mesh_stack[,2],z_mesh_stack[,1])
                        local_thetas_np <- atan2(exp_mesh_stack[,2],exp_mesh_stack[,1]) - atan2(z_mesh_np_stack[,2],z_mesh_np_stack[,1])			

                        local_diffs <- exp_rels - z_rels
                        local_diffs_np <- exp_rels - z_rels_np




                        rhos <- c(rhos,local_rhos)
			rhos_np <- c(rhos_np,local_rhos_np)

                        thetas <- c(thetas,local_thetas)
			thetas_np <- c(thetas_np,local_thetas_np)

                        diffntls <- c(diffntls,c(local_diffs))
			diffntls_np <- c(diffntls_np,c(local_diffs_np))

                        row_names <- c(row_names,paste0("e_",exptl_surfaces[e],"::s_",sim_surfaces[s]))        
                        }
                }




# reshape the output data into six nice matrices and save them

        num_rows <- length(row_names)
        row_names_np <- paste0(row_names,"_np")
	
	num_cols <- as.integer(length(rhos)/num_rows)
	ray_names <- paste0("r",seq(0,num_rays-1,1),":")
	pt_names <- outer(ray_names,seq(1,num_points,1), paste0)
	col_names <- c(t(pt_names))

        diffn <- c(t(pt_names[,2:num_points]))
	diffn_names <- paste0("delta ",diffn)
        diffn_cols <- length(diffn_names)



        rho_matrix <- matrix(rhos,num_rows,num_cols,byrow=TRUE)
	rownames(rho_matrix) <- row_names
        colnames(rho_matrix) <- col_names

        rho_matrix_np <- matrix(rhos_np,num_rows,num_cols,byrow=TRUE)
	rownames(rho_matrix_np) <- row_names_np
        colnames(rho_matrix_np) <- col_names



        theta_matrix <- matrix(thetas,num_rows,num_cols,byrow=TRUE)
	rownames(theta_matrix) <- row_names
        colnames(theta_matrix) <- col_names

        theta_matrix_np <- matrix(thetas_np,num_rows,num_cols,byrow=TRUE)
	rownames(theta_matrix_np) <- row_names_np
        colnames(theta_matrix_np) <- col_names




        diffntls_matrix <- matrix(diffntls,num_rows,diffn_cols,byrow=TRUE)
	rownames(diffntls_matrix) <- row_names
        colnames(diffntls_matrix) <- diffn_names

        diffntls_matrix_np <- matrix(diffntls_np,num_rows,diffn_cols,byrow=TRUE)
	rownames(diffntls_matrix_np) <- row_names_np
        colnames(diffntls_matrix_np) <- diffn_names


        save(rho_matrix,file=paste0(out_dir,"rho_matrix.rda"))
        save(rho_matrix_np,file=paste0(out_dir,"rho_matrix_np.rda"))
        save(theta_matrix,file=paste0(out_dir,"theta_matrix.rda"))
        save(theta_matrix_np,file=paste0(out_dir,"theta_matrix_np.rda"))
        save(diffntls_matrix,file=paste0(out_dir,"diffntls_matrix.rda"))
        save(diffntls_matrix_np,file=paste0(out_dir,"diffntls_matrix_np.rda"))
        }










# uh, compare the experimental surfaces to each other in the same way as a control!
#
#
# remember, the meshes are of relative points, and the "_rel" is the discrete differentials in z
# along each ray, relative to the total z of that ray
#
# Stapleton and Kazic, 23.7.2017


check_experimental_surfs_vs_each_other <- function(exptl_surfaces,out_file_stem,num_points) {

        timestamp <- as.integer(as.POSIXct(Sys.time()))
        results_root <- toggle_partition(0)
	out_dir  <- paste0(results_root,"/results/parameter_sweeps/ann_replots/")
        out_file <- paste0(out_dir,out_file_stem)


        rhos <- vector(length=0)
        thetas <- vector(length=0)
        diffntls <- vector(length=0)
        row_names <- vector(length=0)


        x=seq(-42,42,0.5) ;  y=seq(-7.5,7.5,0.5)

        num_rays <- 6




# no scootching, which would give a distance relative to the experimental
# surface


        for ( e in 1:length(exptl_surfaces) ) {

                load(paste0(out_dir,"s_",exptl_surfaces[e],"_",num_points,"pts_exp_mesh.rda"))
                load(paste0(out_dir,"s_",exptl_surfaces[e],"_",num_points,"pts_exp_rels.rda"))

                exp_mesh <- get(paste0("s_",exptl_surfaces[e],"_",num_points,"pts_exp_mesh"))
                exp_rels <- get(paste0("s_",exptl_surfaces[e],"_",num_points,"pts_exp_rels"))
                exp_mesh_stack <- stack_matrix(exp_mesh)



                for ( s in 1:length(exptl_surfaces) ) {

                        load(paste0(out_dir,"s_",exptl_surfaces[s],"_",num_points,"pts_exp_mesh.rda"))
                        load(paste0(out_dir,"s_",exptl_surfaces[s],"_",num_points,"pts_exp_rels.rda"))


# compare each surface to the others
			
                        z_mesh <- get(paste0("s_",exptl_surfaces[s],"_",num_points,"pts_exp_mesh"))
                        z_rels <- get(paste0("s_",exptl_surfaces[s],"_",num_points,"pts_exp_rels"))
                        z_mesh_stack <- stack_matrix(z_mesh)





# compute stuff and write out results, same as in
# check_lm_fit_quality/4

                        local_rhos <- sqrt(rowSums(exp_mesh_stack - z_mesh_stack)^2)
                        local_thetas <- atan2(exp_mesh_stack[,2],exp_mesh_stack[,1]) - atan2(z_mesh_stack[,2],z_mesh_stack[,1])
                        local_diffs <- exp_rels - z_rels



                        rhos <- c(rhos,local_rhos)
                        thetas <- c(thetas,local_thetas)
                        diffntls <- c(diffntls,c(local_diffs))
                        row_names <- c(row_names,paste0("e_",exptl_surfaces[e],"::e_",exptl_surfaces[s]))        
                        }
                }




# reshape the output data into three nice matrices and save them

        num_rows <- length(row_names)

	
	num_cols <- as.integer(length(rhos)/num_rows)
	ray_names <- paste0("r",seq(0,num_rays-1,1),":")
	pt_names <- outer(ray_names,seq(1,num_points,1), paste0)
	col_names <- c(t(pt_names))

        diffn <- c(t(pt_names[,2:num_points]))
	diffn_names <- paste0("delta ",diffn)
        diffn_cols <- length(diffn_names)



        rho_matrix_e_by_e <- matrix(rhos,num_rows,num_cols,byrow=TRUE)
	rownames(rho_matrix_e_by_e) <- row_names
        colnames(rho_matrix_e_by_e) <- col_names


        theta_matrix_e_by_e <- matrix(thetas,num_rows,num_cols,byrow=TRUE)
	rownames(theta_matrix_e_by_e) <- row_names
        colnames(theta_matrix_e_by_e) <- col_names


        diffntls_matrix_e_by_e <- matrix(diffntls,num_rows,diffn_cols,byrow=TRUE)
	rownames(diffntls_matrix_e_by_e) <- row_names
        colnames(diffntls_matrix_e_by_e) <- diffn_names

        save(rho_matrix_e_by_e,file=paste0(out_dir,"rho_matrix_e_by_e.rda"))
        save(theta_matrix_e_by_e,file=paste0(out_dir,"theta_matrix_e_by_e.rda"))
        save(diffntls_matrix_e_by_e,file=paste0(out_dir,"diffntls_matrix_e_by_e.rda"))
        }























# given a matrix of (x,y,z) values for each of n rays in 3n columns, stack the rays together into
# a 3 column matrix of (x,y,z) for all rays

stack_matrix <- function(amatrix) {
        cols <- ncol(amatrix)
	rows <- nrow(amatrix)	
	stack_rows <- as.integer(rows*cols/3)
        stacked <- matrix(NA,stack_rows,3)


	counter <- 1
	next_row <- 1

	for ( i in seq(1,cols,3) ) {
	        last_row <- counter*rows
		last_col <- i + 2
	        stacked[next_row:last_row,] <- amatrix[,i:last_col]
		next_row <- last_row + 1
		counter <- counter + 1
		}

        return(stacked)
        }











## # need to finish the mo17 case . . . 
## #
## #                                mo17_ray0 <- get_mo17_ray0(x,y,z)
## #                                mo17_rays <- get_mo17_rays(c(-42,7.5),mo17_ray0$coefficients[2],mo17_ray0$coefficients[1],-30)
## #                                get_mo17_intersectns(mo17_rays,.....
## #                         if ( exptl_surfaces[e] == "mo17" ) { contour_set <- "mo" } else { contour_set <- "sim_domed" }







## ## # now figure out how to translate z_mesh to the origin, and compute theta in spherical coordinates
## ## # on the translated vectors
## ## #
## ## # well, translation is just sim - exptl for the (x,y,z) columns:
## ## # http://www.mathplanet.com/education/geometry/transformations/vectors 

## ##                         zs <- as.numeric(exptl_meshes[top:rows,"z"])
## ## 			rows_shifted <- rows - top + 1
## ##                         zs_np <- as.numeric(exptl_meshes[top_np:rows,"z"])
## ## 			rows_shifted_np <- rows - top_np + 1
			
## ##                         shifted_vectors <- matrix(c(zs,exptl_meshes[top:rows,c(exptl_xys[1],exptl_xys[2])]),rows_shifted,3) - z_mesh
## ##                         shifted_vectors_np <- matrix(c(zs_np,exptl_meshes[top_np:rows,c(exptl_xys[1],exptl_xys[2])]),rows_shifted_np,3) - z_np_mesh


## ## 


## ## 			thetas <- atan2(shifted_vectors[,3],shifted_vectors[,2])
## ## 			thetas_np <- atan2(shifted_vectors_np[,3],shifted_vectors_np[,2])


## ## # angles of inclination between normals at simulated and experimental mesh points
## ## #
## ## # actually, calculate this all at once for each simulated surface, just like for the
## ## # experimental surfaces, and then just repeat the computations with those matrices
## ## #

## ## # compute the angle between the normals at the simulated mesh point and the experimental
## ## # mesh point, for each such pair.  If the two surfaces are coincident there,
## ## # then the angle is 0.


## ##                         sim_normals <- give_us_normals(x,y,xs,ys,z_mesh,adj_z)
## ##                         sim_normals_np <- give_us_normals(x,y,xs,ys,z_np_mesh,adj_z_np)



## ## # then compute the angles given these normals and the global set of experimental normals found above
## ## # ,exptl_meshes[top_np:rows,c(exptl_xys,"z")])










## # for each of the surfaces, load its simulations (w/ and w/o peaks), and
## # then for each, compare it to the derived data from each of the
## # experimental surfaces: compute vectors, magnitudes,
## # and directions between the mesh points of the experimental surfaces
## # (always the local origin of the vectors) and the mesh points of the
## # simulated surfaces (the tips of the vectors).


##         for ( s in 1:length(sim_surfaces) ) {


## # load the surfaces simulated using the two versions of the linear model

##                 load(paste0(results_root,"/results/parameter_sweeps/ann_replots/z_from_lm_",sim_surfaces[s],".rda"))
##                 load(paste0(results_root,"/results/parameter_sweeps/ann_replots/z_np_from_lm_",sim_surfaces[s],".rda"))


##                 z <- get(paste0("z_from_lm_",sim_surfaces[s]))
##                 z_np <- get(paste0("z_np_from_lm_",sim_surfaces[s]))

##                 z_peak <- get_peak(z)
## 		z_np_peak <- get_peak(z_np)


##                 foo <- find_relative_diffs_along_rays(x,y,xs,ys,z,num_points)
##                 foo_np <- find_relative_diffs_along_rays(x,y,xs,ys,z_np,num_points) 		


##                 xyz_mesh <- foo[[1]]
##          	z_rels_diff_mat <- foo[[2]]
		 

##                 xyz_mesh_np <- foo_np[[1]]
## 		z_rels_diff_mat_np <- foo_np[[2]]
		 


## # three coordinates/mesh point, for all 65 possible mesh points
## # ignore the mo17 issues for now


##         ## cols <- 3 * nrow(exptl_meshes)
## 	## rows_z_matches <- length(sim_surfaces) * length(exptl_surfaces)
##         ## z_matches_mat <- matrix(z_matches,rows_z_matches,cols,byrow=TRUE)
##         ## z_np_matches_mat <- matrix(z_matches_np,rows_z_matches,cols,byrow=TRUE)

##         ## tmp <- as.vector(t(outer("s_",sim_surfaces,paste0)))
## 	## first <- paste0(tmp,":")
##         ## rownames(z_matches_mat) <- as.vector(t(outer(first,exptl_surfaces,paste0)))
## 	## rownames(z_np_matches_mat) <- rownames(z_matches_mat)


##         ## colnames(z_matches_mat) <- as.vector(t(outer(rownames(exptl_meshes),c(":dirho",":phi",":theta"),paste0)))
##         ## colnames(z_np_matches_mat) <- colnames(z_matches_mat)


##         ## analysis_dir <- paste0(results_root,"/results/parameter_sweeps/analysis/")
##         ## save(z_matches_mat,file=paste0(analysis_dir,"z_matches_mat_",timestamp,".rda"))
## 	## save(z_np_matches_mat,file=paste0(analysis_dir,"z_np_matches_mat_",timestamp,".rda"))

##         }







# if a mesh point is outside the evaluation plane, set it to NA

rm_out_of_plane_points <- function(lower_x,upper_x,lower_y,upper_y,mesh) {
        outies <- which(mesh[,2] < lower_x | mesh[,2] > upper_x | mesh[,3] < lower_y | mesh[,3] > upper_y)
        mesh[outies,2:3] <- NA
	return(mesh)
        }







# find suitable normals around each point in a matrix of mesh points
# by "suitable", we mean normals that point AWAY from the surface, meaning
# upward for a convex (domed) surface.


give_us_normals <- function(x,y,xs,ys,mesh,surface) {

#        print(paste0("dim mesh is: ",dim(mesh)))


         normals <- vector(length=0)


# just grab a rectangular grid of points around the mesh point, but find
# their zs by interpolation on the surface once they're moved to the location of
# the mesh point.  The reason for this is that if the grid points are at a very
# different z than the mesh point, then the normal will point in a different direction,
# and that of course is what we want to find!
#
# So we can't be sloppy and assume we have a local plane that is perpendicular to the
# z axis.
#
# Stapleton and Kazic, 25.5.2017


# two cases can arise.
#
#    1.  the mesh point is near an edge, so that incrementing an arbitrary
#        amount takes us outside the evaluation (x,y) plane.  Then interpolation
#        fails.
#
#    2.  the mesh point is outside the evaluation plane, so interpolation will fail.
#        This case is avoided by first substituting NAs for those points, in
#        rm_out_of_plane_points/3.
#



        for ( i in 1:nrow(mesh) ) {

#                print(paste0("i: ",i))


# figure out what to do when the mesh point is near, but not over, the edge.
# If the mesh point is too close to an edge, just get the corresponding grid points
# from the evaluation plane.
#
# stopped here --- obsolete?



                 wee_xs <- c(mesh[i,2] - 0.5, mesh[i,2] + 0.5)
		 wee_ys <- c(mesh[i,3] - 0.5, mesh[i,3] + 0.5)
                 surf_pts <- expand.grid(wee_xs,wee_ys)
                 interpolated_zs <- interpp(xs,ys,c(surface),surf_pts[,1],surf_pts[,2],linear=TRUE,extrap=FALSE,duplicate="error")
                 xyz_surf_pts <- matrix(c(surf_pts,interpolated_zs),4,3,byrow=FALSE)             





#		grid_pts <- cbind(surf_pts,mesh[i,1])
#                print(grid_pts)


# and, we want the correct z as we will have moved the normal back to the mesh point


                normal <- pick_normal(grid_pts,c(mesh[i,2:3],mesh[i,1]))
                print(paste0("i: ",i," normal: ",normal))		
                normals <- c(normals,normal)
                }


        normals <- matrix(normals,nrow(mesh),3,byrow=TRUE)
        return(normals)
        }
















# given a 4 x 2 matrix of grid points around a mesh point, pick the two
# adjacent tangent vectors that we'll use to calculate the normal to the
# surface at the mesh point.  By convention, we want this normal
# to point away from the surface towards +z.
#
# We subtract the vectors to get them rooted at the mesh point,
# but be sure we do it in the right order to get the
# normals going the right direction --- parallelogram rule.
#
# Kazic, 18.5.2017




pick_normal <- function(grid_pts,mesh_pt) {


# angles is a vector of angles in radians, such that angles[i] is between trans_aug[i] and trans_aug[i+1]

        angles <- vector(length=0)


# We want the normal to point away from the surface, not towards its
# interior.  All the points, and therefore their Euclidea vectors,
# are defined relative to the origin of the evaluation
# space at (0,0,0).  So we must first compute the vectors between the mesh point
# and each of the grid points.  To do this, we subtract the points in the order
#
# grid point - mesh point
# if the grid point is to the left of the mesh point, and
#
# mesh point - grid point
# if the grid point is to the right of the mesh point.
#
# We determine left and right by comparing the x-coordinates of the points.
#
# See https://en.wikipedia.org/wiki/Euclidean_vector (section "Addition and Subtraction") and
# for generalizations, https://en.wikipedia.org/wiki/Vector_space.


	trans_aug <- rbind(grid_pts,grid_pts[1,])
        dots <- diag((as.matrix(trans_aug) %*% as.matrix(t(trans_aug)))[2:5,1:4])
        angles <- acos(dots/(sqrt(sum(as.matrix(trans_aug[1:4,])^2)) * sqrt(sum(as.matrix(trans_aug[2:5,])^2))))




        if ( is.na(angles) == TRUE ) {
                normal <- c(NA,NA,NA)
		print("warning: set normal to NA!")
                } else {

                        max_angle_index <- which(angles == max(angles) & max(angles) <= pi/2)
#			print(paste0(" max_angle_index: ",max_angle_index))
                        if ( ( max_angle_index >= 1 ) && ( max_angle_index <= 5 ) ) {
                                next_index <- max_angle_index + 1
				
		                
# extprod3d is the cross product we know from physics, defined in https://en.wikipedia.org/wiki/Cross_product
#
# stopped here! double-check this next!  We want to compute the normal at the mesh point, not at the origin,
# then move it back to the mesh point

                                local_normal <- extprod3d(as.matrix(trans_aug[max_angle_index,] - mesh_pt),as.matrix(trans_aug[next_index,] - mesh_pt))
                		normal <- local_normal + mesh_pt
                                } else {
                        	        normal <- c(NA,NA,NA)
		                	print("warning: set normal to NA!")
                                        }
			        }
			
	return(normal)
        }













# gotta plot:  one for each of rho, phi, theta
#
# need to color ranges consistently across types of data:
#
#    rho: 0 -> 285      (quite different)
#  theta: -3.0 -> 6.25   (already ok)
#  diffntls -0.19xxxx -> 0.26
#
# see https://www.rdocumentation.org/packages/plotrix/versions/3.6-4/topics/color2D.matplot
# and http://www.statmethods.net/advgraphs/parameters.html
#
#
# for the experimental surfaces against each other
#
#    rho: 0 -> 80      (quite different)
#  theta: -6.5 -> 6.5 
#  diffntls -0.2 -> 0.2
#
#


# color by quantiles?
#
# use Rebecca Barter's superheat package, quite lovely
# http://rlbarter.github.io/superheat/index.html
#
# DON'T use scale=TRUE, or too many values below the maximum disappear
#
# Kazic, 19.7.2017


plot_heatmaps <- function(dir) {

       params <- c("rho","theta","diffntls")
       ranges <- matrix(c(0,300,-6.5,6.5,-0.2,0.3),3,2,byrow=TRUE)
       rownames(ranges) <- params

       for ( i in params ) {


# for just redoing the diffntls
#
# Kazic, 201.11.2017
#
#              i <- "diffntls"

              plot_limits <- ranges[i,]

              mat_name <- paste0(i,"_matrix")
              mat_np_name <- paste0(i,"_matrix_np")
       	
	      load(paste0(dir,"/",mat_name,".rda"))
	      load(paste0(dir,"/",mat_np_name,".rda"))
	
       	      mat <- get(mat_name)
       	      mat_np <- get(mat_np_name)

	      
	      mat <- reorder_by_classificatn_category(0,mat)
 	      mat_np <- reorder_by_classificatn_category(1,mat_np)



#    reverse the order of the rows, so that b73 is on top;

	      mat <- mat[rev(rownames(mat)),]
	      mat_np <- mat_np[rev(rownames(mat_np)),]



# order columns by points, not rays; clarifies true pattern
#
# fixed to reorder diffntls correctly
#
# Kazic, 20.11.2017

              if ( i != "diffntls") { 
                      re_mat <- reorder_matrix_cols_for_superheat(mat)
                      re_mat_np <- reorder_matrix_cols_for_superheat(mat_np)	      
                      } else {
                              re_mat <- reorder_diffntls_matrix_cols_for_superheat(mat)
                              re_mat_np <- reorder_diffntls_matrix_cols_for_superheat(mat_np)	      
                              }


#             cat(i," ",plot_limits,"\n")
#             cat("range of ",mat_name," is ",range(mat),"; NAs are ",any(is.na(mat)),"\n")
#             cat("quantiles of ",mat_name," are: ",quantile(sort(mat),seq(0,1,0.1),na.rm=FALSE,names=FALSE,type=7),"\n")
#             cat("range of ",mat_np_name," is ",range(mat_np),"; NAs are ",any(is.na(mat_np)),"\n")
#             cat("quantiles of ",mat_np_name," are: ",quantile(sort(mat_np),seq(0,1,0.1),na.rm=FALSE,names=FALSE,type=7),"\n")


              plot_superheats(dir,mat_name,i,"peak",plot_limits,mat)
              plot_superheats(dir,mat_np_name,i,"",plot_limits,mat_np)	      

              plot_superheats(dir,paste0(mat_name,"_reordered"),i,"peak",plot_limits,re_mat)
              plot_superheats(dir,paste0(mat_np_name,"_reordered"),i,"",plot_limits,re_mat_np)	      
              }
        }











# reorder the rows so that they follow the table
# lda_categories in ../../first_manuscript/working_draft/stuff_for_figs.org
#
# Stapleton and Kazic, 23.7.2017


reorder_by_classificatn_category <- function(np_flag,amatrix) {

        foo <- rownames(amatrix)
	new_order <- c("e_b73::s_b73",
"e_b73::s_x1a",
"e_b73::s_x2b",
"e_b73::s_x1b",
"e_b73::s_x2a",
"e_b73::s_x3a",
"e_x1a::s_x1a",
"e_x1a::s_b73",
"e_x1a::s_x2b",
"e_x1a::s_x1b",
"e_x1a::s_x2a",
"e_x1a::s_x3a",
"e_x2b::s_x2b",
"e_x2b::s_b73",
"e_x2b::s_x1a",
"e_x2b::s_x1b",
"e_x2b::s_x2a",
"e_x2b::s_x3a",
"e_x3b::s_b73",
"e_x3b::s_x1a",
"e_x3b::s_x1b",
"e_x3b::s_x2a",
"e_x3b::s_x2b",
"e_x3b::s_x3a",
"e_x1b::s_x1b",
"e_x1b::s_x2a",
"e_x1b::s_x3a",
"e_x1b::s_b73",
"e_x1b::s_x1a",
"e_x1b::s_x2b",
"e_x2a::s_x2a",
"e_x2a::s_x1b",
"e_x2a::s_x3a",
"e_x2a::s_b73",
"e_x2a::s_x1a",
"e_x2a::s_x2b",
"e_x3a::s_x3a",
"e_x3a::s_x1b",
"e_x3a::s_x2a",
"e_x3a::s_b73",
"e_x3a::s_x1a",
"e_x3a::s_x2b",
"e_mo17::s_x1b",
"e_mo17::s_x2a",
"e_mo17::s_x3a",
"e_mo17::s_b73",
"e_mo17::s_x1a",
"e_mo17::s_x2b")


        if ( as.logical(np_flag) == TRUE ) { new_order <- paste0(new_order,"_np") }

        amatrix <- amatrix[new_order,]
        amatrix <- rename_rows(amatrix)

        return(amatrix)
        }




# rename the x#{A,B} to the \qtlish names; but call these
# Q#-B for the B73 alleles and Q#-M for the Mo17 alleles.
#
# This is only for publication, of course.


rename_rows <- function(amatrix) {

        new_names <- vector(length=0)
        old_names <- rownames(amatrix)

        for ( i in 1:length(old_names) ) {
                first <- gsub("(x)","QTL",old_names[i],perl=TRUE)
		second <- gsub("(a::)","-B73::",first,perl=TRUE)
		third <- gsub("(a$)","-B73",second,perl=TRUE)
		fourth <-  gsub("(b::)","-Mo17::",third,perl=TRUE)
		fifth <-  gsub("(b$)","-Mo17",fourth,perl=TRUE)
		sixth <-  gsub("(b73)","B73",fifth,perl=TRUE)
		seventh <-  gsub("(mo17)","Mo17",sixth,perl=TRUE)				

#		cat("start :",old_names[i],": new :",seventh,":\n")

                new_names <- c(new_names,seventh)
                }

        rownames(amatrix) <- new_names
        return(amatrix)

        }









# reorder the rho, theta, and diffntls matrices:
#
#    permute the order of the columns, so that points nearer the peak are
#    closer to the top
#
# Kazic, 20.7.2017

# ah, this is fine EXCEPT for the differentials; order comes out screwy
#
# Kazic, 20.11.2017

reorder_matrix_cols_for_superheat <- function(amatrix) {
        cols <- colnames(amatrix)
	new_col_order <- c(matrix(matrix(cols,6,10,byrow=TRUE),1,60,byrow=FALSE))
        reordered <- amatrix[,new_col_order]

        return(reordered)
        }





# modified to fix order for diffntls_matrix
#
# Kazic, 20.11.2017

reorder_diffntls_matrix_cols_for_superheat <- function(amatrix) {
        cols <- colnames(amatrix)
#	new_col_order <- c(matrix(matrix(cols,6,10,byrow=TRUE),1,60,byrow=FALSE))
	new_col_order <- c(outer(paste0("delta r",seq(0,5,1)),paste0(":",seq(2,10,1)),FUN="paste0"))
        reordered <- amatrix[,new_col_order]

        return(reordered)
        }













# plot the heatmaps using superheat for prettiness

plot_superheats <- function(dir,mat_name,var_name,peak_flag,plot_limits,amatrix) {

       if ( peak_flag == "peak" ) { suffix <- "w/ peak" } else { suffix <- "w/o peak" }

       png(paste0(dir,"/",mat_name,"_heatmap.png"), height = 900, width = 800)

#       superheat(amatrix,heat.lim=plot_limits,heat.na.col="white",left.label.size = 0.125,left.label.text.size = 3,grid.hline = FALSE,grid.vline = FALSE,bottom.label.text.angle=90,bottom.label.size=0.1,bottom.label.text.size=3,title=paste0("comparison of experimental to simulated surfaces for ",var_name,", linear model ",suffix))

       superheat(amatrix,heat.lim=plot_limits,heat.na.col="white",left.label.size = 0.185,left.label.text.size = 3,left.label.text.alignment = "left",grid.hline = FALSE,grid.vline = FALSE,bottom.label.text.angle=90,bottom.label.size=0.075,bottom.label.text.size=3,title=paste0("comparison of experimental to simulated surfaces for ",var_name,", linear model ",suffix))

       dev.off()
       }














# Compute the simulated surfaces by solving the matrix equation with and
# without the peak.  The idea is that the peak is an outlier that distorts
# the solution.
#
# Matrix equation:
#
#         X^{-1} %*% p = Z',
#
# where X is the *_lm[,1:4], the values of [X^2 Y^2 X Y]; p is the vector
# of fitted parameter values obtained from either limSolve, Solve, or just
# solving:
#
#         X^{-1} %*% Z = p  (where Z is *_lm[,5]);
#
# and Z' is the expected value of surface at the [X Y] points.
#
# Note that the first row in every *_lm is the peak.  The inverses are the
# Moore-Penrose generalized pseudo-inverses.
#
# Kazic, 2.4.2017

# results_root <- toggle_partition(0); surfaces <- c("b73","mo17","x1a","x1b","x2a","x2b","x3a","x3b") ; simulate_surfaces_using_lm(surfaces,results_root)




# modified variable and file names to use linear models calculated from new relative mesh points
#
# Kazic, 19.7.2017

simulate_surfaces_using_lm <- function(surfaces,results_root) {

        x=seq(-42,42,0.5) ;  y=seq(-7.5,7.5,0.5)


        for  ( k in 1:length(surfaces) ) {	
                if ( surfaces[k] == "mo17" ) { c <- 1 } else { c <- -1 }

                z <- matrix(NA,nrow=length(x),ncol=length(y))
                z_np <- matrix(NA,nrow=length(x),ncol=length(y))


#                load(paste0(results_root,"/results/parameter_sweeps/ann_replots/",surfaces[k],"_lm.rda"))
#                x_mat <- get(paste0(surfaces[k],"_lm"))

                load(paste0(results_root,"/results/parameter_sweeps/ann_replots/",surfaces[k],"_rel_mesh_lm.rda"))
                x_mat <- get(paste0(surfaces[k],"_rel_mesh_lm"))
		b <- nrow(x_mat)
		


# solve for p, including the peak of the experimental surface

		x_inv <- ginv(x_mat[,1:4])
		p <- x_inv %*% x_mat[,5]


# solve for p_np, excluding the peak of the experimental surface

		x_inv_np <- ginv(x_mat[2:b,1:4])
		p_np <- x_inv_np %*% x_mat[2:b,5]




# compute sim surfaces using p and p_np

                for ( i in seq_along(x) ) { for ( j in seq_along(y) ) { z[i,j]= c*(p[1] * (x[i]**2) + p[2] * (y[j]**2))  + (p[3]*x[i] + p[4]*y[j]) } }


                for ( i in seq_along(x) ) { for ( j in seq_along(y) ) { z_np[i,j]= c*(p_np[1] * (x[i]**2) + p_np[2] * (y[j]**2))  + (p_np[3]*x[i] + p_np[4]*y[j]) } }


#                z_mat <- paste0("z_from_lm_",surfaces[k])
                z_mat <- paste0("z_from_lm_",surfaces[k],"_rel_mesh")
                assign(z_mat,z)
#                save(list=z_mat,file=paste0(results_root,"/results/parameter_sweeps/ann_replots/z_from_lm_",surfaces[k],".rda"))
                save(list=z_mat,file=paste0(results_root,"/results/parameter_sweeps/ann_replots/z_from_lm_",surfaces[k],"_rel_mesh.rda"))


#                z_mat_np <- paste0("z_np_from_lm_",surfaces[k])
                z_mat_np <- paste0("z_np_from_lm_",surfaces[k],"_rel_mesh")
                assign(z_mat_np,z_np)
#                save(list=z_mat_np,file=paste0(results_root,"/results/parameter_sweeps/ann_replots/z_np_from_lm_",surfaces[k],".rda"))
                save(list=z_mat_np,file=paste0(results_root,"/results/parameter_sweeps/ann_replots/z_np_from_lm_",surfaces[k],"_rel_mesh.rda"))

                }
        }













# this is the same as the peak-finding in sweep_fcns.r:extrema/3, but I just
# refactored it here
#
# Kazic, 4.2.2017

get_peak <- function(surface) {

        maxl = which(surface == max(surface), arr.ind = TRUE) 
        if ( identical(nrow(maxl),as.integer(1),num.eq=TRUE) ) {
	        surface_max <- max(surface)
	} else {
	        m = floor(nrow(maxl)/2)              # middle of peak if a very flat one
                midmax <- t(as.matrix(maxl[m,]))
                surface_max <- surface[midmax]
                }
		
	return(surface_max)
        }





                
# grab_mesh <- function(x,y,z,contour_set) {
#        sim_mesh_list <- define_mesh(x,y,z,contour_set,15)


grab_mesh <- function(mesh_list) {
	mesh_rows <- nrow(mesh_list$mesh)
        zs_mesh <- rownames(mesh_list$mesh)
	cz_levels <- unlist(regmatches(zs_mesh,gregexpr("c[0-9]{1,2}",zs_mesh,perl=TRUE)))
	z_levels <- unlist(regmatches(cz_levels,gregexpr("[0-9]{1,2}",cz_levels,perl=TRUE)))		
        mesh <- matrix(c(as.numeric(z_levels),mesh_list$mesh[,1],mesh_list$mesh[,2]),mesh_rows,3,byrow=FALSE)

        return(mesh)	
        }





















# hmmm, lotta hard-wiring here, not sure how it will work
#
# Kazic, 26.3.2017


get_mo17_ray0 <- function(x,y,z) {

        contours <- contourLines(x,y,z,levels=get("mo"))
        lsl <- lsfit(x=c(-42,contours[[9]]$x[111],contours[[10]]$x[24]),y=c(7.5,contours[[9]]$y[111],contours[[10]]$y[24]))

        left_bot <- lsfit(x=contours[[9]]$x[1:97],y=contours[[9]]$y[1:97])
        rt_top <- lsfit(x=contours[[10]]$x[38:89],y=contours[[10]]$y[38:89])
        left_top <- lsfit(x=contours[[9]]$x[124:128],y=contours[[9]]$y[124:128])
        rt_bot <- lsfit(x=contours[[10]]$x[1:4],y=contours[[10]]$y[1:4])
        new_x <- (left_bot$coefficients[1] - left_top$coefficients[1])/(left_top$coefficients[2] - left_bot$coefficients[2])
        new_y <- left_top$coefficients[2]*new_x + left_top$coefficients[1]
        xr <- (rt_bot$coefficients[1] - rt_top$coefficients[1])/(rt_top$coefficients[2] - rt_bot$coefficients[2])
        yr <- rt_top$coefficients[2]*xr + rt_top$coefficients[1]
        ray0 <- lsfit(x=c(-42,xl,xr),y=c(7.5,yl,yr))
        return(ray0)
        }












# for the frst paper, plot the 3d surfaces and images with a consistent
# scale, using the new viridis color scheme found at
# https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html
#
# note that x, y, zlimits, ticks, etc. are all loaded into the global environment.
# pal is too, but I return it since that's what the load_standard_plotting_parameters/1
# function does.
#
# Kazic, 21.7.2017

# just plot the images (don't contour lines), but all on the same z scale



# plot_surfaces_n_images("../results/parameter_sweeps/ann_replots/",c("b73","mo17","x1a","x1b","x2a","x2b","x3a","x3b"))

plot_surfaces_n_images <- function(dir,surfaces) {
        pal <- load_standard_plotting_parameters("..")
        pal_len_plus <- length(pal) + 1

        breaks <- seq(-250,150,0.4)
        nlevel <- 7


        vars <- c(outer(c("z","z_from_lm_","z_np_from_lm_"),surfaces,paste0))

        for ( i in vars ) {
	        file <- paste0(dir,i,".rda")
                if ( file.exists(file) ) {
		        load(file)
			if ( exists(i) ) { z <- get(i) }
        

## # plot the surfaces

##                         persp3d(x,y,z,col=pal,xlim=c(-45,45),ylim=c(-8,8),zlim=zlimits,forceClipregion=TRUE,xlab="",ylab="",zlab="",axes=FALSE,box=FALSE,lit=TRUE,specular="black")
		
##                         segments3d(x=as.vector(t(m_for_z[,c(1,4)])),y=as.vector(t(m_for_z[,c(2,5)])),z=as.vector(t(m_for_z[,c(3,6)])))
##                         axis3d('x--',pos=c(-45,-8,0)) ; axis3d('y--',pos=c(-45,-8,0),at=yticks) ; axis3d('z-+',pos=c(-45,8,zlimits[1]),at=zticks)
##                         mtext3d("water",'x--', line = 0.75, at = NULL, pos = NA) ; mtext3d("nitrogen",'y--', line = 0.5, at = NULL, pos = NA) ; mtext3d("z",'z-+', line = 1.25, at = NULL, pos = NA)
## 		        title3d(i)
		
		
##                         writeWebGL(file=paste0(dir,"/",i,".html"),prefix="",font="Arial",width=700,height=700)
## 		        snapshot3d(file=paste0(dir,"/",i,"_snapshot.png"),fmt="png",top=TRUE)
## 		        clear3d()


# plot the images over the same range of z

                        png(paste0(dir,i,"_same_range_image.png"))
                        image.plot(x,y,z,add=FALSE,col=pal,breaks=breaks,nlevel=nlevel,graphics.reset=FALSE,horizontal=TRUE,xlab="water",ylab="nitrogen",main=paste0(i," projected into xy"))
                        dev.off()



# plot the images using the range of each z

                        
			zbreaks <- seq(min(z),max(z),length.out=pal_len_plus)

                        
                        if ( length(zbreaks) > length(pal) ) {
                                png(paste0(dir,i,"_diff_range_image.png"))
                                image.plot(x,y,z,add=FALSE,col=pal,breaks=zbreaks,nlevel=nlevel,graphics.reset=FALSE,horizontal=TRUE,xlab="water",ylab="nitrogen",main=paste0(i," projected into xy"))
                                dev.off()
                                } else { cat(length(zbreaks)," for surface ",i,"\n") }


                        } else { cat("file ",file," does not exist!\n") }
                }
        }








############### stuff to functionalize ###############




# pull together routine analysis of surface fits

# fit_n_test <- function(local,sweep_num,swing) {
#         td <- paste0(local,"artistry/papers/current/ann_stapleton_ms/our_parts/results/parameter_sweeps",sweep_num)
# 	system(paste0("df -kh >> ",file.path(td,"tmp",sweep_num,".r"))
#         system(paste0("./modified_eqn.r ",local,sweep_num,swing))
#         system(paste0("./discover_best_fit_surfaces.r ",local,sweep_num))
# stopped here

# }



# hmm, not quite right.  In the z matrix, we want the peak to be at low
# n and low water; so, the upper left quadrant of z.  But the value for
# col below pushes the search elsewhere.
#
# Kazic, 13.12.2015


# s1 <- which(results["curvature",1:data_res] > -2.5,arr.ind=TRUE)
# s2 <- which(results["col",1:data_res] > 21,arr.ind=TRUE)
# s3 <- which(results["row",1:data_res] < 110,arr.ind=TRUE)
# low_resources <- intersect(s1,s2)
# low_res_curv <- intersect(low_resources,s3)


# # peak position and volume, sizing points by bins of curvature
# # for sizing,
# # http://stackoverflow.com/questions/10341963/3d-scatterplot-in-r-using-rgl-plot3d-different-size-for-each-data-point
# #
# # for multiple conditions in the with statement (stopped here),
# # https://statmethods.wordpress.com/2012/01/30/getting-fancy-with-3-d-scatterplots/

# sizes <- as.numeric(cut(results["curvature",1:data_res], 7))*10
# bins <- split(results["curvature",1:data_res], sizes)
# open3d()
# plot3d(as.numeric(results["col",1:data_res]),as.numeric(results["row",1:data_res]),as.numeric(results["surface_vol",1:data_res]),main="peak postn and volume, colored by curvature",ylim=rev(range(results["row",1:data_res])),xlab="postn nitrogen peak in surface",ylab="postn water peak in surface",zlab="surface volume",pch=19,col="hotpink",size=0)
# with(as.data.frame(results),text3d(results["col",1:data_res],results["row",1:data_res],results["surface_vol",1:data_res],combos,adj = c(1,1)))

# for ( i in seq_along(bins) ) {
#            with(bins[[i]],points3d(results["col",1:data_res],results["row",1:data_res],results["surface_vol",1:data_res],size=i))
#            }







############ debugger calls ##################


# a convenient wrapper to test functions, modify as needed

dummy <- function() {

        date <- Sys.time()
        start <- as.integer(date)

	log_file <- "test.org"
	surfaces_timestamp <- "1111"
        surfaces_file <- "dummy"
	contour_set <- "domed"


	test_set <- c("s55z4691","s55z4692","s55z4693","s55z4694","s55z4695","s55z5161","s55z5162","s55z5163","s55z5331","s55z5501","s55z5502","s55z5681","s55z5682","s55z5683")



        sweep_root <- toggle_partition(0)
        load(paste0(sweep_root,"/results/parameter_sweeps/ann_replots/exptl_meshes.rda"))
        load(paste0(sweep_root,"/results/parameter_sweeps/ann_replots/exptl_meshes_tops.rda"))
        analysis_dir <- paste0(sweep_root,"/results/parameter_sweeps/analysis/")
        b <- nrow(exptl_meshes)


        load(paste0(analysis_dir,"meshes.rda"))
        surfaces <- the_meshes[[1]]
        contours <- the_meshes[[2]]


        summary_matrix <- compare_surfaces_by_z_dists(test_set,surfaces,contours,contour_set,exptl_meshes,exptl_meshes_tops,"1111",analysis_dir)


        stop <- as.integer(Sys.time())

        write_summary_matches(analysis_dir,log_file,date,start,stop,surfaces_timestamp,surfaces_file,surfaces,contour_set,summary_matrix)




#        load("../results/parameter_sweeps/analysis/shape_dists_1477870914.rda")
#        labelled_mat <- dists_matrix_1477870914
#        mat_stem <- paste0(analysis_dir,"shape_dists_",surfaces_timestamp)
#        summary_matrix <- partition_n_stripchart(mat_stem,labelled_mat,surfaces,b)
        }











############ debugging and breakpoint calls


# debug(compare_shapes_by_z_dists)
# debug(compare_surfaces_by_z_dists)
# debug(partition_n_stripchart)
# debug(write_summary_matches)






# for timestamp 1477870914, x2a has no variances; x3a has a few

# dummy("1477870914","domed","shapes.org")
# dummy()


# compare_shapes_by_z_dists("1478728934","domed","shapes.org")
# compare_shapes_by_z_dists("1477870914","domed","shapes.org")



# source("std_view.r")





# setBreakpoint("analysis_fcns.r#6545")



surfaces <- c("b73","mo17","x1a","x1b","x2a","x2b","x3a","x3b")

exptl_surfaces <- c("b73","x1a","x2b","x3b","x1b","x2a","x3a","mo17")
sim_surfaces <- c("b73","x1a","x2b","x1b","x2a","x3a")




# make_n_save_relative_meshes(c("mo17"),"test.org",10)
# make_n_save_relative_meshes(c("b73","x1a","x1b","x2a","x2b","x3a","x3b","mo17"),"test.org",10)
# simulate_surfaces_using_lm(surfaces,"..")
# check_lm_fit_quality(exptl_surfaces,sim_surfaces,"foo.org",10)


# check_experimental_surfs_vs_each_other(exptl_surfaces,"foo.org",10)

plot_heatmaps("../results/parameter_sweeps/ann_replots")

# plot_surfaces_n_images("../results/parameter_sweeps/ann_replots/",c("b73","mo17","x1a","x1b","x2a","x2b","x3a","x3b"))



