# this is ../artistry/papers/current/ann_stapleton_ms/our_parts/simulations/estimate_exptl_parameters.r



# Given the experimental surfaces for the alleles, estimate ac, bc, d, and
# e using a linear model.  This simplification lets us use a constant c and
# permits the use of a linear, rather than a nonlinear system.  Try both
# qr.solve and limSolve, which should give the same results within
# numerical errors.  Then factor out an arbitrary c so that the a and b
# values are reasonable, given our previous simulation results.
#
#
# The modified equation is z = c(ax^2 + by^2) + dx + ey.  Writing this as an
# inhomogeneous system and combining parameters, we have
#
#    caX^2 + cbY^2 + dX + eY = Z, 
#
# where X, Y, and Z are matrices.  For the domed surfaces, we have values for
# these variables at each of the mesh points calculated using
# analysis_fcns.r:define_mesh/5, with thetadeg = 15.  These are stored in
# [[file:../results/parameter_sweeps/ann_replots/exptl_meshes.rda][exptl_meshes.rda]].
# (I have to calculate the mesh points for Mo17 yet. DONE!!!)
#
# So given the vector of mesh points, X, Y, and Z become matrices of
# constants, and the problem asks for values of ca, cb, d, and e that make
# all the linear equations simultaneously true.  In effect, we swap the
# roles of the coefficients and variables.
#
# We have 4 unknowns and 20 -- 30+ equations, including the equation
# defining the value of the peak of each surface.  So it seems likely that
# each experimental surface's system will be overdetermined.  In principle
# one might simply exclude candidate equations, for example picking one
# equation from each contour.  But both QR factorization and limSolve
# should take care of this problem automatically.  So let's try that,
# asking first for an exact solution.
#
# Notice that we are (finally!) optimizing our fits of the equation to the
# surfaces.  So we might find the model doesn't work very well after all.
# We partially forestall this possibility by ignoring regions of the
# surfaces (the top, left, and bottom edges) where we know the model
# doesn't perform well.  But we could be in for some surprises.
#
#
# [[https://en.m.wikibooks.org/wiki/Linear_Algebra][linear algebra wikibook]]
# [[https://en.m.wikipedia.org/wiki/System_of_linear_equations][linear systems]]
# [[https://en.m.wikipedia.org/wiki/QR_decomposition][outline QR decomposition]]
# [[https://stat/ethz.ch/R-manual/R-devel/library/base/html/qr.html][qr]]
# [[https://cran.r-project.org/web/packages/limSolve/index.html][limSolve]]
#
# Kazic, 10.1.2017



rm(list=ls(all=TRUE))

require(limSolve)
require(graphics)




this_dir <- getwd()
source(paste0(this_dir,"/analysis_fcns.r"))
source(paste0(this_dir,"/sweep_fcns.r"))
local_root <- toggle_partition(0)
ann_dir <- paste0(local_root,"/results/parameter_sweeps/ann_replots")

this_dir <- getwd()
if ( this_dir != paste0(local_root,"/simulations") ) { setwd(paste0(local_root,"our_parts/simulations/")) }




surfaces <- c("mo17","b73","x1a","x2b","x3b","x2a","x1b","x3a")






################## setting up equation systems


# grab the peaks for each surface and save to a data structure, rather like
# that for the mesh points
#
# grabbing computation from sweep_fcns.r:extrema/3 and calls to it, a bit
# excessive but saves refactoring.
#
# Need do this only once.
#
# Kazic, 10.1.2017
#
#
# Well, need to get the values for x and y, not just their positions in the
# matrix, so ended up doing twice.  And then rearranged the matrices to fit
# the changed order of the surfaces, since that surface order groups the
# experimental surfaces nicely.
#
# Kazic, 12.1.2017



grab_maxima <- function(surfaces) {

        maxima <- matrix(NA,nrow=3,ncol=length(surfaces))
	rownames(maxima) <- c("x","y","z")
	colnames(maxima) <- surfaces
        x=seq(-42,42,0.5) ;  y=seq(-7.5,7.5,0.5)


        for ( i in 1:length(surfaces) ) {
                load(paste0(ann_dir,"/z",surfaces[i],".rda"))


# yes, silly, but safer than assuming all matrices will be the same

                rows <- nrow(z)
	        cols <- ncol(z)
                maxl = which(z == max(z), arr.ind = TRUE)
                corners <- cbind(z[1,1], z[1,cols], z[rows,1], z[rows,cols])

                results <- extrema(z,maxl,corners)

# oops, want the actual values for x and y, not their position in the
# matrix
#
# Kazic, 12.1.2017
#
#               maxima[1,i] <- results["row",1]
#		maxima[2,i] <- results["col",1]

                maxima[1,i] <- x[results["row",1]]
		maxima[2,i] <- y[results["col",1]]

		maxima[3,i] <- results["z_max",1]
		rm(z)
                }
		
        exptl_surfaces_maxima <- maxima
	save(exptl_surfaces_maxima,file=paste0(ann_dir,"/exptl_surfaces_maxima.rda"))
        }





# now covers mo17 case
#
# Kazic, 18.1.2017


form_linear_systems <- function(surfaces) {

        for ( i in 1:length(surfaces) ) {

                if ( surfaces[i] == "mo17" ) {
          	        mesh_data <- mo17_mesh
                        } else { mesh_data <- exptl_meshes }


                total_rows <- nrow(mesh_data)
	        colx <- paste0(surfaces[i],":x")
	        coly <- paste0(surfaces[i],":y")		
                top <- exptl_meshes_tops[colx]
		rows_needed <- total_rows - top + 2
                lin_mat <- matrix(NA,nrow=rows_needed,ncol=5)


                lin_mat[1,] <- c(exptl_surfaces_maxima[1,surfaces[i]]^2,exptl_surfaces_maxima[2,surfaces[i]]^2,exptl_surfaces_maxima[1,surfaces[i]],exptl_surfaces_maxima[2,surfaces[i]],exptl_surfaces_maxima[3,surfaces[i]])
                lin_mat[2:rows_needed,] <- c(mesh_data[top:total_rows,colx]^2,mesh_data[top:total_rows,coly]^2,mesh_data[top:total_rows,colx],mesh_data[top:total_rows,coly],mesh_data[top:total_rows,1])

                colnames(lin_mat) <- c("x**2","y**2","x","y","z")
		rownames(lin_mat) <- c("peak",rownames(mesh_data)[top:total_rows])
		mat <- paste0(surfaces[i],"_lm")
		assign(mat,lin_mat)
		save(list=mat,file=paste0(ann_dir,"/",mat,".rda"))
                }
        }








# form linear models using the relative mesh points calculated
# using analysis_fcns.r:make_n_save_relative_meshes/4
#
# Kazic, 19.7.2017


form_linear_systems_relmesh <- function(surfaces) {

        ann_dir <- "../results/parameter_sweeps/ann_replots/"
        load(paste0(ann_dir,"exptl_surfaces_maxima.rda"))


        for ( i in 1:length(surfaces) ) {

                mesh_name <- paste0("s_",surfaces[i],"_10pts_exp_mesh")
                load(paste0(ann_dir,mesh_name,".rda"))
                unstacked_mesh <- get(mesh_name)
                mesh <- stack_matrix(unstacked_mesh)
                mesh_rows <- nrow(mesh)
		total_rows <- mesh_rows + 1
                lin_mat <- matrix(NA,nrow=total_rows,ncol=5)


                lin_mat[1,] <- c(exptl_surfaces_maxima[1,surfaces[i]]^2,exptl_surfaces_maxima[2,surfaces[i]]^2,exptl_surfaces_maxima[1,surfaces[i]],exptl_surfaces_maxima[2,surfaces[i]],exptl_surfaces_maxima[3,surfaces[i]])
                lin_mat[2:total_rows,] <- c(mesh[,1]^2,mesh[,2]^2,mesh[,1],mesh[,2],mesh[,3])



# nasty; never named the rows of the relative mesh, do it here
#
# Kazic, 19.7.2017

                mesh_row_names <- c(t(outer(paste0("r",seq(0,5,1),":"),seq(1,10,1),paste0)))

                colnames(lin_mat) <- c("x**2","y**2","x","y","z")
		rownames(lin_mat) <- c("peak",mesh_row_names)
		mat <- paste0(surfaces[i],"_rel_mesh_lm")
		assign(mat,lin_mat)
		save(list=mat,file=paste0(ann_dir,"/",mat,".rda"))
                }
        }


















# notice use of get() to fetch value of matrix and assign it to a variable


fit_linear_models <- function(surfaces,outfile) {
        file <- paste0(ann_dir,"/",outfile)
        file_handle <- file(file,"a")
        date <- Sys.time()
        timestamp <- as.integer(date)


        parms <- matrix(NA,nrow=length(surfaces),ncol=5)
        rownames(parms) <- surfaces
	colnames(parms) <- c("a","b","c","d","e")


        cat("\n\n\n* Results of Linear Model Fitting for Timestamp",paste0(timestamp),"\n\nRun on",as.character(date),"\n\n",file=file_handle)
        cat("\n\n** fit data\n\n",file=file_handle)
	cat("#+tblname: fit_data_",as.character(timestamp),"\n",file=file_handle)
        cat("| surface | method | ca | cb | d  | e | residualNorm | solutionNorm | IsError | type | sqrt(solutionNorm) |\n",file=file_handle)
        cat("|-----+----+------+---+----+---+---+--+---+-----+-----|\n",file=file_handle)		        


        for ( i in 1:length(surfaces) ) {
        
                load(paste0(ann_dir,"/",surfaces[i],"_lm.rda"))
		lm_mat <- get(paste0(surfaces[i],"_lm"))
                f <- Solve(lm_mat[,1:4],lm_mat[,5])
                g <- lsei(A=lm_mat[,1:4],B=lm_mat[,5],fulloutput=TRUE,verbose=FALSE)
                root <- sqrt(g$solutionNorm)
		
                cat("| ",surfaces[i]," | Solve | ",f[1]," | ",f[2]," | ",f[3]," | ",f[4]," | \n",file=file_handle)		        
		cat("| | lsei | ",g$X[1]," | ",g$X[2]," | ",g$X[3]," | ",g$X[4]," | ",g$residualNorm," | ",g$solutionNorm," | ",g$IsError," | ",g$type," | ",root," |\n",file=file_handle)


                if ( surfaces[i] == "mo17" ) { c <- 1 } else { c <- -1 }
                a <- c * f[1]
		b <- c * f[2]
                parms[i,] <- c(a,b,c,f[3],f[4])
                }

        cat("\n\n** summary table\n\n",file=file_handle)
	cat("Manually inserted the column of leading |s below the table for nice formatting.\n\n\n",file=file_handle) 
	cat("#+tblname: parm_summary_",as.character(timestamp),"\n",file=file_handle) 
        cat("| surface | a | b | c | d | e |\n",file=file_handle) 
        cat("|-----+----+------+---+----+--|\n",file=file_handle)		        
        write.table(parms,file=file_handle,sep=" | ",row.names=TRUE,col.names=FALSE,eol=" |\n",quote=FALSE)
        cat("\n\n|\n|\n|\n|\n|\n|\n|\n|\n",file=file_handle)		        
        close(file_handle)
        }









# revised to use the linear models calculated from the relative meshes;
# only revisions are to data file and variable names.
#
# Kazic, 19.7.2017



fit_linear_models_relmesh <- function(surfaces,outfile) {
        file <- paste0(ann_dir,"/",outfile)
        file_handle <- file(file,"a")
        date <- Sys.time()
        timestamp <- as.integer(date)


        parms <- matrix(NA,nrow=length(surfaces),ncol=5)
        rownames(parms) <- surfaces
	colnames(parms) <- c("a","b","c","d","e")


        cat("\n\n\n* Results of Linear Model Fitting using Relative Meshes for Timestamp",paste0(timestamp),"\n\nRun on",as.character(date),"\n\n",file=file_handle)
        cat("\n\n** fit data\n\n",file=file_handle)
	cat("#+tblname: fit_data_",as.character(timestamp),"\n",file=file_handle)
        cat("| surface | method | ca | cb | d  | e | residualNorm | solutionNorm | IsError | type | sqrt(solutionNorm) |\n",file=file_handle)
        cat("|-----+----+------+---+----+---+---+--+---+-----+-----|\n",file=file_handle)		        


        for ( i in 1:length(surfaces) ) {
        
                load(paste0(ann_dir,"/",surfaces[i],"_rel_mesh_lm.rda"))
		lm_mat <- get(paste0(surfaces[i],"_rel_mesh_lm"))
                f <- Solve(lm_mat[,1:4],lm_mat[,5])
                g <- lsei(A=lm_mat[,1:4],B=lm_mat[,5],fulloutput=TRUE,verbose=FALSE)
                root <- sqrt(g$solutionNorm)
		
                cat("| ",surfaces[i]," | Solve | ",f[1]," | ",f[2]," | ",f[3]," | ",f[4]," | \n",file=file_handle)		        
		cat("| | lsei | ",g$X[1]," | ",g$X[2]," | ",g$X[3]," | ",g$X[4]," | ",g$residualNorm," | ",g$solutionNorm," | ",g$IsError," | ",g$type," | ",root," |\n",file=file_handle)


                if ( surfaces[i] == "mo17" ) { c <- 1 } else { c <- -1 }
                a <- c * f[1]
		b <- c * f[2]
                parms[i,] <- c(a,b,c,f[3],f[4])
                }

        cat("\n\n** summary table\n\n",file=file_handle)
	cat("Manually inserted the column of leading |s below the table for nice formatting.\n\n\n",file=file_handle) 
	cat("#+tblname: parm_summary_",as.character(timestamp),"\n",file=file_handle) 
        cat("| surface | a | b | c | d | e |\n",file=file_handle) 
        cat("|-----+----+------+---+----+--|\n",file=file_handle)		        
        write.table(parms,file=file_handle,sep=" | ",row.names=TRUE,col.names=FALSE,eol=" |\n",quote=FALSE)
        cat("\n\n|\n|\n|\n|\n|\n|\n|\n|\n",file=file_handle)		        
        close(file_handle)
        }






































############ bunch of directives here


# grab_maxima(surfaces)
# load(paste0(ann_dir,"/exptl_surfaces_maxima.rda"))
# load(paste0(ann_dir,"/exptl_meshes_tops.rda"))
# load(paste0(ann_dir,"/exptl_meshes.rda"))
# load(paste0(ann_dir,"/mo17_mesh.rda"))



# setBreakpoint("estimate_exptl_parameters.r#206")

# form_linear_systems_relmesh(surfaces)
# fit_linear_models_relmesh(surfaces,"fits.org")
