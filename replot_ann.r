#!/usr/local/bin/Rscript


rm(list=ls(all=TRUE))
require(rgl)
setwd("~/me/artistry/papers/current/ann_stapleton_ms/first_manuscript/working_draft")
exptl <- read.csv("ann_qtl_restructured.csv")
wd <- "../../our_parts/results/parameter_sweeps/ann_replots/"
swing <- 0.01


source("../../our_parts/simulations/sweep_fcns.r")
x=seq(-42,42,0.5) ;  y=seq(-7.5,7.5,0.5)
proxies_matrix = matrix(NA,nrow=length(proxies),ncol=nrow(exptl))



# don't clear the lights this time

open3d(scale=c(16/90,16/16,16/45))            
bg3d(color="white")
pp <- dget("../../our_parts/simulations/tilted_ann6.r")
par3d(pp)







# max of her surfaces is <40, but it looks like the peak is higher, so extend zlimits

zlimits <- c(-5,45)

m <- matrix(c(-45,-8,0,45,-8,0,
              -45,-8,0,-45,8,0,
              -45,8,zlimits[1],-45,8,zlimits[2]),nrow=3,byrow=TRUE)


# get some extra ticks on z
# from http://stackoverflow.com/questions/10811669/r-rgl-distance-between-axis-ticks-and-tick-labels
#
# my.ticks <- pretty(z, n=5)
# axis3d('y', at=my.ticks, labels=rep("", 5))
# mtext3d(paste(my.ticks), at=my.ticks, edge='y', line=.6)

yticks <- c(-6,-3,0,3,6)
zticks <- c(-5,5,15,25,35,45)



# adjust colors: black light makes the tones darker
#
# old B73 too dark:  
# 90 darkorange #FF8C00 255 140 0
#
# tan too light: 620 tan #D2B48C 210 180 140
#
# orange still too dark 498 orange #FFA500 255 165 0
#
# too dark: 621 tan1 #FFA54F 255 165 79
#
# about the best: 76 darkgoldenrod1 #FFB90F 255 185 15
#
#
# for mo17, not 616 steelblue1 #63B8FF 99 184 255
# too dark: 121 deepskyblue #00BFFF 0 191 255
# 68 cyan #00FFFF 0 255 255

# b73col <- rgb(255,185,15,0.5, max = 255)
# mo17col <- rgb(0,255,255,0.5, max = 255)


                                        # b73col <- "orange1"
b73col <- "darkorange"
mo17col <- "dodgerblue"









for ( r in 1:nrow(exptl) ) {

    
        z=matrix(NA,nrow=length(x),ncol=length(y))

        for ( i in seq_along(x) ) { for ( j in seq_along(y) ) { z[i,j]= exptl[r,"inter"] + exptl[r,"water"]*x[i] + exptl[r,"water2"]*x[i]**2 + exptl[r,"nlevel"]*y[j]  + exptl[r,"nlevel2"]*y[j]**2 + exptl[r,"nbywater"]*x[i]*y[j] } }

        zfile <- paste0("z",as.character(exptl[r,1]),".RData",collapse="")
        save(z,file=file.path(wd,zfile))


        z_proxies <- compute_proxies(z,length(proxies),swing)
        proxies_matrix[,r] <- z_proxies
        rownames(proxies_matrix) <- rownames(z_proxies)              # ugh, but I can't figure out a better way


        
        if ( r %% 2 == 0 ) { coll <- mo17col } else { coll <- b73col }
        print(paste0("r: ",r," col: ",coll))

        persp3d(x,y,z,col=coll,xlim=c(-45,45),ylim=c(-8,8),zlim=zlimits,forceClipregion=TRUE,xlab="",ylab="",zlab="",axes=FALSE,box=FALSE,lit=TRUE,specular="black")

        segments3d(x=as.vector(t(m[,c(1,4)])),y=as.vector(t(m[,c(2,5)])),z=as.vector(t(m[,c(3,6)])))
        axis3d('x--',pos=c(-45,-8,0))
        axis3d('y--',pos=c(-45,-8,0),at=yticks)
        axis3d('z-+',pos=c(-45,8,zlimits[1]),at=zticks)
    
        mtext3d("water",'x--', line = 0.75, at = NULL, pos = NA)
        mtext3d("nitrogen",'y--', line = 0.5, at = NULL, pos = NA)
        mtext3d("z",'z-+', line = 1.25, at = NULL, pos = NA)

    
        writeWebGL(file=file.path(wd,paste0("z",as.character(exptl[r,1]),".html")),prefix="",font="Arial",width=700,height=700)
        snapshot3d(file=file.path(wd,paste0("z",as.character(exptl[r,1]),"_snapshot.png")),fmt="png",top=TRUE)

        
        clear3d()
        rm(z)
        }



summaries <- summary_statistics(proxies_matrix,wd)




