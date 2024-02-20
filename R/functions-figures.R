###############
# DEPENDENCIES
###############
library(extrafont)
library(fontcm)
loadfonts(quiet = TRUE)
library(wesanderson)
#library(latex2exp)

#######################
# AUXILLIARY FUNCTIONS
#######################

toPdf <- function(expr, filename, ...) {
  toDev(expr, pdf, filename, ...)
}

figPath  <-  function(name) {
  file.path('./figures', name)
}

toDev <- function(expr, dev, filename, ..., verbose=TRUE) {
  if ( verbose )
    cat(sprintf('Creating %s\n', filename))
    dev(filename, family='CM Roman', ...)
#    dev(filename, family='Arial', ...)
    on.exit(dev.off())
    eval.parent(substitute(expr))
}



####################
# PLOTTING FUNCTIONS
####################

#' Plot text or points according to relative axis position.
#'
#' @title Plot text or points according to relative axis position
#' @param px Relative x-axis position (in proportion) where character is to be plotted.
#' @param py Relative y-axis position (in proportion) where character is to be plotted.
#' @param lab Plotted text. Works if argument \code{\link[graphics]{text}} is \code{TRUE}.
#' @param adj See argument of same name in R base function \code{\link[graphics]{par}}.
#' @param text Logical. Should text or points be plotted?
#' @param log Used if the original plot uses the argument log, e.g. \code{log='x'}, \code{log='y'} or \code{log='xy'}.
#' @param ... Additional arguments to R base function \code{\link[graphics]{text}}.
#' @export
proportionalLabel <- function(px, py, lab, adj=c(0, 1), text=TRUE, log=FALSE, ...) {
    usr  <-  par('usr')
    x.p  <-  usr[1] + px*(usr[2] - usr[1])
    y.p  <-  usr[3] + py*(usr[4] - usr[3])
    if(log=='x') {
        x.p<-10^(x.p)
    }
    if(log=='y') {
        y.p<-10^(y.p)
    }
    if(log=='xy') {
        x.p<-10^(x.p)
        y.p<-10^(y.p)
    }
    if(text){
        text(x.p, y.p, lab, adj=adj, ...)
    } else {
        points(x.p, y.p, ...)
    }
}



proportionalArrows <- function(px1, py1, px2, py2, adj=c(0, 1), log=FALSE, length=length, ...) {
    usr  <-  par('usr')
    x.p1  <-  usr[1] + px1*(usr[2] - usr[1])
    y.p1  <-  usr[3] + py1*(usr[4] - usr[3])
    x.p2  <-  usr[1] + px2*(usr[2] - usr[1])
    y.p2  <-  usr[3] + py2*(usr[4] - usr[3])
    if(log=='x') {
        x.p1  <-  10^(x.p1)
        x.p2  <-  10^(x.p2)
    }
    if(log=='y') {
        y.p1  <-  10^(y.p1)
        y.p2  <-  10^(y.p2)
    }
    if(log=='xy') {
        x.p1  <-  10^(x.p1)
        y.p1  <-  10^(y.p1)
        x.p2  <-  10^(x.p2)
        y.p2  <-  10^(y.p2)
    }
    arrows(x0=x.p1, y0=y.p1, x1=x.p2, y1=y.p2, length=length,...)
}

#' Draw equally-spaced white lines on plot window.
#'
#' @title Equally-spaced white lines on plot window
#' @param ... Additional arguments to internal function \code{\link{proportionalLabel}}.
#' @author Diego Barneche
#' @export
plotGrid  <-  function(lineCol='white',...) {
    proportionalLabel(rep(0.2, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(rep(0.4, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(rep(0.6, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(rep(0.8, 2), c(0,1), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.2, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.4, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.6, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
    proportionalLabel(c(0,1), rep(0.8, 2), text=FALSE, type='l', col=lineCol, lwd=0.5, ...)
}


#' Internal. Create nice rounded numbers for plotting.
#'
#' @title Rounded numbers for plotting
#' @param value A numeric vector.
#' @param precision Number of rounding digits.
#' @return A character vector.
#' @author Diego Barneche.
rounded  <-  function(value, precision=1) {
  sprintf(paste0('%.', precision, 'f'), round(value, precision))
}


#' Creates transparent colours
#'
#' @title Creates transparent colours
#' @param col Colour.
#' @param opacity equivalent to alpha transparency parameter
#' @export
transparentColor <- function(col, opacity=0.5) {
    if (length(opacity) > 1 && any(is.na(opacity))) {
        n        <-  max(length(col), length(opacity))
        opacity  <-  rep(opacity, length.out=n)
        col      <-  rep(col, length.out=n)
        ok       <-  !is.na(opacity)
        ret      <-  rep(NA, length(col))
        ret[ok]  <-  Recall(col[ok], opacity[ok])
        ret
    } else {
        tmp  <-  col2rgb(col)/255
        rgb(tmp[1,], tmp[2,], tmp[3,], alpha=opacity)
    }
}

fibonacci.scale  <-  function(n) {
    fibs  <-  c(0,1)
    for(i in 2:n) {
        fibs  <-  c(fibs, (fibs[i] + fibs[i-1]))
    }
    (fibs/max(fibs))[-2]
}



##############################################################
##############################################################
##  2-locus model figures

#' Heterozygote Advantage Figure -- 2 locus model
#' Reproducing Fig. 1a from Charlesworth & Wall (1999)
heterozygoteAdvantageFig  <-  function(df = "./data/hetAdvSimData_r0.5.csv") {

    ## import data.frames
    dat  <-  read.csv(file = df)
    d1   <-  dat[dat$s == 0.1,]
    d1   <-  d1[d1$t == 0.1,]
    d2   <-  dat[dat$s == 0.2,]
    d2   <-  d2[d2$t == 0.2,]
    d3   <-  dat[dat$s == 0.1,]
    d3   <-  d3[d3$t == 0.2,]
    d4   <-  dat[dat$s == 0.2,]
    d4   <-  d4[d4$t == 0.1,]

    # set plot layout
    layout.mat <- matrix(c(1), nrow=1, ncol=1, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

    # make plot
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(4,5,5,2), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(-0.03,0.08), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(invSelAdv ~ beta, lwd=1.5, lty=1, col=transparentColor('#252525', opacity=1), data=d1)
        points(invSelAdv ~ beta, pch=21, bg=transparentColor('#252525', opacity=0.75), data=d1)
        lines(invSelAdv ~ beta, lwd=1.5, lty=2, col=transparentColor('#252525', opacity=1), data=d2)
        points(invSelAdv ~ beta, pch=22, bg=transparentColor('#252525', opacity=0.75), data=d2)
        lines(invSelAdv[1:9] ~ beta[1:9], lwd=1.5, lty=3, col=transparentColor('#252525', opacity=1), data=d3)
        points(invSelAdv[1:9] ~ beta[1:9], pch=23, bg=transparentColor('#252525', opacity=0.75), data=d3)
        lines(invSelAdv[1:9] ~ beta[1:9], lwd=1.5, lty=4, col=transparentColor('#252525', opacity=1), data=d4)
        points(invSelAdv[1:9] ~ beta[1:9], pch=24, bg=transparentColor('#252525', opacity=0.75), data=d4)
        abline(h=0, lwd=2)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel( 0.5,  1.1,   expression(paste("Heterozygote Advantage")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.2,  0.5,   expression(paste(italic(s[I]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.2,  expression(paste("Frequency of sib matings (", italic(alpha),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2]*0.37,
               y       =  usr[4]*0.95,
               legend  =  c(
                            expression(paste(italic(s), " = ", italic(t), " = ", 0.1)),
                            expression(paste(italic(s), " = ", italic(t), " = ", 0.2)),
                            expression(paste(italic(s), " = ", 0.1, "; ", italic(t), " = ", 0.2)),
                            expression(paste(italic(s), " = ", 0.2, "; ", italic(t), " = ", 0.1))),
               pch     =  c(21,22,23,24),
               pt.bg      =  c(transparentColor('#252525', opacity=0.75),
                               transparentColor('#252525', opacity=0.75),
                               transparentColor('#252525', opacity=0.75),
                               transparentColor('#252525', opacity=0.75)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
}







shelteringFig  <-  function() {

    ## import data.frames
    dat    <-  read.csv(file = "./data/shelteringData.csv")
    r1  <-  dat[dat$r == 0.001,]
    d1  <-  r1[r1$h == 0.0,]
    d1  <-  d1[d1$s == 0.01,]
    d2  <-  r1[r1$h == 0.0,]
    d2  <-  d2[d2$s == 0.05,]
    d3  <-  r1[r1$h == 0.1,]
    d3  <-  d3[d3$s == 0.01,]
    d4  <-  r1[r1$h == 0.1,]
    d4  <-  d4[d4$s == 0.05,]

    # set plot layout
    layout.mat <- matrix(c(1:4), nrow=2, ncol=2, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

    # Panel A
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(-0.006,0.001), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        abline(h=0, lwd=2)
        lines(invSelAdv ~ beta, lwd=1.5, lty=1, col=transparentColor('#252525', opacity=1), data=d1)
        points(invSelAdv ~ beta, pch=21, bg=transparentColor('grey60', opacity=1), data=d1)
        lines(invSelAdv ~ beta, lwd=1.5, lty=2, col=transparentColor('#252525', opacity=1), data=d2)
        points(invSelAdv ~ beta, pch=22, bg=transparentColor('grey60', opacity=1), data=d2)
        lines(invSelAdv ~ beta, lwd=1.5, lty=3, col=transparentColor('#252525', opacity=1), data=d3)
        points(invSelAdv ~ beta, pch=23, bg=transparentColor('grey60', opacity=1), data=d3)
        lines(invSelAdv ~ beta, lwd=1.5, lty=4, col=transparentColor('#252525', opacity=1), data=d4)
        points(invSelAdv ~ beta, pch=24, bg=transparentColor('grey60', opacity=1), data=d4)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, at=c(0.001, 0, -0.002, -0.004, -0.006))
        # Plot labels etc.
        proportionalLabel(0.04, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(r), " = ", 0.001)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(paste(italic(s[I]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
#        proportionalLabel( 0.5,  -0.4,  expression(paste("Frequency Full-Sib. Mating (", italic(beta),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  -usr[4]*1.25,
               legend  =  c(
                            expression(paste(italic(h), " = ", 0, "; ", italic(s), " = ", 0.01)),
                            expression(paste(italic(h), " = ", 0, "; ", italic(s), " = ", 0.05)),
                            expression(paste(italic(h), " = ", 0.1, "; ", italic(s), " = ", 0.01)),
                            expression(paste(italic(h), " = ", 0.1, "; ", italic(s), " = ", 0.05))),
               pch     =  c(21,22,23,24),
               pt.bg      =  c(transparentColor('grey60', opacity=1),
                               transparentColor('grey60', opacity=1),
                               transparentColor('grey60', opacity=1),
                               transparentColor('grey60', opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


    # Panel B: r = 0.01
    r2  <-  dat[dat$r == 0.01,]
    d1  <-  r2[r2$h == 0.0,]
    d1  <-  d1[d1$s == 0.01,]
    d2  <-  r2[r2$h == 0.0,]
    d2  <-  d2[d2$s == 0.05,]
    d3  <-  r2[r2$h == 0.1,]
    d3  <-  d3[d3$s == 0.01,]
    d4  <-  r2[r2$h == 0.1,]
    d4  <-  d4[d4$s == 0.05,]
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(-0.006,0.001), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        abline(h=0, lwd=2)
        lines(invSelAdv ~ beta, lwd=1.5, lty=1, col=transparentColor('#252525', opacity=1), data=d1)
        points(invSelAdv ~ beta, pch=21, bg=transparentColor('grey60', opacity=1), data=d1)
        lines(invSelAdv ~ beta, lwd=1.5, lty=2, col=transparentColor('#252525', opacity=1), data=d2)
        points(invSelAdv ~ beta, pch=22, bg=transparentColor('grey60', opacity=1), data=d2)
        lines(invSelAdv ~ beta, lwd=1.5, lty=3, col=transparentColor('#252525', opacity=1), data=d3)
        points(invSelAdv ~ beta, pch=23, bg=transparentColor('grey60', opacity=1), data=d3)
        lines(invSelAdv ~ beta, lwd=1.5, lty=4, col=transparentColor('#252525', opacity=1), data=d4)
        points(invSelAdv ~ beta, pch=24, bg=transparentColor('grey60', opacity=1), data=d4)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, at=c(0.001, 0, -0.002, -0.004, -0.006), labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.075, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(r), " = ", 0.01)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(-0.4,  0.5,   expression(paste(italic(s[I]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
#        proportionalLabel( 0.5,  -0.4,  expression(paste("Frequency Full-Sib. Mating (", italic(beta),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        

    # Panel C
    r3  <-  dat[dat$r == 0.1,]
    d1  <-  r3[r3$h == 0.0,]
    d1  <-  d1[d1$s == 0.01,]
    d2  <-  r3[r3$h == 0.0,]
    d2  <-  d2[d2$s == 0.05,]
    d3  <-  r3[r3$h == 0.1,]
    d3  <-  d3[d3$s == 0.01,]
    d4  <-  r3[r3$h == 0.1,]
    d4  <-  d4[d4$s == 0.05,]
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(-0.005,0.001), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        abline(h=0, lwd=2)
        lines(invSelAdv ~ beta, lwd=1.5, lty=1, col=transparentColor('#252525', opacity=1), data=d1)
        points(invSelAdv ~ beta, pch=21, bg=transparentColor('grey60', opacity=1), data=d1)
        lines(invSelAdv ~ beta, lwd=1.5, lty=2, col=transparentColor('#252525', opacity=1), data=d2)
        points(invSelAdv ~ beta, pch=22, bg=transparentColor('grey60', opacity=1), data=d2)
        lines(invSelAdv ~ beta, lwd=1.5, lty=3, col=transparentColor('#252525', opacity=1), data=d3)
        points(invSelAdv ~ beta, pch=23, bg=transparentColor('grey60', opacity=1), data=d3)
        lines(invSelAdv ~ beta, lwd=1.5, lty=4, col=transparentColor('#252525', opacity=1), data=d4)
        points(invSelAdv ~ beta, pch=24, bg=transparentColor('grey60', opacity=1), data=d4)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.075, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(r), " = ", 0.1)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(paste(italic(s[I]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5, -0.3,  expression(paste("Frequency of sib matings (", italic(beta),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        


    # Panel D: r = 0.5
    r4  <-  dat[dat$r == 0.5,]
    d1  <-  r4[r4$h == 0.0,]
    d1  <-  d1[d1$s == 0.01,]
    d2  <-  r4[r4$h == 0.0,]
    d2  <-  d2[d2$s == 0.05,]
    d3  <-  r4[r4$h == 0.1,]
    d3  <-  d3[d3$s == 0.01,]
    d4  <-  r4[r4$h == 0.1,]
    d4  <-  d4[d4$s == 0.05,]
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(-0.005,0.001), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        abline(h=0, lwd=2)
        lines(invSelAdv ~ beta, lwd=1.5, lty=1, col=transparentColor('#252525', opacity=1), data=d1)
        points(invSelAdv ~ beta, pch=21, bg=transparentColor('grey60', opacity=1), data=d1)
        lines(invSelAdv ~ beta, lwd=1.5, lty=2, col=transparentColor('#252525', opacity=1), data=d2)
        points(invSelAdv ~ beta, pch=22, bg=transparentColor('grey60', opacity=1), data=d2)
        lines(invSelAdv ~ beta, lwd=1.5, lty=3, col=transparentColor('#252525', opacity=1), data=d3)
        points(invSelAdv ~ beta, pch=23, bg=transparentColor('grey60', opacity=1), data=d3)
        lines(invSelAdv ~ beta, lwd=1.5, lty=4, col=transparentColor('#252525', opacity=1), data=d4)
        points(invSelAdv ~ beta, pch=24, bg=transparentColor('grey60', opacity=1), data=d4)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.075, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(r), " = ", 0.5)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.3,  expression(paste("Frequency of sib matings (", italic(beta),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
}






shelteringSexSpecificFig  <-  function() {

    ## import data.frames
    dat    <-  read.csv(file = "./data/shelteringDataSexSpecific.csv")
    r1  <-  dat[dat$r == 0.001,]
    d1  <-  r1[r1$h == 0.0,]
    d1  <-  d1[d1$sf == 0.01,]
    d2  <-  r1[r1$h == 0.0,]
    d2  <-  d2[d2$sf == 0.05,]
    d3  <-  r1[r1$h == 0.1,]
    d3  <-  d3[d3$sf == 0.01,]
    d4  <-  r1[r1$h == 0.1,]
    d4  <-  d4[d4$sf == 0.05,]

    # set plot layout
    layout.mat <- matrix(c(1:4), nrow=2, ncol=2, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

    # Panel A
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(-0.006,0.001), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        abline(h=0, lwd=2)
        lines(invSelAdv ~ beta, lwd=1.5, lty=1, col=transparentColor('#252525', opacity=1), data=d1)
        points(invSelAdv ~ beta, pch=21, bg=transparentColor('grey60', opacity=1), data=d1)
        lines(invSelAdv ~ beta, lwd=1.5, lty=2, col=transparentColor('#252525', opacity=1), data=d2)
        points(invSelAdv ~ beta, pch=22, bg=transparentColor('grey60', opacity=1), data=d2)
        lines(invSelAdv ~ beta, lwd=1.5, lty=3, col=transparentColor('#252525', opacity=1), data=d3)
        points(invSelAdv ~ beta, pch=23, bg=transparentColor('grey60', opacity=1), data=d3)
        lines(invSelAdv ~ beta, lwd=1.5, lty=4, col=transparentColor('#252525', opacity=1), data=d4)
        points(invSelAdv ~ beta, pch=24, bg=transparentColor('grey60', opacity=1), data=d4)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, at=c(0.001, 0, -0.002, -0.004, -0.006))
        # Plot labels etc.
        proportionalLabel(0.04, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(r), " = ", 0.001)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(paste(italic(s[I]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
#        proportionalLabel( 0.5,  -0.4,  expression(paste("Frequency Full-Sib. Mating (", italic(beta),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  -usr[4]*1.25,
               legend  =  c(
                            expression(paste(italic(h), " = ", 0, "; ", italic(s[i]), " = ", 0.01, ", ", italic(s[m]), " = ", 0.05)),
                            expression(paste(italic(h), " = ", 0, "; ", italic(s[i]), " = ", 0.05, ", ", italic(s[m]), " = ", 0.01)),
                            expression(paste(italic(h), " = ", 0.1, "; ", italic(s[i]), " = ", 0.01, ", ", italic(s[m]), " = ", 0.05)),
                            expression(paste(italic(h), " = ", 0.1, "; ", italic(s[i]), " = ", 0.05, ", ", italic(s[m]), " = ", 0.01))),
               pch     =  c(21,22,23,24),
               pt.bg      =  c(transparentColor('grey60', opacity=1),
                               transparentColor('grey60', opacity=1),
                               transparentColor('grey60', opacity=1),
                               transparentColor('grey60', opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


    # Panel B: r = 0.01
    r2  <-  dat[dat$r == 0.01,]
    d1  <-  r2[r2$h == 0.0,]
    d1  <-  d1[d1$sf == 0.01,]
    d2  <-  r2[r2$h == 0.0,]
    d2  <-  d2[d2$sf == 0.05,]
    d3  <-  r2[r2$h == 0.1,]
    d3  <-  d3[d3$sf == 0.01,]
    d4  <-  r2[r2$h == 0.1,]
    d4  <-  d4[d4$sf == 0.05,]
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(-0.006,0.001), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        abline(h=0, lwd=2)
        lines(invSelAdv ~ beta, lwd=1.5, lty=1, col=transparentColor('#252525', opacity=1), data=d1)
        points(invSelAdv ~ beta, pch=21, bg=transparentColor('grey60', opacity=1), data=d1)
        lines(invSelAdv ~ beta, lwd=1.5, lty=2, col=transparentColor('#252525', opacity=1), data=d2)
        points(invSelAdv ~ beta, pch=22, bg=transparentColor('grey60', opacity=1), data=d2)
        lines(invSelAdv ~ beta, lwd=1.5, lty=3, col=transparentColor('#252525', opacity=1), data=d3)
        points(invSelAdv ~ beta, pch=23, bg=transparentColor('grey60', opacity=1), data=d3)
        lines(invSelAdv ~ beta, lwd=1.5, lty=4, col=transparentColor('#252525', opacity=1), data=d4)
        points(invSelAdv ~ beta, pch=24, bg=transparentColor('grey60', opacity=1), data=d4)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, at=c(0.001, 0, -0.002, -0.004, -0.006), labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.075, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(r), " = ", 0.01)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(-0.4,  0.5,   expression(paste(italic(s[I]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
#        proportionalLabel( 0.5,  -0.4,  expression(paste("Frequency Full-Sib. Mating (", italic(beta),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        

    # Panel C
    r3  <-  dat[dat$r == 0.1,]
    d1  <-  r3[r3$h == 0.0,]
    d1  <-  d1[d1$sf == 0.01,]
    d2  <-  r3[r3$h == 0.0,]
    d2  <-  d2[d2$sf == 0.05,]
    d3  <-  r3[r3$h == 0.1,]
    d3  <-  d3[d3$sf == 0.01,]
    d4  <-  r3[r3$h == 0.1,]
    d4  <-  d4[d4$sf == 0.05,]
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(-0.005,0.001), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        abline(h=0, lwd=2)
        lines(invSelAdv ~ beta, lwd=1.5, lty=1, col=transparentColor('#252525', opacity=1), data=d1)
        points(invSelAdv ~ beta, pch=21, bg=transparentColor('grey60', opacity=1), data=d1)
        lines(invSelAdv ~ beta, lwd=1.5, lty=2, col=transparentColor('#252525', opacity=1), data=d2)
        points(invSelAdv ~ beta, pch=22, bg=transparentColor('grey60', opacity=1), data=d2)
        lines(invSelAdv ~ beta, lwd=1.5, lty=3, col=transparentColor('#252525', opacity=1), data=d3)
        points(invSelAdv ~ beta, pch=23, bg=transparentColor('grey60', opacity=1), data=d3)
        lines(invSelAdv ~ beta, lwd=1.5, lty=4, col=transparentColor('#252525', opacity=1), data=d4)
        points(invSelAdv ~ beta, pch=24, bg=transparentColor('grey60', opacity=1), data=d4)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.075, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(r), " = ", 0.1)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(paste(italic(s[I]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5, -0.3,  expression(paste("Frequency of sib matings (", italic(beta),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        


    # Panel D: r = 0.5
    r4  <-  dat[dat$r == 0.5,]
    d1  <-  r4[r4$h == 0.0,]
    d1  <-  d1[d1$sf == 0.01,]
    d2  <-  r4[r4$h == 0.0,]
    d2  <-  d2[d2$sf == 0.05,]
    d3  <-  r4[r4$h == 0.1,]
    d3  <-  d3[d3$sf == 0.01,]
    d4  <-  r4[r4$h == 0.1,]
    d4  <-  d4[d4$sf == 0.05,]
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(-0.005,0.001), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        abline(h=0, lwd=2)
        lines(invSelAdv ~ beta, lwd=1.5, lty=1, col=transparentColor('#252525', opacity=1), data=d1)
        points(invSelAdv ~ beta, pch=21, bg=transparentColor('grey60', opacity=1), data=d1)
        lines(invSelAdv ~ beta, lwd=1.5, lty=2, col=transparentColor('#252525', opacity=1), data=d2)
        points(invSelAdv ~ beta, pch=22, bg=transparentColor('grey60', opacity=1), data=d2)
        lines(invSelAdv ~ beta, lwd=1.5, lty=3, col=transparentColor('#252525', opacity=1), data=d3)
        points(invSelAdv ~ beta, pch=23, bg=transparentColor('grey60', opacity=1), data=d3)
        lines(invSelAdv ~ beta, lwd=1.5, lty=4, col=transparentColor('#252525', opacity=1), data=d4)
        points(invSelAdv ~ beta, pch=24, bg=transparentColor('grey60', opacity=1), data=d4)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.075, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(r), " = ", 0.5)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.3,  expression(paste("Frequency of sib matings (", italic(beta),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        

}







##############################################################
##############################################################
##  3-locus model figures

#' Heterozygote Advantage Figure -- 2 locus model
#' Reproducing Fig. 1a from Charlesworth & Wall (1999)
HetAdv3Locus  <-  function(df = "./data/3LocusHetAdvSimData_q0.5_r0.5.csv") {

    ## import data.frames
    dat  <-  read.csv(file = df)
    d1   <-  dat[dat$s1 == 0.1,]
    d1   <-  d1[d1$t1 == 0.1,]
    d2   <-  dat[dat$s1 == 0.2,]
    d2   <-  d2[d2$t1 == 0.2,]
    d3   <-  dat[dat$s1 == 0.1,]
    d3   <-  d3[d3$t1 == 0.2,]
    d4   <-  dat[dat$s1 == 0.2,]
    d4   <-  d4[d4$t1 == 0.1,]

    # set plot layout
    layout.mat <- matrix(c(1), nrow=1, ncol=1, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

    # make plot
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(4,5,5,2), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(0,0.08), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(invSelAdv ~ beta, lwd=1.5, lty=1, col=transparentColor('#252525', opacity=1), data=d1)
        points(invSelAdv ~ beta, pch=21, bg=transparentColor('#252525', opacity=0.75), data=d1)
        lines(invSelAdv ~ beta, lwd=1.5, lty=2, col=transparentColor('#252525', opacity=1), data=d2)
        points(invSelAdv ~ beta, pch=22, bg=transparentColor('#252525', opacity=0.75), data=d2)
        lines(invSelAdv[1:9] ~ beta[1:9], lwd=1.5, lty=3, col=transparentColor('#252525', opacity=1), data=d3)
        points(invSelAdv[1:9] ~ beta[1:9], pch=23, bg=transparentColor('#252525', opacity=0.75), data=d3)
        lines(invSelAdv[1:9] ~ beta[1:9], lwd=1.5, lty=4, col=transparentColor('#252525', opacity=1), data=d4)
        points(invSelAdv[1:9] ~ beta[1:9], pch=24, bg=transparentColor('#252525', opacity=0.75), data=d4)
        abline(h=0, lwd=2)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel( 0.5,  1.1,   expression(paste("Heterozygote Advantage, ", bold(B), " locus neutral")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.2,  0.5,   expression(paste(italic(s[I]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.2,  expression(paste("Frequency of sib matings (", italic(alpha),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2]*0.37,
               y       =  usr[4]*0.95,
               legend  =  c(
                            expression(paste(italic(s[1]), " = ", italic(t[1]), " = ", 0.1)),
                            expression(paste(italic(s[1]), " = ", italic(t[1]), " = ", 0.2)),
                            expression(paste(italic(s[1]), " = ", 0.1, "; ", italic(t[1]), " = ", 0.2)),
                            expression(paste(italic(s[1]), " = ", 0.2, "; ", italic(t[1]), " = ", 0.1))),
               pch     =  c(21,22,23,24),
               pt.bg      =  c(transparentColor('#252525', opacity=0.75),
                               transparentColor('#252525', opacity=0.75),
                               transparentColor('#252525', opacity=0.75),
                               transparentColor('#252525', opacity=0.75)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

}


#' Heterozygote Advantage Figure -- 2 locus model
#' Effect of recombination
HetAdv3LocusQR  <-  function(df1 = "./data/3LocusHetAdvSimData_q0.1_r0.5.csv",
                             df2 = "./data/3LocusHetAdvSimData_q0.05_r0.5.csv",
                             df3 = "./data/3LocusHetAdvSimData_q0.1_r0.05.csv",
                             df4 = "./data/3LocusHetAdvSimData_q0.05_r0.05.csv") {

    ## import data.frames
    dat  <-  read.csv(file = df1)
    d1   <-  dat[dat$s1 == 0.05,]
    d1   <-  d1[d1$t1 == 0.05,]
    d2   <-  dat[dat$s1 == 0.1,]
    d2   <-  d2[d2$t1 == 0.1,]
    d3   <-  dat[dat$s1 == 0.05,]
    d3   <-  d3[d3$t1 == 0.1,]
    d4   <-  dat[dat$s1 == 0.1,]
    d4   <-  d4[d4$t1 == 0.05,]

    # set plot layout
    layout.mat <- matrix(c(1,2,3,4), nrow=2, ncol=2, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

    # make plot
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(4,5,5,2), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(-0.01,0.08), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(invSelAdv ~ beta, lwd=1.5, lty=1, col=transparentColor('#252525', opacity=1), data=d1)
        points(invSelAdv ~ beta, pch=21, bg=transparentColor('#252525', opacity=0.75), data=d1)
        lines(invSelAdv ~ beta, lwd=1.5, lty=2, col=transparentColor('#252525', opacity=1), data=d2)
        points(invSelAdv ~ beta, pch=22, bg=transparentColor('#252525', opacity=0.75), data=d2)
        lines(invSelAdv[1:9] ~ beta[1:9], lwd=1.5, lty=3, col=transparentColor('#252525', opacity=1), data=d3)
        points(invSelAdv[1:9] ~ beta[1:9], pch=23, bg=transparentColor('#252525', opacity=0.75), data=d3)
        lines(invSelAdv[1:9] ~ beta[1:9], lwd=1.5, lty=4, col=transparentColor('#252525', opacity=1), data=d4)
        points(invSelAdv[1:9] ~ beta[1:9], pch=24, bg=transparentColor('#252525', opacity=0.75), data=d4)
        abline(h=0, lwd=2)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel( 1.4,  1.5,   expression(paste("Heterozygote Advantage -- w/ linkage")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(q), "=", 0.1, "; ",italic(r), "=", 0.5)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.2,  0.5,   expression(paste(italic(s[I]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.2,  expression(paste("Frequency of sib matings (", italic(alpha),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2]*0.65,
               y       =  usr[4]*0.95,
               legend  =  c(
                            expression(paste(italic(s[1]), " = ", italic(t[1]), " = ",italic(s[2]), " = ", italic(t[2]), " = ", 0.05)),
                            expression(paste(italic(s[1]), " = ", italic(t[1]), " = ",italic(s[2]), " = ", italic(t[2]), " = ", 0.1)),
                            expression(paste(italic(s[1]), " = ", italic(s[2]), " = ", 0.05, "; ", italic(t[1]), " = ", italic(t[2]), " = ", 0.1)),
                            expression(paste(italic(s[1]), " = ", italic(s[2]), " = ", 0.1, "; ", italic(t[1]), " = ", italic(t[2]), " = ", 0.05))),
               pch     =  c(21,22,23,24),
               pt.bg      =  c(transparentColor('#252525', opacity=0.75),
                               transparentColor('#252525', opacity=0.75),
                               transparentColor('#252525', opacity=0.75),
                               transparentColor('#252525', opacity=0.75)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

    # Panel B -- 
    ## import data.frames
    dat  <-  read.csv(file = df2)
    d1   <-  dat[dat$s1 == 0.05,]
    d1   <-  d1[d1$t1 == 0.05,]
    d2   <-  dat[dat$s1 == 0.1,]
    d2   <-  d2[d2$t1 == 0.1,]
    d3   <-  dat[dat$s1 == 0.05,]
    d3   <-  d3[d3$t1 == 0.1,]
    d4   <-  dat[dat$s1 == 0.1,]
    d4   <-  d4[d4$t1 == 0.05,]
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(-0.01,0.08), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(invSelAdv ~ beta, lwd=1.5, lty=1, col=transparentColor('#252525', opacity=1), data=d1)
        points(invSelAdv ~ beta, pch=21, bg=transparentColor('#252525', opacity=0.75), data=d1)
        lines(invSelAdv ~ beta, lwd=1.5, lty=2, col=transparentColor('#252525', opacity=1), data=d2)
        points(invSelAdv ~ beta, pch=22, bg=transparentColor('#252525', opacity=0.75), data=d2)
        lines(invSelAdv[1:9] ~ beta[1:9], lwd=1.5, lty=3, col=transparentColor('#252525', opacity=1), data=d3)
        points(invSelAdv[1:9] ~ beta[1:9], pch=23, bg=transparentColor('#252525', opacity=0.75), data=d3)
        lines(invSelAdv[1:9] ~ beta[1:9], lwd=1.5, lty=4, col=transparentColor('#252525', opacity=1), data=d4)
        points(invSelAdv[1:9] ~ beta[1:9], pch=24, bg=transparentColor('#252525', opacity=0.75), data=d4)
        abline(h=0, lwd=2)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(q), "=", 0.05, "; ",italic(r), "=", 0.5)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(-0.2,  0.5,   expression(paste(italic(s[I]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.2,  expression(paste("Frequency of sib matings (", italic(alpha),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        

    # Panel C -- 
    ## import data.frames
    dat  <-  read.csv(file = df3)
    d1   <-  dat[dat$s1 == 0.05,]
    d1   <-  d1[d1$t1 == 0.05,]
    d2   <-  dat[dat$s1 == 0.1,]
    d2   <-  d2[d2$t1 == 0.1,]
    d3   <-  dat[dat$s1 == 0.05,]
    d3   <-  d3[d3$t1 == 0.1,]
    d4   <-  dat[dat$s1 == 0.1,]
    d4   <-  d4[d4$t1 == 0.05,]
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(-0.01,0.08), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(invSelAdv ~ beta, lwd=1.5, lty=1, col=transparentColor('#252525', opacity=1), data=d1)
        points(invSelAdv ~ beta, pch=21, bg=transparentColor('#252525', opacity=0.75), data=d1)
        lines(invSelAdv ~ beta, lwd=1.5, lty=2, col=transparentColor('#252525', opacity=1), data=d2)
        points(invSelAdv ~ beta, pch=22, bg=transparentColor('#252525', opacity=0.75), data=d2)
        lines(invSelAdv[1:9] ~ beta[1:9], lwd=1.5, lty=3, col=transparentColor('#252525', opacity=1), data=d3)
        points(invSelAdv[1:9] ~ beta[1:9], pch=23, bg=transparentColor('#252525', opacity=0.75), data=d3)
        lines(invSelAdv[1:9] ~ beta[1:9], lwd=1.5, lty=4, col=transparentColor('#252525', opacity=1), data=d4)
        points(invSelAdv[1:9] ~ beta[1:9], pch=24, bg=transparentColor('#252525', opacity=0.75), data=d4)
        abline(h=0, lwd=2)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(q), "=", 0.1, "; ",italic(r), "=", 0.05)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.2,  0.5,   expression(paste(italic(s[I]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.2,  expression(paste("Frequency of sib matings (", italic(alpha),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        

    # Panel C -- 
    ## import data.frames
    dat  <-  read.csv(file = df4)
    d1   <-  dat[dat$s1 == 0.05,]
    d1   <-  d1[d1$t1 == 0.05,]
    d2   <-  dat[dat$s1 == 0.1,]
    d2   <-  d2[d2$t1 == 0.1,]
    d3   <-  dat[dat$s1 == 0.05,]
    d3   <-  d3[d3$t1 == 0.1,]
    d4   <-  dat[dat$s1 == 0.1,]
    d4   <-  d4[d4$t1 == 0.05,]
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(-0.01,0.08), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(invSelAdv ~ beta, lwd=1.5, lty=1, col=transparentColor('#252525', opacity=1), data=d1)
        points(invSelAdv ~ beta, pch=21, bg=transparentColor('#252525', opacity=0.75), data=d1)
        lines(invSelAdv ~ beta, lwd=1.5, lty=2, col=transparentColor('#252525', opacity=1), data=d2)
        points(invSelAdv ~ beta, pch=22, bg=transparentColor('#252525', opacity=0.75), data=d2)
        lines(invSelAdv[1:9] ~ beta[1:9], lwd=1.5, lty=3, col=transparentColor('#252525', opacity=1), data=d3)
        points(invSelAdv[1:9] ~ beta[1:9], pch=23, bg=transparentColor('#252525', opacity=0.75), data=d3)
        lines(invSelAdv[1:9] ~ beta[1:9], lwd=1.5, lty=4, col=transparentColor('#252525', opacity=1), data=d4)
        points(invSelAdv[1:9] ~ beta[1:9], pch=24, bg=transparentColor('#252525', opacity=0.75), data=d4)
        abline(h=0, lwd=2)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(q), "=", 0.05, "; ",italic(r), "=", 0.05)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(-0.2,  0.5,   expression(paste(italic(s[I]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.2,  expression(paste("Frequency of sib matings (", italic(alpha),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
}




sheltering3LocusSuppFig  <-  function() {

    ## import data.frames
    dat    <-  read.csv(file = "./data/3LocusHetShelteringData_qrhsCompare.csv")
    
    # Subset data for plotting
    Qs  <-  unique(dat$q)
    Rs  <-  unique(dat$r)
    Hs  <-  unique(dat$hfA) 
    # unique(dat$hmA) == unique(dat$hmB)
    q1    <-  dat[dat$q == Qs[1],]
    q1r1  <-  q1[q1$r == Rs[1],]
    d1    <-  q1r1[q1r1$hfA == Hs[1],]
    d2    <-  q1r1[q1r1$hfA == Hs[2],]

 
    # set plot layout
    layout.mat <- matrix(c(1:4), nrow=2, ncol=2, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

    # Panel A
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(-0.003,0.001), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        abline(h=0, lwd=2)
        lines(invSelAdv ~ beta, lwd=1.5, lty=1, col=transparentColor('#252525', opacity=1), data=d1)
        points(invSelAdv ~ beta, pch=22, bg=transparentColor('grey60', opacity=1), data=d1)
        lines(invSelAdv ~ beta, lwd=1.5, lty=1, col=transparentColor('#252525', opacity=1), data=d2)
        points(invSelAdv ~ beta, pch=21, bg=transparentColor('grey60', opacity=1), data=d2)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(q)==0.5, ", ", italic(r)==0.5, "")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(paste(italic(s[I]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(paste(italic(h[i]), " = 0.0")),
                            expression(paste(italic(h[i]), " = 0.1"))),
               pch     =  c(21,22),
               pt.bg      =  c(transparentColor('grey60', opacity=1),
                               transparentColor('grey60', opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


    # Panel B: r = 0.1
    # Subset data for plotting
    q1r2  <-  q1[q1$r == Rs[2],]
    d1    <-  q1r2[q1r2$hfA == Hs[1],]
    d2    <-  q1r2[q1r2$hfA == Hs[2],]

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(-0.003,0.001), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        abline(h=0, lwd=2)
        lines(invSelAdv ~ beta, lwd=1.5, lty=1, col=transparentColor('#252525', opacity=1), data=d1)
        points(invSelAdv ~ beta, pch=22, bg=transparentColor('grey60', opacity=1), data=d1)
        lines(invSelAdv ~ beta, lwd=1.5, lty=1, col=transparentColor('#252525', opacity=1), data=d2)
        points(invSelAdv ~ beta, pch=21, bg=transparentColor('grey60', opacity=1), data=d2)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.075, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(q)==0.5, ", ", italic(r)==0.001, "")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)

    # Panel C
    q2    <-  dat[dat$q == Qs[2],]
    q2r1  <-  q2[q2$r == Rs[1],]
    d1    <-  q2r1[q2r1$hfA == Hs[1],]
    d2    <-  q2r1[q2r1$hfA == Hs[2],]
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(-0.003,0.001), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        abline(h=0, lwd=2)
        lines(invSelAdv ~ beta, lwd=1.5, lty=1, col=transparentColor('#252525', opacity=1), data=d1)
        points(invSelAdv ~ beta, pch=22, bg=transparentColor('grey60', opacity=1), data=d1)
        lines(invSelAdv ~ beta, lwd=1.5, lty=1, col=transparentColor('#252525', opacity=1), data=d2)
        points(invSelAdv ~ beta, pch=21, bg=transparentColor('grey60', opacity=1), data=d2)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.075, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(q)==0.001, ", ", italic(r)==0.5, "")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(paste(italic(s[I]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5, -0.3,  expression(paste("Frequency of sib matings (", italic(alpha),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        


    # Panel D: r = 0.5
    q2r2  <-  q2[q2$r == Rs[2],]
    d1    <-  q2r2[q2r2$hfA == Hs[1],]
    d2    <-  q2r2[q2r2$hfA == Hs[2],]
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(-0.003,0.001), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        abline(h=0, lwd=2)
        lines(invSelAdv ~ beta, lwd=1.5, lty=1, col=transparentColor('#252525', opacity=1), data=d1)
        points(invSelAdv ~ beta, pch=22, bg=transparentColor('grey60', opacity=1), data=d1)
        lines(invSelAdv ~ beta, lwd=1.5, lty=1, col=transparentColor('#252525', opacity=1), data=d2)
        points(invSelAdv ~ beta, pch=21, bg=transparentColor('grey60', opacity=1), data=d2)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.075, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(q)==0.001, ", ", italic(r)==0.001, "")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.3,  expression(paste("Frequency of sib matings (", italic(alpha),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        

}






sheltering3LocusFig  <-  function() {

    ## import data.frames
    dat    <-  read.csv(file = "./data/3LocusHetShelteringData_u1e-05_v1e-05.csv")

    q1  <-  dat[dat$q == 0.5,]
    d1  <-  q1[q1$hfA == 0.0,]
    d1  <-  d1[d1$sfA == 0.01,]
    d2  <-  q1[q1$hfA == 0.0,]
    d2  <-  d2[d2$sfA == 0.05,]
    d3  <-  q1[q1$hfA == 0.1,]
    d3  <-  d3[d3$sfA == 0.01,]
    d4  <-  q1[q1$hfA == 0.1,]
    d4  <-  d4[d4$sfA == 0.05,]

    # set plot layout
    layout.mat <- matrix(c(1:4), nrow=2, ncol=2, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

    # Panel A
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(-0.006,0.001), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        abline(h=0, lwd=2)
        lines(invSelAdv ~ beta, lwd=1.5, lty=1, col=transparentColor('#252525', opacity=1), data=d1)
        points(invSelAdv ~ beta, pch=21, bg=transparentColor('grey60', opacity=1), data=d1)
        lines(invSelAdv ~ beta, lwd=1.5, lty=2, col=transparentColor('#252525', opacity=1), data=d2)
        points(invSelAdv ~ beta, pch=22, bg=transparentColor('grey60', opacity=1), data=d2)
        lines(invSelAdv ~ beta, lwd=1.5, lty=3, col=transparentColor('#252525', opacity=1), data=d3)
        points(invSelAdv ~ beta, pch=23, bg=transparentColor('grey60', opacity=1), data=d3)
        lines(invSelAdv ~ beta, lwd=1.5, lty=4, col=transparentColor('#252525', opacity=1), data=d4)
        points(invSelAdv ~ beta, pch=24, bg=transparentColor('grey60', opacity=1), data=d4)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, at=c(0.001, 0, -0.002, -0.004, -0.006))
        # Plot labels etc.
        proportionalLabel(0.04, 1.075, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(q), " = ", 0.5)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(paste(italic(s[I]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  -usr[4]*1.25,
               legend  =  c(
                            expression(paste(italic(h[i]), " = 0.0; ", italic(s[i]), " = ", 0.01)),
                            expression(paste(italic(h[i]), " = 0.0; ", italic(s[i]), " = ", 0.05)),
                            expression(paste(italic(h[i]), " = ", 0.1, "; ", italic(s[i]), " = ", 0.01)),
                            expression(paste(italic(h[i]), " = ", 0.1, "; ", italic(s[i]), " = ", 0.05))),
               pch     =  c(21,22,23,24),
               pt.bg      =  c(transparentColor('grey60', opacity=1),
                               transparentColor('grey60', opacity=1),
                               transparentColor('grey60', opacity=1),
                               transparentColor('grey60', opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


    # Panel B: r = 0.1
    q2  <-  dat[dat$q == 0.1,]
    d1  <-  q2[q2$hfA == 0.0,]
    d1  <-  d1[d1$sfA == 0.01,]
    d2  <-  q2[q2$hfA == 0.0,]
    d2  <-  d2[d2$sfA == 0.05,]
    d3  <-  q2[q2$hfA == 0.1,]
    d3  <-  d3[d3$sfA == 0.01,]
    d4  <-  q2[q2$hfA == 0.1,]
    d4  <-  d4[d4$sfA == 0.05,]
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(-0.006,0.001), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        abline(h=0, lwd=2)
        lines(invSelAdv ~ beta, lwd=1.5, lty=1, col=transparentColor('#252525', opacity=1), data=d1)
        points(invSelAdv ~ beta, pch=21, bg=transparentColor('grey60', opacity=1), data=d1)
        lines(invSelAdv ~ beta, lwd=1.5, lty=2, col=transparentColor('#252525', opacity=1), data=d2)
        points(invSelAdv ~ beta, pch=22, bg=transparentColor('grey60', opacity=1), data=d2)
        lines(invSelAdv ~ beta, lwd=1.5, lty=3, col=transparentColor('#252525', opacity=1), data=d3)
        points(invSelAdv ~ beta, pch=23, bg=transparentColor('grey60', opacity=1), data=d3)
        lines(invSelAdv ~ beta, lwd=1.5, lty=4, col=transparentColor('#252525', opacity=1), data=d4)
        points(invSelAdv ~ beta, pch=24, bg=transparentColor('grey60', opacity=1), data=d4)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, at=c(0.001, 0, -0.002, -0.004, -0.006), labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.075, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(q), " = ", 0.1)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)

    # Panel C
    q3  <-  dat[dat$q == 0.01,]
    d1  <-  q3[q3$hfA == 0.0,]
    d1  <-  d1[d1$sfA == 0.01,]
    d2  <-  q3[q3$hfA == 0.0,]
    d2  <-  d2[d2$sfA == 0.05,]
    d3  <-  q3[q3$hfA == 0.1,]
    d3  <-  d3[d3$sfA == 0.01,]
    d4  <-  q3[q3$hfA == 0.1,]
    d4  <-  d4[d4$sfA == 0.05,]
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(-0.006,0.001), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        abline(h=0, lwd=2)
        lines(invSelAdv ~ beta, lwd=1.5, lty=1, col=transparentColor('#252525', opacity=1), data=d1)
        points(invSelAdv ~ beta, pch=21, bg=transparentColor('grey60', opacity=1), data=d1)
        lines(invSelAdv ~ beta, lwd=1.5, lty=2, col=transparentColor('#252525', opacity=1), data=d2)
        points(invSelAdv ~ beta, pch=22, bg=transparentColor('grey60', opacity=1), data=d2)
        lines(invSelAdv ~ beta, lwd=1.5, lty=3, col=transparentColor('#252525', opacity=1), data=d3)
        points(invSelAdv ~ beta, pch=23, bg=transparentColor('grey60', opacity=1), data=d3)
        lines(invSelAdv ~ beta, lwd=1.5, lty=4, col=transparentColor('#252525', opacity=1), data=d4)
        points(invSelAdv ~ beta, pch=24, bg=transparentColor('grey60', opacity=1), data=d4)
        # axes
        axis(1, las=1)
        axis(2, las=1, at=c(0.001, 0, -0.002, -0.004, -0.006))
        # Plot labels etc.
        proportionalLabel(0.04, 1.075, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(q), " = ", 0.01)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(paste(italic(s[I]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5, -0.3,  expression(paste("Frequency of sib matings (", italic(beta),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        


    # Panel D: r = 0.5
    q4  <-  dat[dat$q == 0.001,]
    d1  <-  q4[q4$hfA == 0.0,]
    d1  <-  d1[d1$sfA == 0.01,]
    d2  <-  q4[q4$hfA == 0.0,]
    d2  <-  d2[d2$sfA == 0.05,]
    d3  <-  q4[q4$hfA == 0.1,]
    d3  <-  d3[d3$sfA == 0.01,]
    d4  <-  q4[q4$hfA == 0.1,]
    d4  <-  d4[d4$sfA == 0.05,]
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(-0.006,0.001), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        abline(h=0, lwd=2)
        lines(invSelAdv  ~ beta, lwd=1.5, lty=1, col=transparentColor('#252525', opacity=1), data=d1)
        points(invSelAdv ~ beta, pch=21, bg=transparentColor('grey60', opacity=1), data=d1)
        lines(invSelAdv  ~ beta, lwd=1.5, lty=2, col=transparentColor('#252525', opacity=1), data=d2)
        points(invSelAdv ~ beta, pch=22, bg=transparentColor('grey60', opacity=1), data=d2)
        lines(invSelAdv  ~ beta, lwd=1.5, lty=3, col=transparentColor('#252525', opacity=1), data=d3)
        points(invSelAdv ~ beta, pch=23, bg=transparentColor('grey60', opacity=1), data=d3)
        lines(invSelAdv  ~ beta, lwd=1.5, lty=4, col=transparentColor('#252525', opacity=1), data=d4)
        points(invSelAdv ~ beta, pch=24, bg=transparentColor('grey60', opacity=1), data=d4)
        # axes
        axis(1, las=1)
        axis(2, las=1, at=c(0.001, 0, -0.002, -0.004, -0.006), labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.075, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(q), " = ", 0.001)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.3,  expression(paste("Frequency of sib matings (", italic(beta),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        

}




sheltering3Locus_hsCompareFig  <-  function() {

    ## import data.frames
    dat    <-  read.csv(file = "./data/3LocusHetShelteringData_q0.5_r0.005_hsCompare.csv")

    h.01  <-  dat[dat$hfA == 0.01,]
    d1    <-  h.01[h.01$sfA == 0.01,]
    d2    <-  h.01[h.01$sfA == 0.05,]
    h.25  <-  dat[dat$hfA == 0.25,]
    d3    <-  h.25[h.25$sfA == 0.01,]
    d4    <-  h.25[h.25$sfA == 0.05,]

    # Panel A
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(5,5,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,1), ylim = c(-0.015,0.0045), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        abline(h=0, lwd=2)
        lines(invSelAdv ~ beta, lwd=1.5, lty=1, col=transparentColor('#252525', opacity=1), data=d1)
        points(invSelAdv ~ beta, pch=21, bg=transparentColor('grey60', opacity=1), data=d1)
        lines(invSelAdv ~ beta, lwd=1.5, lty=2, col=transparentColor('#252525', opacity=1), data=d2)
        points(invSelAdv ~ beta, pch=22, bg=transparentColor('grey60', opacity=1), data=d2)
        lines(invSelAdv ~ beta, lwd=1.5, lty=3, col=transparentColor('#252525', opacity=1), data=d3)
        points(invSelAdv ~ beta, pch=23, bg=transparentColor('grey60', opacity=1), data=d3)
        lines(invSelAdv ~ beta, lwd=1.5, lty=4, col=transparentColor('#252525', opacity=1), data=d4)
        points(invSelAdv ~ beta, pch=24, bg=transparentColor('grey60', opacity=1), data=d4)
        # axes
        axis(1, las=1)
        axis(2, las=1)#, at=c(0.001, 0, -0.002, -0.004, -0.006))
        # Plot labels etc.
        proportionalLabel( 0.5,  1.1,   expression(paste(italic(q), " = 0.5, ", italic(r), " = 0.005")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.2,  0.5,   expression(paste(italic(s[I]))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.15,  expression(paste("Frequency of sib matings (", italic(alpha),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(paste(italic(h), " = 0.01; ", italic(s), " = ", 0.01)),
                            expression(paste(italic(h), " = 0.01; ", italic(s), " = ", 0.05)),
                            expression(paste(italic(h), " = 0.25; ", italic(s), " = ", 0.01)),
                            expression(paste(italic(h), " = 0.25; ", italic(s), " = ", 0.05))),
               pch     =  c(21,22,23,24),
               pt.bg      =  c(transparentColor('grey60', opacity=1),
                               transparentColor('grey60', opacity=1),
                               transparentColor('grey60', opacity=1),
                               transparentColor('grey60', opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


}
