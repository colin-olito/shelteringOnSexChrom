#  Functions to create figures for inversions expanding SLR 
#  in the presence of deleterious mutational variation
source('./R/functions-figures.R')
source('./R/simulations-inversions-delMut-RECODE.R')
source('./R/simulations-inversions-delMut-RECODE-Nei.R')
source('./R/Auto-Inversions-NeiMyVersion.R')



#' General solution for the frequency trajectory 
#' of an initially unloaded SLR-expanding inversion
YI.genSol  <-  function(YI.init, U, x, s, h, t) {
     (exp(((1 - exp(-h*s*t))*U*x) / (1 - exp(-h*s)))*YI.init)/(1 - YI.init + exp(((1 - exp(-h*s*t))*U*x) / (1 - exp(-h*s)))*YI.init)
 }

#' Effective initial frequency of unloaded inveresions.
#' Equal to the limit of YI.genSol as t --> infinity
qstar  <-  function(q.0, U, x, s, h) {
     (q.0*exp((U*x)/(1 - exp(-s*h)))) / (1 - q.0*(1 - exp((U*x)/(1 - exp(-h*s)))))
 }



#' FIGURE 1. Illustration of deterministic dynamics
#'           for different dominance values
deterministicDominanceIllustration  <-  function(wHap = FALSE, x=0.2) {

    # set plot layout
    layout.mat <- matrix(c(1:12), nrow=4, ncol=3, byrow=FALSE)
    layout     <- layout(layout.mat,respect=TRUE)

    # Make data for plotting 
    s        <-  0.01
    h        <-  0.25
    nTot     <-  10^4 
    U        <-  0.02
    u        <-  U/nTot
    x        <-  0.2

    ###########
    # h = 0.25

    # make data for plotting
    r.0.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 0, x = x, h = h, generations = 2500)
    r.1.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 1, x = x, h = h, generations = 2500)
    r.0.rel.w.YI.t  <-  r.0.Dat$rel.w.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.1.rel.w.YI.t  <-  r.1.Dat$rel.w.YI.t[round(r.1.Dat$YI.t,digits=4) > 0]

    # Old-version
    r.0.DatOld  <-  makeDeterministicFigSimData(U = U, s = s, r = 0, x = x, h = h, generations = 2500)
    r.1.DatOld  <-  makeDeterministicFigSimData(U = U, s = s, r = 1, x = x, h = h, generations = 2500)
    r.0.rel.w.YI.t.Nei  <-  r.0.DatOld$rel.w.YI.t[round(r.0.DatOld$YI.t,digits=4) < 1]
    r.1.rel.w.YI.t.Nei  <-  r.1.DatOld$rel.w.YI.t[round(r.1.DatOld$YI.t,digits=4) > 0]

    # Haploid Model
    r.0.Hap.Dat  <-  makeDeterministicSimHaploidData(U = U, s = s, r = 0, x = x, h = h, generations = 2500)
    r.1.Hap.Dat  <-  makeDeterministicSimHaploidData(U = U, s = s, r = 1, x = x, h = h, generations = 2500)


    # Colors
    COLS  <- colorRampPalette(wes_palette("Zissou1"))(4)

    # Axes limits
    relFitMin  <-  0.996
    relFitMax  <-  1.005


## Panel A
# Effect of initial del. mut. load on inversion relative fitness 
    # make plot
    par(omi=c(0.5, 0.25, 0.5, 0.5), mar = c(3,6,1,0), bty='o', xaxt='s', yaxt='s')
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(relFitMin, relFitMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Inversion relative fitness
        lines(r.0.rel.w.YI.t ~ seq_along(r.0.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1))
        lines(r.1.rel.w.YI.t ~ seq_along(r.1.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1))
        # Inversion relative fitness -- Haploid model
        if(wHap) {
            lines(rel.w.YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.0.Hap.Dat)
            lines(rel.w.YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.1.Hap.Dat)
        }
        # Benchmark for relative fitness of 1
        abline(h=1, lwd=1, lty=2, col=1)
        # Points depicting (arbitrary) extinction frequency of 0.00001
        points(r.1.rel.w.YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        # Annotate 1 + Ux approximation
        points((1 + U*x) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        proportionalArrows(0.127, 0.86, 0.06, 0.86, cex=1, adj=c(0.5, 0.5), xpd=NA, log=FALSE, length=0.05)
        proportionalLabel(0.275, 0.86, expression(1+italic(Ux)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, expression(paste(italic(Ux/(hs)==1.6))), cex=2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.45,  0.5,   expression(paste("Inversion Relative Fitness")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(d)==0),
                            expression(italic(d)==1)),
               lty     =  1,
               lwd     =  2,
               col     =  c(transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

qstr       <-  qstar(U=U, q.0=r.0.Dat$YI.t[1], x=x, s=s, h=h)
t          <-  seq_along(r.0.rel.w.YI.t)
YI.gensol  <-  YI.genSol(U=U, YI.init=r.0.Dat$YI.t[1], x=x, s=s, h=h, t=t)

## Panel D
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,(qstr*1.25)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Inversion Frequencies
        lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(YI.gensol ~ t, lwd=2, lty=3)
        if(wHap) {
            lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.0.Hap.Dat)
            lines(YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.1.Hap.Dat)
        }
        # Points indicating inversion extinctions
        points(YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        # benchmark qstar
        abline(h=qstr, lty=2)
        proportionalLabel(0.9, 0.83, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.45,  0.5,   expression(paste("Inversion Frequency (", italic(Y[I]), ")")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        # Legend
        legend(
               x       =  usr[2]*0.95,
               y       =  usr[4],
               legend  =  c(
                            expression(paste("Gen. sol. for ", italic(Y[I]), " when ", italic(d)==0))),
               lty     =  3,
               lwd     =  2,
               col     =  1,
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

    # reset Colors
     COLS  <-  c('red', '#252525')

## Panel G
# deleterious allele frequencies for r = 1
     # Make plot
    delMax  <-  r.1.Dat$qt.Xf.W[1]*1.05
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines for deleterious allele frequencies
        lines(qt.Xf.W[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.Xf.W[1:length(r.1.rel.w.YI.t)]), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qt.Xf.D[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.Xf.D[1:length(r.1.rel.w.YI.t)]), lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qt.YI.W[1:length(r.1.rel.w.YI.t)]/YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.YI.W[1:length(r.1.rel.w.YI.t)]), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
        # Points indicating inversion extinctions
        points(qt.Xf.W[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(qt.Xf.D[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(qt.YI.W[length(r.1.rel.w.YI.t)]/YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.1.Dat)
        # axes
        axis(1, las=1, cex=1.3, labels=NA)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'G', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(d)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.45,  0.5,   expression(paste("Deleterious Allele Frequency")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        
#        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend = c(expression(italic(q)[italic(X[f])]^italic(W)), 
                          expression(italic(q)[italic(Y[I])]^italic(W)), 
                          expression(italic(q)[italic(X[f])]^italic(D))),
               lty     =  c(1,1,3),
               lwd     =  2,
               col     =  c(transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


## Panel J
# deleterious allele frequencies for r = 1
     # Make plot
    delMax  <-  max(r.1.Dat$qt.Xf.W)
    delMin  <-  min(r.1.Dat$qt.Xf.D)
    delRange  <-  delMax - delMin
    delMax  <-  delMax + 0.05*delRange
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(delMin,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines for deleterious allele frequencies
        lines(qt.Xf.W[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.Xf.W[1:length(r.1.rel.w.YI.t)]), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qt.Xf.D[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.Xf.D[1:length(r.1.rel.w.YI.t)]), lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
#        lines(qt.YI.W[1:length(r.1.rel.w.YI.t)]/YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.YI.W[1:length(r.1.rel.w.YI.t)]), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
        # Points indicating inversion extinctions
        points(qt.Xf.W[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(qt.Xf.D[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
#        points(qt.YI.W[length(r.1.rel.w.YI.t)]/YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.1.Dat)
        # axes
        axis(1, las=1, cex=1.3)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'J', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(d)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.45,  0.5,   expression(paste("Deleterious Allele Frequency")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend = c(expression(italic(q)[italic(X[f])]^italic(W)), 
#                          expression(italic(q)[italic(Y[I])]^italic(W)), 
                          expression(italic(q)[italic(X[f])]^italic(D))),
               lty     =  c(1,3),
               lwd     =  2,
               col     =  c(transparentColor(COLS[2], opacity=1),
#                            transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )



    ###########
    # h = 0.1
    # make data for plotting
    h  <-  0.1

    # Nei-style model
    r.0.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 0, x = x, h = h, generations = 8000)
    r.1.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 1, x = x, h = h, generations = 8000)
    r.2.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 2, x = x, h = h, generations = 8000)
    r.3.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 3, x = x, h = h, generations = 8000)

    r.0.rel.w.YI.t  <-  r.0.Dat$rel.w.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.1.rel.w.YI.t  <-  r.1.Dat$rel.w.YI.t[round(r.1.Dat$YI.t,digits=4) > 0]
    r.2.rel.w.YI.t  <-  r.2.Dat$rel.w.YI.t[round(r.2.Dat$YI.t,digits=4) > 0]
    r.3.rel.w.YI.t  <-  r.3.Dat$rel.w.YI.t[round(r.3.Dat$YI.t,digits=4) > 0]

    r.0.Hap.Dat  <-  makeDeterministicSimHaploidData(U = U, s = s, r = 0, x = x, h = h, generations = 8000)
    r.1.Hap.Dat  <-  makeDeterministicSimHaploidData(U = U, s = s, r = 1, x = x, h = h, generations = 8000)
    r.2.Hap.Dat  <-  makeDeterministicSimHaploidData(U = U, s = s, r = 2, x = x, h = h, generations = 8000)
    r.3.Hap.Dat  <-  makeDeterministicSimHaploidData(U = U, s = s, r = 3, x = x, h = h, generations = 8000)

   # Colors
    COLS     <-  colorRampPalette(wes_palette("Zissou1"))(4)

## Panel B
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(relFitMin,relFitMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Inversion relative fitnesses
        lines(r.3.rel.w.YI.t ~ seq_along(r.3.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1))
        lines(r.2.rel.w.YI.t ~ seq_along(r.2.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.1.rel.w.YI.t ~ seq_along(r.1.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r.0.rel.w.YI.t ~ seq_along(r.0.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1))
        # Inversion relative fitness -- Haploid model
        if(wHap) {
            lines(rel.w.YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.0.Hap.Dat)
            lines(rel.w.YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.1.Hap.Dat)
            lines(rel.w.YI.t[1:length(r.2.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.2.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.2.Hap.Dat)
            lines(rel.w.YI.t[1:length(r.3.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.3.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.3.Hap.Dat)
        }
        # Benchmark inversion relative fitness of 1
        abline(h=1, lwd=1, lty=2, col=1)
        # Points depicting (arbitrary) extinction frequency of 0.00001
        points(r.1.rel.w.YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r.2.rel.w.YI.t[length(r.2.rel.w.YI.t)] ~ length(r.2.rel.w.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
        points(r.3.rel.w.YI.t[length(r.3.rel.w.YI.t)] ~ length(r.3.rel.w.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        # Annotate 1 + Ux approximation
        points((1 + U*x) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, expression(paste(italic(Ux/(hs)==4))), cex=2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(d)==0),
                            expression(italic(d)==1),
                            expression(italic(d)==2),
                            expression(italic(d)==3)),
               lty     =  1,
               lwd     =  2,
               col     =  c(transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[4], opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )



qstr       <-  qstar(U=U, q.0=r.0.Dat$YI.t[1], x=x, s=s, h=h)
t          <-  seq_along(r.0.rel.w.YI.t)
YI.gensol  <-  YI.genSol(U=U, YI.init=r.0.Dat$YI.t[1], x=x, s=s, h=h, t=t)

## Panel E
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,(qstr*1.25)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Inversion frequencies
        lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(YI.t[1:length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.2.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.2.Dat)
        lines(YI.t[1:length(r.3.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.3.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.3.Dat)
        lines(YI.gensol ~ t, lwd=2, lty=3)
        if(wHap) {
            lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.0.Hap.Dat)
            lines(YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.1.Hap.Dat)
            lines(YI.t[1:length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.2.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.2.Hap.Dat)
            lines(YI.t[1:length(r.3.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.3.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.3.Hap.Dat)
        }
        # Points
        points(YI.t[length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.1.rel.w.YI.t)], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(YI.t[length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.2.rel.w.YI.t)], pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.2.Dat)
        points(YI.t[length(r.3.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.3.rel.w.YI.t)], pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.3.Dat)
        # benchmark qstar
        abline(h=qstr, lty=2)
        proportionalLabel(0.9, 0.83, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'E', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # reset Colors
     COLS  <-  c('red', '#252525')

## Panel H
# deleterious allele frequencies for r = 1
     # Make plot
    delMax  <-  r.1.Dat$qt.Xf.W[1]*1.05
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines for deleterious allele frequencies
        lines(qt.Xf.W[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.Xf.W[1:length(r.1.rel.w.YI.t)]), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qt.Xf.D[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.Xf.D[1:length(r.1.rel.w.YI.t)]), lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qt.YI.W[1:length(r.1.rel.w.YI.t)]/YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.YI.W[1:length(r.1.rel.w.YI.t)]), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
        # Points indicating inversion extinctions
        points(qt.Xf.W[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(qt.Xf.D[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(qt.YI.W[length(r.1.rel.w.YI.t)]/YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.1.Dat)
        # axes
        axis(1, las=1, cex=1.3, labels=NA)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'H', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(d)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(-0.3,  0.5,   expression(paste("Deleterious Allele Frequency")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        
#        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        

     # Make plot
    delMax  <-  max(r.1.Dat$qt.Xf.W)
    delMin  <-  min(r.1.Dat$qt.Xf.D)
    delRange  <-  delMax - delMin
    delMax  <-  delMax + 0.05*delRange

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(delMin,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Deleterious allele frequencies
        lines(qt.Xf.W[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.Xf.W)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qt.Xf.D[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.Xf.D)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
#        lines(qt.YI.W[1:length(r.1.rel.w.YI.t)]/YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.YI.W)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
        # Points indicating inversion extinctions
        points(qt.Xf.W[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(qt.Xf.D[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
#        points(qt.YI.W[length(r.1.rel.w.YI.t)]/YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.1.Dat)
        # axes
        axis(1, las=1, cex=1.2)
        axis(2, las=1, cex=1.2)        
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'K', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(d)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        


    ###########
    # h = 0.05

    # make data for plotting
    h  <-  0.05
    # Nei-style model
    r.0.Dat   <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 0,  x = x, h = h, generations = 20000)
    r.1.Dat   <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 1,  x = x, h = h, generations = 20000)
    r.2.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 2, x = x, h = h, generations = 20000)
    r.3.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 3, x = x, h = h, generations = 20000)

    r.0.rel.w.YI.t   <-  r.0.Dat$rel.w.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.1.rel.w.YI.t   <-  r.1.Dat$rel.w.YI.t[round(r.1.Dat$YI.t,digits=4) > 0]
    r.2.rel.w.YI.t  <-  r.2.Dat$rel.w.YI.t[round(r.2.Dat$YI.t,digits=4) > 0]
    r.3.rel.w.YI.t  <-  r.3.Dat$rel.w.YI.t[round(r.3.Dat$YI.t,digits=4) > 0]

    # Haploid model
    r.0.Hap.Dat  <-  makeDeterministicSimHaploidData(U = U, s = s, r = 0, x = x, h = h, generations = 20000)
    r.1.Hap.Dat  <-  makeDeterministicSimHaploidData(U = U, s = s, r = 1, x = x, h = h, generations = 20000)
    r.2.Hap.Dat  <-  makeDeterministicSimHaploidData(U = U, s = s, r = 2, x = x, h = h, generations = 20000)
    r.3.Hap.Dat  <-  makeDeterministicSimHaploidData(U = U, s = s, r = 3, x = x, h = h, generations = 20000)

   # Colors
    COLS     <-  colorRampPalette(wes_palette("Zissou1"))(4)

## Panel C
# Effect of initial del. mut. load on inversion relative fitness 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(relFitMin,relFitMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Inversion relative fitnesses
        lines(r.0.rel.w.YI.t ~ seq_along(r.0.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1))
        lines(r.1.rel.w.YI.t ~ seq_along(r.1.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r.2.rel.w.YI.t ~ seq_along(r.2.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.3.rel.w.YI.t ~ seq_along(r.3.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1))
        # Inversion relative fitness -- Haploid model
        if(wHap) {
            lines(rel.w.YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.0.Hap.Dat)
            lines(rel.w.YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.1.Hap.Dat)
            lines(rel.w.YI.t[1:length(r.2.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.2.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.2.Hap.Dat)
            lines(rel.w.YI.t[1:length(r.3.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.3.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.3.Hap.Dat)
        }
        # Benchmark Inversion relative fitness of 1
        abline(h=1, lwd=1, lty=2, col=1)
        # Points depicting (arbitrary) extinction frequency of 0.00001
        #points(r.0.rel.w.YI.t[length(r.0.rel.w.YI.t)] ~ length(r.0.rel.w.YI.t), pch=8, col='#252525')
        points(r.1.rel.w.YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r.2.rel.w.YI.t[length(r.2.rel.w.YI.t)] ~ length(r.2.rel.w.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
        points(r.3.rel.w.YI.t[length(r.3.rel.w.YI.t)] ~ length(r.3.rel.w.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        # Annotate 1 + Ux approximation
        points((1 + U*x) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, expression(paste(italic(Ux/(hs)==8))), cex=2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(d)==0),
                            expression(italic(d)==1),
                            expression(italic(d)==2),
                            expression(italic(d)==3)),
               lty     =  1,
               lwd     =  2,
               col     =  c(transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[4], opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

qstr       <-  qstar(U=U, q.0=r.0.Dat$YI.t[1], x=x, s=s, h=h)
t          <-  seq_along(r.0.rel.w.YI.t)
YI.gensol  <-  YI.genSol(U=U, YI.init=r.0.Dat$YI.t[1], x=x, s=s, h=h, t=t)

## Panel F
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,(qstr*1.25)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Inversion Frequencies
        lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(YI.t[1:length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.2.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.2.Dat)
        lines(YI.t[1:length(r.3.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.3.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.3.Dat)
        lines(YI.gensol ~ t, lwd=2, lty=3)
        if(wHap) {
            lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.0.Hap.Dat)
            lines(YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.1.Hap.Dat)
            lines(YI.t[1:length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.2.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.2.Hap.Dat)
            lines(YI.t[1:length(r.3.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.3.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.3.Hap.Dat)
        }
        # Points indicating inversion extinctions
        #points(YI.t[length(r.0.rel.w.YI.t)] ~ length(r.0.rel.w.YI.t), pch=8, col='#252525', data=r.0.Dat)
        points(YI.t[length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.1.rel.w.YI.t)], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(YI.t[length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.2.rel.w.YI.t)], pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.2.Dat)
        points(YI.t[length(r.3.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.3.rel.w.YI.t)], pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.3.Dat)
       # benchmark qstar
        abline(h=qstr, lty=2)
        proportionalLabel(0.9, 0.83, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'F', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # reset Colors
     COLS  <-  c('red', '#252525')

## Panel I

    # reset Colors
     COLS  <-  c('red', '#252525')

## Panel I
# deleterious allele frequencies for r = 1
     # Make plot
    delMax  <-  r.1.Dat$qt.Xf.W[1]*1.05
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines for deleterious allele frequencies
        lines(qt.Xf.W[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.Xf.W[1:length(r.1.rel.w.YI.t)]), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qt.Xf.D[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.Xf.D[1:length(r.1.rel.w.YI.t)]), lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qt.YI.W[1:length(r.1.rel.w.YI.t)]/YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.YI.W[1:length(r.1.rel.w.YI.t)]), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
        # Points indicating inversion extinctions
        points(qt.Xf.W[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(qt.Xf.D[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(qt.YI.W[length(r.1.rel.w.YI.t)]/YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.1.Dat)
        # axes
        axis(1, las=1, cex=1.3, labels=NA)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'I', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(d)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(-0.3,  0.5,   expression(paste("Deleterious Allele Frequency")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        
#        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        
     
     # Make plot
    delMax  <-  max(r.1.Dat$qt.Xf.W)
    delMin  <-  min(r.1.Dat$qt.Xf.D)
    delRange  <-  delMax - delMin
    delMax  <-  delMax + 0.05*delRange

     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(delMin,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Deleterious allele frequencies
        lines(qt.Xf.W[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.Xf.W)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qt.Xf.D[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.Xf.D)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
#        lines(qt.YI.W[1:length(r.1.rel.w.YI.t)]/YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.YI.W)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
        # Points indicating inversion extinctions
        points(qt.Xf.W[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(qt.Xf.D[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
#        points(qt.YI.W[length(r.1.rel.w.YI.t)]/YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.1.Dat)
        # axes
        axis(1, las=1, cex=1.2)
        axis(2, las=1, cex=1.2)        
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'L', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(d)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        


}

#####################################









# Supp FIGURE . Illustration of deterministic dynamics
#           for different dominance values
deterministicMutRateIllustration  <-  function(wHap=FALSE) {

    # set plot layout
    layout.mat <- matrix(c(1:12), nrow=4, ncol=3, byrow=FALSE)
    layout     <- layout(layout.mat,respect=TRUE)

    # Make data for plotting 
    s        <-  0.01
    h        <-  0.25
    nTot     <-  10^4 

    ###########
    # U = 0.02
    U        <-  0.02
    u        <-  U/nTot
    x        <-  0.

    # make data for plotting
    r.0.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 0, x = x, h = h, generations = 2500)
    r.1.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 1, x = x, h = h, generations = 2500)
    r.0.rel.w.YI.t  <-  r.0.Dat$rel.w.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.1.rel.w.YI.t  <-  r.1.Dat$rel.w.YI.t[round(r.1.Dat$YI.t,digits=4) > 0]

    # Haploid model
    r.0.Hap.Dat  <-  makeDeterministicSimHaploidData(U = U, s = s, r = 0, x = x, h = h, generations = 2500)
    r.1.Hap.Dat  <-  makeDeterministicSimHaploidData(U = U, s = s, r = 1, x = x, h = h, generations = 2500)

   # Colors
    COLS  <- colorRampPalette(wes_palette("Zissou1"))(4)

    # Axes limits
    relFitMin  <-  0.997
    relFitMax  <-  1 + U*x*1.25

## Panel A
# Effect of initial del. mut. load on inversion relative fitness 
    # make plot
    par(omi=c(0.5, 0.5, 0.5, 0.0), mar = c(3,6,1,0), bty='o', xaxt='s', yaxt='s')
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(relFitMin, relFitMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Inversion relative frequencies
        lines(r.0.rel.w.YI.t ~ seq_along(r.0.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1))
        lines(r.1.rel.w.YI.t ~ seq_along(r.1.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1))
        # Inversion relative fitness -- Haploid model
        if(wHap) {
            lines(rel.w.YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("tomato", opacity=1), data=r.0.Hap.Dat)
            lines(rel.w.YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("tomato", opacity=1), data=r.1.Hap.Dat)
        }
        # Benchmark inversion relative fitness of 1
        abline(h=1, lwd=1, lty=2, col=1)
        # Points depicting (arbitrary) extinction frequency of 0.00001
        points(r.1.rel.w.YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        # Annotate 1 + Ux approximation
        points((1 + U*x) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        proportionalArrows(0.125, 0.8475, 0.07, 0.8475, cex=1, adj=c(0.5, 0.5), xpd=NA, log=FALSE, length=0.05)
        proportionalLabel(0.25, 0.85, expression(1+italic(Ux)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, expression(paste(italic(Ux/(hs)==1.6))), cex=2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.45,  0.5,   expression(paste("Inversion Relative Fitness")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(d)==0),
                            expression(italic(d)==1)),
               lty     =  1,
               lwd     =  2,
               col     =  c(transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
    if(wHap){
        legend(
               x       =  usr[2],
               y       =  usr[4]*0.998,
               legend  =  c(expression("Haploid Model")),
               lty     =  3,
               lwd     =  2,
               col     =  'tomato',
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
    }


qstr       <-  qstar(U=U, q.0=r.0.Dat$YI.t[1], x=x, s=s, h=h)
t          <-  seq_along(r.0.rel.w.YI.t)
YI.gensol  <-  YI.genSol(U=U, YI.init=r.0.Dat$YI.t[1], x=x, s=s, h=h, t=t)

## Panel D
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,(qstr*1.25)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Inversion Frequencies
        lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        # General solution for r = 0
        lines(YI.gensol ~ t, lwd=2, lty=3, col="#252525")
        if(wHap) {
            # Inversion frequencies -- Haploid model
            lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("tomato", opacity=1), data=r.0.Hap.Dat)
            lines(YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("tomato", opacity=1), data=r.1.Hap.Dat)
        }
        # Points
        points(YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        # benchmark qstar
        abline(h=qstr, lty=2)
        proportionalLabel(0.9, 0.85, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.45,  0.5,   expression(paste("Inversion Frequency (", italic(Y[I]), ")")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        # Legend
        legend(
               x       =  usr[2]*0.95,
               y       =  usr[4],
               legend  =  c(
                            expression(paste("General solution for ", italic(Y[I]), " when ", italic(d)==0))),
               lty     =  3,
               lwd     =  2,
               col     =  "#252525",
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
        if(wHap) {
        legend(
               x       =  usr[2]*0.6,
               y       =  usr[4]*0.9,
               legend  =  c(
                            expression(paste("Haploid Model for ", italic(Y[I])))),
               lty     =  3,
               lwd     =  2,
               col     =  '#252525',
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
        }

    # reset Colors
     COLS  <-  c('red', '#252525')

## Panel G
# deleterious allele frequencies for r = 1
     # Make plot
    delMax  <-  r.1.Dat$qt.Xf.W[1]*1.05
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(qt.Xf.W[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.Xf.W[1:length(r.1.rel.w.YI.t)]), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qt.Xf.D[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.Xf.D[1:length(r.1.rel.w.YI.t)]), lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qt.YI.W[1:length(r.1.rel.w.YI.t)]/YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.YI.W[1:length(r.1.rel.w.YI.t)]), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
        # Points indicating inversion extinctions
        points(qt.Xf.W[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(qt.Xf.D[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(qt.YI.W[length(r.1.rel.w.YI.t)]/YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.1.Dat)
        # axes
        axis(1, las=1, cex=1.3, labels=NA)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'G', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(d)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.45,  0.5,   expression(paste("Deleterious Allele Frequency")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        
#        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend = c(expression(italic(q)[italic(X[f])]^italic(W)), 
                          expression(italic(q)[italic(Y[I])]^italic(W)), 
                          expression(italic(q)[italic(X[f])]^italic(D))),
               lty     =  c(1,1,3),
               lwd     =  2,
               col     =  c(transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


## Panel J
    delMax  <-  max(r.1.Dat$qt.Xf.W)
    delMin  <-  min(r.1.Dat$qt.Xf.D)
    delRange  <-  delMax - delMin
    delMax  <-  delMin + 1.05*delRange
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(delMin,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(qt.Xf.W[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.Xf.W[1:length(r.1.rel.w.YI.t)]), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qt.Xf.D[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.Xf.D[1:length(r.1.rel.w.YI.t)]), lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
#        lines(qt.YI.W[1:length(r.1.rel.w.YI.t)]/YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.YI.W[1:length(r.1.rel.w.YI.t)]), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
        # Points indicating inversion extinctions
        points(qt.Xf.W[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(qt.Xf.D[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
#        points(qt.YI.W[length(r.1.rel.w.YI.t)]/YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.1.Dat)
        # axes
        axis(1, las=1, cex=1.3)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'J', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(d)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.45,  0.5,   expression(paste("Deleterious Allele Frequency")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend = c(expression(italic(q)[italic(X[f])]^italic(W)), 
#                          expression(italic(q)[italic(Y[I])]^italic(W)), 
                          expression(italic(q)[italic(X[f])]^italic(D))),
               lty     =  c(1,3),
               lwd     =  2,
               col     =  c(transparentColor(COLS[2], opacity=1),
#                            transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )



    ###########
    # U = 0.05
    U  <-  0.05

    # Nei-style SLR model
    r.0.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 0, x = x, h = h, generations = 2500)
    r.1.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 1, x = x, h = h, generations = 2500)
    r.2.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 2, x = x, h = h, generations = 2500)
    r.3.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 3, x = x, h = h, generations = 2500)

    r.0.rel.w.YI.t  <-  r.0.Dat$rel.w.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.1.rel.w.YI.t  <-  r.1.Dat$rel.w.YI.t[round(r.1.Dat$YI.t,digits=4) > 0]
    r.2.rel.w.YI.t  <-  r.2.Dat$rel.w.YI.t[round(r.2.Dat$YI.t,digits=4) > 0]
    r.3.rel.w.YI.t  <-  r.3.Dat$rel.w.YI.t[round(r.3.Dat$YI.t,digits=4) > 0]

    # Haploid Model
    r.0.Hap.Dat  <-  makeDeterministicSimHaploidData(U = U, s = s, r = 0, x = x, h = h, generations = 2500)
    r.1.Hap.Dat  <-  makeDeterministicSimHaploidData(U = U, s = s, r = 1, x = x, h = h, generations = 2500)
    r.2.Hap.Dat  <-  makeDeterministicSimHaploidData(U = U, s = s, r = 2, x = x, h = h, generations = 2500)
    r.3.Hap.Dat  <-  makeDeterministicSimHaploidData(U = U, s = s, r = 3, x = x, h = h, generations = 2500)

   # Colors
    COLS     <-  colorRampPalette(wes_palette("Zissou1"))(4)

    relFitMin  <-  0.9925
    relFitMax  <-  1 + U*x*1.25

## Panel B
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(relFitMin,relFitMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(r.3.rel.w.YI.t ~ seq_along(r.3.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1))
        lines(r.2.rel.w.YI.t ~ seq_along(r.2.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.1.rel.w.YI.t ~ seq_along(r.1.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r.0.rel.w.YI.t ~ seq_along(r.0.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1))
        # Inversion relative fitness -- Haploid model
        if(wHap) {
            lines(rel.w.YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("tomato", opacity=1), data=r.0.Hap.Dat)
            lines(rel.w.YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("tomato", opacity=1), data=r.1.Hap.Dat)
            lines(rel.w.YI.t[1:length(r.2.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.2.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("tomato", opacity=1), data=r.2.Hap.Dat)
            lines(rel.w.YI.t[1:length(r.3.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.3.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("tomato", opacity=1), data=r.3.Hap.Dat)
        }
        # Benchmark inversion relative fitness of 1
        abline(h=1, lwd=1, lty=2, col=1)
        # Points depicting (arbitrary) extinction frequency of 0.00001
        points(r.1.rel.w.YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r.2.rel.w.YI.t[length(r.2.rel.w.YI.t)] ~ length(r.2.rel.w.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
        points(r.3.rel.w.YI.t[length(r.3.rel.w.YI.t)] ~ length(r.3.rel.w.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        # Annotate 1 + Ux approximation
        points((1 + U*x) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, expression(paste(italic(Ux/(hs)==4))), cex=2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(d)==0),
                            expression(italic(d)==1),
                            expression(italic(d)==2),
                            expression(italic(d)==3)),
               lty     =  1,
               lwd     =  2,
               col     =  c(transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[4], opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )



qstr       <-  qstar(U=U, q.0=r.0.Dat$YI.t[1], x=x, s=s, h=h)
t          <-  seq_along(r.0.rel.w.YI.t)
YI.gensol  <-  YI.genSol(U=U, YI.init=r.0.Dat$YI.t[1], x=x, s=s, h=h, t=t)

## Panel E
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
#    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,(qstr*1.25)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(YI.t[1:length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.2.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.2.Dat)
        lines(YI.t[1:length(r.3.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.3.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.3.Dat)
        # General solution for r = 0
        lines(YI.gensol ~ t, lwd=2, lty=3, col="#252525")
        if(wHap) {
            # Inversion frequencies -- Haploid model
            lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.0.Hap.Dat)
            lines(YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.1.Hap.Dat)
            lines(YI.t[1:length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.2.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.2.Hap.Dat)
            lines(YI.t[1:length(r.3.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.3.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.3.Hap.Dat)
        }
        # Points
        points(YI.t[length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.1.rel.w.YI.t)], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(YI.t[length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.2.rel.w.YI.t)], pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.2.Dat)
        points(YI.t[length(r.3.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.3.rel.w.YI.t)], pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.3.Dat)
        # benchmark qstar
        abline(h=qstr, lty=2)
        proportionalLabel(0.9, 0.85, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'E', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # reset Colors
     COLS  <-  c('red', '#252525')

## Panel H
# deleterious allele frequencies for r = 5
     # Make plot
     delMax  <-  r.1.Dat$qt.Xf.W[1]*1.05
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(qt.Xf.W[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.Xf.W[1:length(r.1.rel.w.YI.t)]), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qt.Xf.D[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.Xf.D[1:length(r.1.rel.w.YI.t)]), lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qt.YI.W[1:length(r.1.rel.w.YI.t)]/YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.YI.W[1:length(r.1.rel.w.YI.t)]), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
        # Points indicating inversion extinctions
        points(qt.Xf.W[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(qt.Xf.D[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(qt.YI.W[length(r.1.rel.w.YI.t)]/YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.1.Dat)
        # axes
        axis(1, las=1, cex=1.3, labels=NA)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'H', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(d)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        


## Panel K
    delMax  <-  max(r.1.Dat$qt.Xf.W)
    delMin  <-  min(r.1.Dat$qt.Xf.D)
    delRange  <-  delMax - delMin
    delMax  <-  delMin + 1.05*delRange
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(delMin,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(qt.Xf.W[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.Xf.W[1:length(r.1.rel.w.YI.t)]), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qt.Xf.D[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.Xf.D[1:length(r.1.rel.w.YI.t)]), lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
#        lines(qt.YI.W[1:length(r.1.rel.w.YI.t)]/YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.YI.W[1:length(r.1.rel.w.YI.t)]), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
        # Points indicating inversion extinctions
        points(qt.Xf.W[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(qt.Xf.D[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
#        points(qt.YI.W[length(r.1.rel.w.YI.t)]/YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.1.Dat)
        # axes
        axis(1, las=1, cex=1.3)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'K', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(d)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(-0.45,  0.5,   expression(paste("Deleterious Allele Frequency")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        




    ###########
    U  <-  0.1

    # First SLR model
    r.0.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 0, x = x, h = h, generations = 2500)
    r.2.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 2, x = x, h = h, generations = 2500)
    r.4.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 4, x = x, h = h, generations = 2500)
    r.6.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 6, x = x, h = h, generations = 2500)

    r.0.rel.w.YI.t  <-  r.0.Dat$rel.w.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.2.rel.w.YI.t  <-  r.2.Dat$rel.w.YI.t[round(r.2.Dat$YI.t,digits=4) > 0]
    r.4.rel.w.YI.t  <-  r.4.Dat$rel.w.YI.t[round(r.4.Dat$YI.t,digits=4) > 0]
    r.6.rel.w.YI.t  <-  r.6.Dat$rel.w.YI.t[round(r.6.Dat$YI.t,digits=4) > 0]

    # Haploid model
    r.0.Hap.Dat  <-  makeDeterministicSimHaploidData(U = U, s = s, r = 0,  x = x, h = h, generations = 2500)
    r.2.Hap.Dat  <-  makeDeterministicSimHaploidData(U = U, s = s, r = 2,  x = x, h = h, generations = 2500)
    r.4.Hap.Dat  <-  makeDeterministicSimHaploidData(U = U, s = s, r = 4, x = x, h = h, generations = 2500)
    r.6.Hap.Dat  <-  makeDeterministicSimHaploidData(U = U, s = s, r = 6, x = x, h = h, generations = 2500)

   # Colors
    COLS     <-  colorRampPalette(wes_palette("Zissou1"))(4)

    relFitMin  <-  0.985
    relFitMax  <-  1 + U*x*1.25

## Panel C
# Effect of initial del. mut. load on inversion relative fitness 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(relFitMin,relFitMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Inversion relative fitnesses
        lines(r.0.rel.w.YI.t ~ seq_along(r.0.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1))
        lines(r.2.rel.w.YI.t ~ seq_along(r.2.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r.4.rel.w.YI.t ~ seq_along(r.4.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.6.rel.w.YI.t ~ seq_along(r.6.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1))
        # Inversion relative fitness -- Haploid model
        if(wHap) {
            lines(rel.w.YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("tomato", opacity=1), data=r.0.Hap.Dat)
            lines(rel.w.YI.t[1:length(r.2.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.2.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("tomato", opacity=1), data=r.2.Hap.Dat)
            lines(rel.w.YI.t[1:length(r.4.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.4.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("tomato", opacity=1), data=r.4.Hap.Dat)
            lines(rel.w.YI.t[1:length(r.6.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.6.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("tomato", opacity=1), data=r.6.Hap.Dat)
        }
        # Benchmakr relative fitness of 1
        abline(h=1, lwd=1, lty=2, col=1)
        # Points depicting (arbitrary) extinction frequency of 0.00001
#        points(r.0.rel.w.YI.t[length(r.0.rel.w.YI.t)] ~ length(r.0.rel.w.YI.t), pch=8, col='#252525')
        points(r.2.rel.w.YI.t[length(r.2.rel.w.YI.t)] ~ length(r.2.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r.4.rel.w.YI.t[length(r.4.rel.w.YI.t)] ~ length(r.4.rel.w.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
        points(r.6.rel.w.YI.t[length(r.6.rel.w.YI.t)] ~ length(r.6.rel.w.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        # Annotate 1 + Ux approximation
        points((1+U*x) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, expression(paste(italic(Ux/(hs)==8))), cex=2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(d)==0),
                            expression(italic(d)==2),
                            expression(italic(d)==4),
                            expression(italic(d)==6)),
               lty     =  1,
               lwd     =  2,
               col     =  c(transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[4], opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


qstr  <-  qstar(U=U, q.0=r.0.Dat$YI.t[1], x=x, s=s, h=h)
t          <-  seq_along(r.0.rel.w.YI.t)
YI.gensol  <-  YI.genSol(U=U, YI.init=r.0.Dat$YI.t[1], x=x, s=s, h=h, t=t)

## Panel F
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,(qstr*1.25)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.2.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.2.Dat)
        lines(YI.t[1:length(r.4.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.4.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.4.Dat)
        lines(YI.t[1:length(r.6.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.6.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.6.Dat)
            # General solution for r = 0
            lines(YI.gensol ~ t, lwd=2, lty=3, col="#252525")
        if(wHap) {
            # Inversion frequencies -- Haploid model
            lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("tomato", opacity=1), data=r.0.Hap.Dat)
            lines(YI.t[1:length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.2.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("tomato", opacity=1), data=r.2.Hap.Dat)
            lines(YI.t[1:length(r.4.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.4.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("tomato", opacity=1), data=r.4.Hap.Dat)
            lines(YI.t[1:length(r.6.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.6.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("tomato", opacity=1), data=r.6.Hap.Dat)
        }
        # Points
#        points(YI.t[length(r.0.rel.w.YI.t)] ~ length(r.0.rel.w.YI.t), pch=8, col='#252525', data=r.0.Dat)
        points(YI.t[length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.2.rel.w.YI.t)], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.2.Dat)
        points(YI.t[length(r.4.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.4.rel.w.YI.t)], pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.4.Dat)
        points(YI.t[length(r.6.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.6.rel.w.YI.t)], pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.6.Dat)
        # benchmark qstar
        abline(h=qstr, lty=2)
        proportionalLabel(0.9, 0.85, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'F', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # reset Colors
     COLS  <-  c('red', '#252525')

## Panel I
# deleterious allele frequencies for r = 5
     # Make plot
    delMax  <-  r.2.Dat$qt.Xf.W[1]*1.05
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(qt.Xf.W[1:length(r.2.rel.w.YI.t)] ~ seq_along(qt.Xf.W[1:length(r.2.rel.w.YI.t)]), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.2.Dat)
        lines(qt.Xf.D[1:length(r.2.rel.w.YI.t)] ~ seq_along(qt.Xf.D[1:length(r.2.rel.w.YI.t)]), lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.2.Dat)
        lines(qt.YI.W[1:length(r.2.rel.w.YI.t)]/YI.t[1:length(r.2.rel.w.YI.t)] ~ seq_along(qt.YI.W[1:length(r.2.rel.w.YI.t)]), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.2.Dat)
        # Points indicating inversion extinctions
        points(qt.Xf.W[length(r.2.rel.w.YI.t)] ~ length(r.2.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.2.Dat)
        points(qt.Xf.D[length(r.2.rel.w.YI.t)] ~ length(r.2.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.2.Dat)
        points(qt.YI.W[length(r.2.rel.w.YI.t)]/YI.t[length(r.2.rel.w.YI.t)] ~ length(r.2.rel.w.YI.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.2.Dat)
        # axes
        axis(1, las=1, cex=1.2, labels=NA)
        axis(2, las=1, cex=1.2)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'I', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(d)==2), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        

## Panel L
    delMax  <-  max(r.2.Dat$qt.Xf.W)
    delMin  <-  min(r.2.Dat$qt.Xf.D)
    delRange  <-  delMax - delMin
    delMax  <-  delMin + 1.05*delRange
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(delMin,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(qt.Xf.W[1:length(r.2.rel.w.YI.t)] ~ seq_along(qt.Xf.W[1:length(r.2.rel.w.YI.t)]), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.2.Dat)
        lines(qt.Xf.D[1:length(r.2.rel.w.YI.t)] ~ seq_along(qt.Xf.D[1:length(r.2.rel.w.YI.t)]), lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.2.Dat)
#        lines(qt.YI.W[1:length(r.2.rel.w.YI.t)]/YI.t[1:length(r.2.rel.w.YI.t)] ~ seq_along(qt.YI.W[1:length(r.2.rel.w.YI.t)]), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.2.Dat)
        # Points indicating inversion extinctions
        points(qt.Xf.W[length(r.2.rel.w.YI.t)] ~ length(r.2.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.2.Dat)
        points(qt.Xf.D[length(r.2.rel.w.YI.t)] ~ length(r.2.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.2.Dat)
#        points(qt.YI.W[length(r.2.rel.w.YI.t)]/YI.t[length(r.2.rel.w.YI.t)] ~ length(r.2.rel.w.YI.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.2.Dat)
        # axes
        axis(1, las=1, cex=1.3)
        axis(2, las=1, cex=1.3)
#        axis(2, las=1, cex=1.3, at = axTicks(2), label = sciNotation(axTicks(2), 5))
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'L', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(d)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(-0.45,  0.5,   expression(paste("Deleterious Allele Frequency")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        

}


sciNotation <- function(x, digits = 1) {
    if (length(x) > 1) {
        return(append(sciNotation(x[1]), sciNotation(x[-1])))
    }
    if (!x) return(0)
    exponent <- floor(log10(x))
    base <- round(x / 10^exponent, digits)
    as.expression(substitute(base %*% 10^exponent, 
            list(base = base, exponent = exponent)))
}


#######################################
#######################################

#' Supplementary Figure -- Autosomal Inversions. Exact Recursions
#' 
#' Reproducing Nei et al. (1967) Fig. 1.
Nei_1967_recessiveMutations_AutoExact  <-  function() {

    # Constant Parameters
    u  <-  10^-5

    # make data for plotting
    case1  <-  makeAutoDeterministicFigSimData_Nei(nTot = 100, U = u*100, I.0 = 10^-5, 
                                                   s = 0.01, r = 0, x = 1, h = 0, generations = 6000)
    case2  <-  makeAutoDeterministicFigSimData_Nei(nTot = 1000, U = u*1000, I.0 = 10^-5, 
                                                   s = 0.01, r = 0, x = 1, h = 0, generations = 6000)
    case3  <-  makeAutoDeterministicFigSimData_Nei(nTot = 100, U = u*100, I.0 = 10^-5, 
                                                   s = 0.1, r = 0, x = 1, h = 0, generations = 6000)
    case4  <-  makeAutoDeterministicFigSimData_Nei(nTot = 1000, U = u*1000, I.0 = 10^-5, 
                                                   s = 0.1, r = 0, x = 1, h = 0, generations = 6000)
    case5  <-  makeAutoDeterministicFigSimData_Nei(nTot = 100, U = u*100, I.0 = 10^-3, 
                                                   s = 0.01, r = 0, x = 1, h = 0, generations = 6000)
    case6  <-  makeAutoDeterministicFigSimData_Nei(nTot = 1000, U = u*1000, I.0 = 10^-3, 
                                                   s = 0.01, r = 0, x = 1, h = 0, generations = 6000)
    case7  <-  makeAutoDeterministicFigSimData_Nei(nTot = 100, U = u*100, I.0 = 10^-3, 
                                                   s = 0.1, r = 0, x = 1, h = 0, generations = 6000)
    case8  <-  makeAutoDeterministicFigSimData_Nei(nTot = 1000, U = u*1000, I.0 = 10^-3, 
                                                   s = 0.1, r = 0, x = 1, h = 0, generations = 6000)

    ## Make Plot
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(case1)), ylim = c(10^-5,1.1), log='y', ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80', log='y')
        box()
        # Inversion Frequencies
        lines(I.t ~ seq_along(I.t), lwd=2, lty=1, col=transparentColor('#252525', opacity=1), data=case1)
        lines(I.t ~ seq_along(I.t), lwd=2, lty=1, col=transparentColor('#252525', opacity=1), data=case2)
        lines(I.t ~ seq_along(I.t), lwd=2, lty=1, col=transparentColor('#252525', opacity=1), data=case3)
        lines(I.t ~ seq_along(I.t), lwd=2, lty=1, col=transparentColor('#252525', opacity=1), data=case4)
        lines(I.t ~ seq_along(I.t), lwd=2, lty=1, col=transparentColor('#252525', opacity=1), data=case5)
        lines(I.t ~ seq_along(I.t), lwd=2, lty=1, col=transparentColor('#252525', opacity=1), data=case6)
        lines(I.t ~ seq_along(I.t), lwd=2, lty=1, col=transparentColor('#252525', opacity=1), data=case7)
        lines(I.t ~ seq_along(I.t), lwd=2, lty=1, col=transparentColor('#252525', opacity=1), data=case8)
        # Annotations
        proportionalLabel(0.75, 0.138, expression(paste("3")), cex=1, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.75, 0.265, expression(paste("1")), cex=1, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.75, 0.51, expression(paste("7")), cex=1, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.75, 0.625, expression(paste("5")), cex=1, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.75, 0.84, expression(paste("4")), cex=1, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.25, 0.8, expression(paste("2")), cex=1, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.2, 0.85, expression(paste("8")), cex=1, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.14, 0.9, expression(paste("6")), cex=1, adj=c(0.5, 0.5), xpd=NA, log='y')
        # axes
        axis(1, las=1, cex=1.3)
        axis(2, las=1, at=c(10^-5, 
                            10^-4, 
                            10^-3, 
                            10^-2, 
                            10^-1,
                            1.0), 
                    labels=c(expression(10^-5), 
                             expression(10^-4), 
                             expression(10^-3), 
                             expression(10^-2), 
                             expression(10^-1),
                             expression("1.0")))
        # Plot labels etc.
        proportionalLabel( 0.5,  1.05,  expression(paste("Reproducing Figure 1 of Nei et al. (1967)")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(-0.15,  0.5,  expression(paste("Inversion Frequency (", italic(I[t]), ")")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90, log='y')
        proportionalLabel( 0.5,  -0.15, expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, log='y')

}






#' Supplementary Figure -- Autosomal Inversions. Exact Recursions
#' 
#' Illustration of deterministic dynamics for different
#' mutation rate values
deterministicMutRateIllustration_AutoExact  <-  function() {

    # Plot Layout
    layout.mat <- matrix(c(1:9), nrow=3, ncol=3, byrow=FALSE)
    layout     <- layout(layout.mat,respect=TRUE)

    # Constant Parameters
    s        <-  0.01
    h        <-  0.25
    nTot     <-  10^4 

    ###########
    # U = 0.02
    U        <-  0.02
    u        <-  U/nTot
    x        <-  0.2

    # make data for plotting
    r.0.Dat  <-  makeAutoDeterministicFigSimData_Nei(U = U, s = s, r = 0, x = x, h = h, generations = 2500)
    r.1.Dat  <-  makeAutoDeterministicFigSimData_Nei(U = U, s = s, r = 1, x = x, h = h, generations = 2500)
    r.0.rel.w.I.t  <-  c()
    r.1.rel.w.I.t  <-  c()
    for(i in 2:length(r.0.Dat$I.t)) {
        r.0.rel.w.I.t[i-1]  <-  r.0.Dat$I.t[i]/r.0.Dat$I.t[i-1]
        r.1.rel.w.I.t[i-1]  <-  r.1.Dat$I.t[i]/r.1.Dat$I.t[i-1]
    }
    r.0.rel.w.I.t  <-  c(r.0.rel.w.I.t,NA)
    r.1.rel.w.I.t  <-  c(r.1.rel.w.I.t,NA)
    r.1.rel.w.I.t  <-  r.1.rel.w.I.t[r.1.Dat$I.t > 5*10^-5]

   # Colors
    COLS  <- colorRampPalette(wes_palette("Zissou1"))(4)

    # Axes limits
    relFitMin  <-  0.997
    relFitMax  <-  1 + U*x*1.25


## Panel A
# Effect of initial del. mut. load on inversion relative fitness 
    # make plot
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(relFitMin, relFitMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Inversion relative frequencies
        lines(r.0.rel.w.I.t ~ seq_along(r.0.rel.w.I.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1))
        lines(r.1.rel.w.I.t ~ seq_along(r.1.rel.w.I.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1))
        # Benchmark inversion relative fitness of 1
        abline(h=1, lwd=1, lty=2, col=1)
        # Points depicting (arbitrary) extinction frequency of 0.00001
        points(r.1.rel.w.I.t[length(r.1.rel.w.I.t)] ~ length(r.1.rel.w.I.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        # Annotate 1 + Ux approximation
        points((1 + U*x) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        proportionalArrows(0.125, 0.86, 0.06, 0.86, cex=1, adj=c(0.5, 0.5), xpd=NA, log=FALSE, length=0.05)
        proportionalLabel(0.25, 0.86, expression(1+italic(Ux)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, expression(paste(italic(Ux/(hs)==1.6))), cex=2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,   expression(paste("Inversion Relative Fitness")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(d)==0),
                            expression(italic(d)==1)),
               lty     =  1,
               lwd     =  2,
               col     =  c(transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

qstr       <-  qstar(U=U, q.0=r.0.Dat$I.t[1], x=x, s=s, h=h)
t          <-  seq_along(r.0.rel.w.I.t)
I.gensol   <-  YI.genSol(U=U, YI.init=r.0.Dat$I.t[1], x=x, s=s, h=h, t=t)

## Panel D
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,(qstr*1.25)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Inversion Frequencies
        lines(I.t[1:length(r.0.rel.w.I.t)] ~ seq_along(I.t)[1:length(r.0.rel.w.I.t)], lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        lines(I.t[1:length(r.1.rel.w.I.t)] ~ seq_along(I.t)[1:length(r.1.rel.w.I.t)], lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        # General solution for r = 0
        lines(I.gensol ~ t, lwd=2, lty=3, col="#252525")
        # Points
        points(I.t[length(r.1.rel.w.I.t)] ~ length(r.1.rel.w.I.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        # benchmark qstar
        abline(h=qstr, lty=2)
        proportionalLabel(0.9, 0.85, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,   expression(paste("Inversion Frequency (", italic(I[t]), ")")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        # Legend
        legend(
               x       =  usr[2]*0.95,
               y       =  usr[4],
               legend  =  c(
                            expression(paste("General solution for ", italic(I[t]), " when ", italic(d)==0))),
               lty     =  3,
               lwd     =  2,
               col     =  "#252525",
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

    # reset Colors
     COLS  <-  c('red', '#252525')

## Panel G
# deleterious allele frequencies for r = 1
     # Make plot
    delMax  <-  r.1.Dat$qt.W[1]*1.05
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines((qt.W/(1-I.t))[1:length(r.1.rel.w.I.t)] ~ seq_along(qt.W[1:length(r.1.rel.w.I.t)]), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines((qt.D/(1-I.t))[1:length(r.1.rel.w.I.t)] ~ seq_along(qt.D[1:length(r.1.rel.w.I.t)]), lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines((qt.I.W/I.t)[1:length(r.1.rel.w.I.t)] ~ seq_along(qt.I.W[1:length(r.1.rel.w.I.t)]), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
        # Points indicating inversion extinctions
        points((qt.W/(1-I.t))[length(r.1.rel.w.I.t)] ~ length(r.1.rel.w.I.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points((qt.D/(1-I.t))[length(r.1.rel.w.I.t)] ~ length(r.1.rel.w.I.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points((qt.I.W/I.t)[length(r.1.rel.w.I.t)] ~ length(r.1.rel.w.I.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.1.Dat)
        # axes
        axis(1, las=1, cex=1.3)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'G', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(d)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,   expression(paste("Deleterious Allele Frequency")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend = c(expression(italic(Q)[italic(N)]^italic(W)), 
                          expression(italic(Q)[italic(I)]^italic(W)), 
                          expression(italic(Q)[italic(N)]^italic(D))),
               lty     =  c(1,1,3),
               lwd     =  2,
               col     =  c(transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )



    ###########
    # U = 0.05
    U  <-  0.05
    r.0.Dat  <-  makeAutoDeterministicFigSimData_Nei(U = U, s = s, r = 0, x = x, h = h, generations = 3000)
    r.1.Dat  <-  makeAutoDeterministicFigSimData_Nei(U = U, s = s, r = 1, x = x, h = h, generations = 3000)
    r.2.Dat  <-  makeAutoDeterministicFigSimData_Nei(U = U, s = s, r = 2, x = x, h = h, generations = 3000)
    r.3.Dat  <-  makeAutoDeterministicFigSimData_Nei(U = U, s = s, r = 3, x = x, h = h, generations = 3000)

    r.0.rel.w.I.t  <-  c()
    r.1.rel.w.I.t  <-  c()
    r.2.rel.w.I.t  <-  c()
    r.3.rel.w.I.t  <-  c()
    for(i in 2:length(r.0.Dat$I.t)) {
        r.0.rel.w.I.t[i-1]  <-  r.0.Dat$I.t[i]/r.0.Dat$I.t[i-1]
        r.1.rel.w.I.t[i-1]  <-  r.1.Dat$I.t[i]/r.1.Dat$I.t[i-1]
        r.2.rel.w.I.t[i-1]  <-  r.2.Dat$I.t[i]/r.2.Dat$I.t[i-1]
        r.3.rel.w.I.t[i-1]  <-  r.3.Dat$I.t[i]/r.3.Dat$I.t[i-1]
    }
    r.0.rel.w.I.t  <-  c(r.0.rel.w.I.t,NA)
    r.1.rel.w.I.t  <-  c(r.1.rel.w.I.t,NA)
    r.2.rel.w.I.t  <-  c(r.2.rel.w.I.t,NA)
    r.3.rel.w.I.t  <-  c(r.3.rel.w.I.t,NA)
    r.1.rel.w.I.t  <-  r.1.rel.w.I.t[r.1.Dat$I.t > 5*10^-5]
    r.2.rel.w.I.t  <-  r.2.rel.w.I.t[r.2.Dat$I.t > 5*10^-5]
    r.3.rel.w.I.t  <-  r.3.rel.w.I.t[r.3.Dat$I.t > 5*10^-5]

   # Colors
    COLS     <-  colorRampPalette(wes_palette("Zissou1"))(4)

    relFitMin  <-  0.9925
    relFitMax  <-  1 + U*x*1.25

## Panel B
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(relFitMin,relFitMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(r.3.rel.w.I.t ~ seq_along(r.3.rel.w.I.t), lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1))
        lines(r.2.rel.w.I.t ~ seq_along(r.2.rel.w.I.t), lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.1.rel.w.I.t ~ seq_along(r.1.rel.w.I.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r.0.rel.w.I.t ~ seq_along(r.0.rel.w.I.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1))
        # Benchmark inversion relative fitness of 1
        abline(h=1, lwd=1, lty=2, col=1)
        # Points depicting (arbitrary) extinction frequency of 0.00001
        points(r.1.rel.w.I.t[length(r.1.rel.w.I.t)] ~ length(r.1.rel.w.I.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r.2.rel.w.I.t[length(r.2.rel.w.I.t)] ~ length(r.2.rel.w.I.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
        points(r.3.rel.w.I.t[length(r.3.rel.w.I.t)] ~ length(r.3.rel.w.I.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        # Annotate 1 + Ux approximation
        points((1 + U*x) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, expression(paste(italic(Ux/(hs)==4))), cex=2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(d)==0),
                            expression(italic(d)==1),
                            expression(italic(d)==2),
                            expression(italic(d)==3)),
               lty     =  1,
               lwd     =  2,
               col     =  c(transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[4], opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )



qstr       <-  qstar(U=U, q.0=r.0.Dat$I.t[1], x=x, s=s, h=h)
t          <-  seq_along(r.0.rel.w.I.t)
I.gensol  <-  YI.genSol(U=U, YI.init=r.0.Dat$I.t[1], x=x, s=s, h=h, t=t)

## Panel E
# Effect of initial del. mut. load on inversion frequency dynamics 
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,(qstr*1.25)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(I.t[1:length(r.0.rel.w.I.t)] ~ seq_along(I.t)[1:length(r.0.rel.w.I.t)], lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        lines(I.t[1:length(r.1.rel.w.I.t)] ~ seq_along(I.t)[1:length(r.1.rel.w.I.t)], lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(I.t[1:length(r.2.rel.w.I.t)] ~ seq_along(I.t)[1:length(r.2.rel.w.I.t)], lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.2.Dat)
        lines(I.t[1:length(r.3.rel.w.I.t)] ~ seq_along(I.t)[1:length(r.3.rel.w.I.t)], lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.3.Dat)
        # General solution for r = 0
        lines(I.gensol ~ t, lwd=2, lty=3, col="#252525")
        # Points
        points(I.t[length(r.1.rel.w.I.t)] ~ seq_along(I.t)[length(r.1.rel.w.I.t)], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(I.t[length(r.2.rel.w.I.t)] ~ seq_along(I.t)[length(r.2.rel.w.I.t)], pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.2.Dat)
        points(I.t[length(r.3.rel.w.I.t)] ~ seq_along(I.t)[length(r.3.rel.w.I.t)], pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.3.Dat)
        # benchmark qstar
        abline(h=qstr, lty=2)
        proportionalLabel(0.9, 0.85, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'E', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # reset Colors
     COLS  <-  c('red', '#252525')

## Panel H
# deleterious allele frequencies for r = 1
     # Make plot
     delMax  <-  r.1.Dat$qt.W[1]*1.05
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines((qt.W/(1-I.t))[1:length(r.1.rel.w.I.t)] ~ seq_along(qt.W[1:length(r.1.rel.w.I.t)]), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines((qt.D/(1-I.t))[1:length(r.1.rel.w.I.t)] ~ seq_along(qt.D[1:length(r.1.rel.w.I.t)]), lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines((qt.I.W/I.t)[1:length(r.1.rel.w.I.t)] ~ seq_along(qt.I.W[1:length(r.1.rel.w.I.t)]), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
        # Points indicating inversion extinctions
        points((qt.W/(1-I.t))[length(r.1.rel.w.I.t)] ~ length(r.1.rel.w.I.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points((qt.D/(1-I.t))[length(r.1.rel.w.I.t)] ~ length(r.1.rel.w.I.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points((qt.I.W/I.t)[length(r.1.rel.w.I.t)] ~ length(r.1.rel.w.I.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.1.Dat)
        # axes
        axis(1, las=1, cex=1.3)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'H', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(d)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        




    ###########
    # U = 0.1
    U  <-  0.1
    r.0.Dat   <-  makeAutoDeterministicFigSimData_Nei(U = U, s = s, r = 0,  x = x, h = h, generations = 3000)
    r.2.Dat   <-  makeAutoDeterministicFigSimData_Nei(U = U, s = s, r = 2,  x = x, h = h, generations = 3000)
    r.4.Dat  <-  makeAutoDeterministicFigSimData_Nei(U = U, s = s, r = 4, x = x, h = h, generations = 3000)
    r.6.Dat  <-  makeAutoDeterministicFigSimData_Nei(U = U, s = s, r = 6, x = x, h = h, generations = 3000)

    r.0.rel.w.I.t  <-  c()
    r.2.rel.w.I.t  <-  c()
    r.4.rel.w.I.t  <-  c()
    r.6.rel.w.I.t  <-  c()
    for(i in 2:length(r.0.Dat$I.t)) {
        r.0.rel.w.I.t[i-1]  <-  r.0.Dat$I.t[i]/r.0.Dat$I.t[i-1]
        r.2.rel.w.I.t[i-1]  <-  r.2.Dat$I.t[i]/r.2.Dat$I.t[i-1]
        r.4.rel.w.I.t[i-1]  <-  r.4.Dat$I.t[i]/r.4.Dat$I.t[i-1]
        r.6.rel.w.I.t[i-1]  <-  r.6.Dat$I.t[i]/r.6.Dat$I.t[i-1]
    }
    r.0.rel.w.I.t  <-  c(r.0.rel.w.I.t,NA)
    r.2.rel.w.I.t  <-  c(r.2.rel.w.I.t,NA)
    r.4.rel.w.I.t  <-  c(r.4.rel.w.I.t,NA)
    r.6.rel.w.I.t  <-  c(r.6.rel.w.I.t,NA)
    r.2.rel.w.I.t  <-  r.2.rel.w.I.t[r.2.Dat$I.t > 5*10^-5]
    r.4.rel.w.I.t  <-  r.4.rel.w.I.t[r.4.Dat$I.t > 5*10^-5]
    r.6.rel.w.I.t  <-  r.6.rel.w.I.t[r.6.Dat$I.t > 5*10^-5]


   # Colors
    COLS     <-  colorRampPalette(wes_palette("Zissou1"))(4)

    relFitMin  <-  0.985
    relFitMax  <-  1 + U*x*1.25

## Panel C
# Effect of initial del. mut. load on inversion relative fitness 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(relFitMin,relFitMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Inversion relative fitnesses
        lines(r.0.rel.w.I.t ~ seq_along(r.0.rel.w.I.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1))
        lines(r.2.rel.w.I.t ~ seq_along(r.2.rel.w.I.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r.4.rel.w.I.t ~ seq_along(r.4.rel.w.I.t), lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.6.rel.w.I.t ~ seq_along(r.6.rel.w.I.t), lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1))
        # Benchmakr relative fitness of 1
        abline(h=1, lwd=1, lty=2, col=1)
        # Points depicting (arbitrary) extinction frequency of 0.00001
        points(r.0.rel.w.I.t[length(r.0.rel.w.I.t)] ~ length(r.0.rel.w.I.t), pch=8, col='#252525')
        points(r.2.rel.w.I.t[length(r.2.rel.w.I.t)] ~ length(r.2.rel.w.I.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r.4.rel.w.I.t[length(r.4.rel.w.I.t)] ~ length(r.4.rel.w.I.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
        points(r.6.rel.w.I.t[length(r.6.rel.w.I.t)] ~ length(r.6.rel.w.I.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        # Annotate 1 + Ux approximation
        points((1+U*x) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, expression(paste(italic(Ux/(hs)==8))), cex=2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(d)==0),
                            expression(italic(d)==2),
                            expression(italic(d)==4),
                            expression(italic(d)==6)),
               lty     =  1,
               lwd     =  2,
               col     =  c(transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[4], opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


qstr      <-  qstar(U=U, q.0=r.0.Dat$I.t[1], x=x, s=s, h=h)
t         <-  seq_along(r.0.rel.w.I.t)
I.gensol  <-  YI.genSol(U=U, YI.init=r.0.Dat$I.t[1], x=x, s=s, h=h, t=t)

## Panel F
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,(qstr*1.25)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        lines(I.t[1:length(r.0.rel.w.I.t)] ~ seq_along(I.t)[1:length(r.0.rel.w.I.t)], lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        lines(I.t[1:length(r.2.rel.w.I.t)] ~ seq_along(I.t)[1:length(r.2.rel.w.I.t)], lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.2.Dat)
        lines(I.t[1:length(r.4.rel.w.I.t)] ~ seq_along(I.t)[1:length(r.4.rel.w.I.t)], lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.4.Dat)
        lines(I.t[1:length(r.6.rel.w.I.t)] ~ seq_along(I.t)[1:length(r.6.rel.w.I.t)], lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.6.Dat)
        # General solution for r = 0
        lines(I.gensol ~ t, lwd=2, lty=3, col="#252525")
        # Points
        points(I.t[length(r.2.rel.w.I.t)] ~ seq_along(I.t)[length(r.2.rel.w.I.t)], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.2.Dat)
        points(I.t[length(r.4.rel.w.I.t)] ~ seq_along(I.t)[length(r.4.rel.w.I.t)], pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.4.Dat)
        points(I.t[length(r.6.rel.w.I.t)] ~ seq_along(I.t)[length(r.6.rel.w.I.t)], pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.6.Dat)
        # benchmark qstar
        abline(h=qstr, lty=2)
        proportionalLabel(0.9, 0.85, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'F', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # reset Colors
     COLS  <-  c('red', '#252525')

## Panel I
# deleterious allele frequencies for r = 5
     # Make plot
    delMax  <-  r.2.Dat$qt.W[1]*1.05
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines((qt.W/(1-I.t))[1:length(r.2.rel.w.I.t)] ~ seq_along(qt.W[1:length(r.2.rel.w.I.t)]), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.2.Dat)
        lines((qt.D/(1-I.t))[1:length(r.2.rel.w.I.t)] ~ seq_along(qt.D[1:length(r.2.rel.w.I.t)]), lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.2.Dat)
        lines((qt.I.W/I.t)[1:length(r.2.rel.w.I.t)] ~ seq_along(qt.I.W[1:length(r.2.rel.w.I.t)]), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.2.Dat)
        # Points indicating inversion extinctions
        points((qt.W/(1-I.t))[length(r.2.rel.w.I.t)] ~ length(r.2.rel.w.I.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.2.Dat)
        points((qt.D/(1-I.t))[length(r.2.rel.w.I.t)] ~ length(r.2.rel.w.I.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.2.Dat)
        points((qt.I.W/I.t)[length(r.2.rel.w.I.t)] ~ length(r.2.rel.w.I.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.2.Dat)
        # axes
        axis(1, las=1, cex=1.2)
        axis(2, las=1, cex=1.2)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'I', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(d)==2), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        

}





# Supp FIGURE . Illustration of deterministic dynamics
#           for different mutation rate values
deterministicDominanceIllustration_AutoExact  <-  function(x=0.2) {

    # set plot layout
    layout.mat <- matrix(c(1:9), nrow=3, ncol=3, byrow=FALSE)
    layout     <- layout(layout.mat,respect=TRUE)

    # Make data for plotting 
    s        <-  0.01
    h        <-  0.25
    nTot     <-  10^4 

    ###########
    # U = 0.02
    U        <-  0.02
    u        <-  U/nTot
    

    # make data for plotting
    r.0.Dat  <-  makeAutoDeterministicFigSimData_Nei(U = U, s = s, r = 0, x = x, h = h, generations = 2500)
    r.1.Dat  <-  makeAutoDeterministicFigSimData_Nei(U = U, s = s, r = 1, x = x, h = h, generations = 2500)
    r.0.rel.w.I.t  <-  c()
    r.1.rel.w.I.t  <-  c()
    for(i in 2:length(r.0.Dat$I.t)) {
        r.0.rel.w.I.t[i-1]  <-  r.0.Dat$I.t[i]/r.0.Dat$I.t[i-1]
        r.1.rel.w.I.t[i-1]  <-  r.1.Dat$I.t[i]/r.1.Dat$I.t[i-1]
    }
    r.0.rel.w.I.t  <-  c(r.0.rel.w.I.t,NA)
    r.1.rel.w.I.t  <-  c(r.1.rel.w.I.t,NA)
    r.1.rel.w.I.t  <-  r.1.rel.w.I.t[r.1.Dat$I.t > 5*10^-5]

   # Colors
    COLS  <- colorRampPalette(wes_palette("Zissou1"))(4)

    # Axes limits
    relFitMin  <-  0.996
    relFitMax  <-  1 + U*x*1.25


## Panel A
# Effect of initial del. mut. load on inversion relative fitness 
    # make plot
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(relFitMin, relFitMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Inversion relative frequencies
        lines(r.0.rel.w.I.t ~ seq_along(r.0.rel.w.I.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1))
        lines(r.1.rel.w.I.t ~ seq_along(r.1.rel.w.I.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1))
        # Benchmark inversion relative fitness of 1
        abline(h=1, lwd=1, lty=2, col=1)
        # Points depicting (arbitrary) extinction frequency of 0.00001
        points(r.1.rel.w.I.t[length(r.1.rel.w.I.t)] ~ length(r.1.rel.w.I.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        # Annotate 1 + Ux approximation
        points((1 + U*x) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        proportionalArrows(0.125, 0.86, 0.06, 0.86, cex=1, adj=c(0.5, 0.5), xpd=NA, log=FALSE, length=0.05)
        proportionalLabel(0.25, 0.86, expression(1+italic(Ux)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, expression(paste(italic(Ux/(hs)==1.6))), cex=2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,   expression(paste("Inversion Relative Fitness")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(d)==0),
                            expression(italic(d)==1)),
               lty     =  1,
               lwd     =  2,
               col     =  c(transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

qstr       <-  qstar(U=U, q.0=r.0.Dat$I.t[1], x=x, s=s, h=h)
t          <-  seq_along(r.0.rel.w.I.t)
I.gensol   <-  YI.genSol(U=U, YI.init=r.0.Dat$I.t[1], x=x, s=s, h=h, t=t)

## Panel D
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,(qstr*1.25)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Inversion Frequencies
        lines(I.t[1:length(r.0.rel.w.I.t)] ~ seq_along(I.t)[1:length(r.0.rel.w.I.t)], lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        lines(I.t[1:length(r.1.rel.w.I.t)] ~ seq_along(I.t)[1:length(r.1.rel.w.I.t)], lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        # General solution for r = 0
        lines(I.gensol ~ t, lwd=2, lty=3, col="#252525")
        # Points
        points(I.t[length(r.1.rel.w.I.t)] ~ length(r.1.rel.w.I.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        # benchmark qstar
        abline(h=qstr, lty=2)
        proportionalLabel(0.9, 0.83, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,   expression(paste("Inversion Frequency (", italic(I[t]), ")")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        # Legend
        legend(
               x       =  usr[2]*0.95,
               y       =  usr[4],
               legend  =  c(
                            expression(paste("General solution for ", italic(I[t]), " when ", italic(d)==0))),
               lty     =  3,
               lwd     =  2,
               col     =  "#252525",
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

    # reset Colors
     COLS  <-  c('red', '#252525')

## Panel G
# deleterious allele frequencies for r = 1
     # Make plot
    delMax  <-  r.1.Dat$qt.W[1]*1.05
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines((qt.W/(1-I.t))[1:length(r.1.rel.w.I.t)] ~ seq_along(qt.W[1:length(r.1.rel.w.I.t)]), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines((qt.D/(1-I.t))[1:length(r.1.rel.w.I.t)] ~ seq_along(qt.D[1:length(r.1.rel.w.I.t)]), lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines((qt.I.W/I.t)[1:length(r.1.rel.w.I.t)] ~ seq_along(qt.I.W[1:length(r.1.rel.w.I.t)]), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
        # Points indicating inversion extinctions
        points((qt.W/(1-I.t))[length(r.1.rel.w.I.t)] ~ length(r.1.rel.w.I.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points((qt.D/(1-I.t))[length(r.1.rel.w.I.t)] ~ length(r.1.rel.w.I.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points((qt.I.W/I.t)[length(r.1.rel.w.I.t)] ~ length(r.1.rel.w.I.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.1.Dat)
        # axes
        axis(1, las=1, cex=1.3)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'G', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(d)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,   expression(paste("Deleterious Allele Frequency")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend = c(expression(italic(Q)[italic(N)]^italic(W)), 
                          expression(italic(Q)[italic(I)]^italic(W)), 
                          expression(italic(Q)[italic(N)]^italic(D))),
               lty     =  c(1,1,3),
               lwd     =  2,
               col     =  c(transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )



    ###########
    # h = 0.1
    h  <-  0.1
    r.0.Dat  <-  makeAutoDeterministicFigSimData_Nei(U = U, s = s, r = 0, x = x, h = h, generations = 8000)
    r.1.Dat  <-  makeAutoDeterministicFigSimData_Nei(U = U, s = s, r = 1, x = x, h = h, generations = 8000)
    r.2.Dat  <-  makeAutoDeterministicFigSimData_Nei(U = U, s = s, r = 2, x = x, h = h, generations = 8000)
    r.3.Dat  <-  makeAutoDeterministicFigSimData_Nei(U = U, s = s, r = 3, x = x, h = h, generations = 8000)

    r.0.rel.w.I.t  <-  c()
    r.1.rel.w.I.t  <-  c()
    r.2.rel.w.I.t  <-  c()
    r.3.rel.w.I.t  <-  c()
    for(i in 2:length(r.0.Dat$I.t)) {
        r.0.rel.w.I.t[i-1]  <-  r.0.Dat$I.t[i]/r.0.Dat$I.t[i-1]
        r.1.rel.w.I.t[i-1]  <-  r.1.Dat$I.t[i]/r.1.Dat$I.t[i-1]
        r.2.rel.w.I.t[i-1]  <-  r.2.Dat$I.t[i]/r.2.Dat$I.t[i-1]
        r.3.rel.w.I.t[i-1]  <-  r.3.Dat$I.t[i]/r.3.Dat$I.t[i-1]
    }
    r.0.rel.w.I.t  <-  c(r.0.rel.w.I.t,NA)
    r.1.rel.w.I.t  <-  c(r.1.rel.w.I.t,NA)
    r.2.rel.w.I.t  <-  c(r.2.rel.w.I.t,NA)
    r.3.rel.w.I.t  <-  c(r.3.rel.w.I.t,NA)
    r.1.rel.w.I.t  <-  r.1.rel.w.I.t[r.1.Dat$I.t > 5*10^-5]
    r.2.rel.w.I.t  <-  r.2.rel.w.I.t[r.2.Dat$I.t > 5*10^-5]
    r.3.rel.w.I.t  <-  r.3.rel.w.I.t[r.3.Dat$I.t > 5*10^-5]

   # Colors
    COLS     <-  colorRampPalette(wes_palette("Zissou1"))(4)

## Panel B
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(relFitMin,relFitMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(r.3.rel.w.I.t ~ seq_along(r.3.rel.w.I.t), lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1))
        lines(r.2.rel.w.I.t ~ seq_along(r.2.rel.w.I.t), lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.1.rel.w.I.t ~ seq_along(r.1.rel.w.I.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r.0.rel.w.I.t ~ seq_along(r.0.rel.w.I.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1))
        # Benchmark inversion relative fitness of 1
        abline(h=1, lwd=1, lty=2, col=1)
        # Points depicting (arbitrary) extinction frequency of 0.00001
        points(r.1.rel.w.I.t[length(r.1.rel.w.I.t)] ~ length(r.1.rel.w.I.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r.2.rel.w.I.t[length(r.2.rel.w.I.t)] ~ length(r.2.rel.w.I.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
        points(r.3.rel.w.I.t[length(r.3.rel.w.I.t)] ~ length(r.3.rel.w.I.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        # Annotate 1 + Ux approximation
        points((1 + U*x) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, expression(paste(italic(Ux/(hs)==4))), cex=2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(d)==0),
                            expression(italic(d)==1),
                            expression(italic(d)==2),
                            expression(italic(d)==3)),
               lty     =  1,
               lwd     =  2,
               col     =  c(transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[4], opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )



qstr       <-  qstar(U=U, q.0=r.0.Dat$I.t[1], x=x, s=s, h=h)
t          <-  seq_along(r.0.rel.w.I.t)
I.gensol   <-  YI.genSol(U=U, YI.init=r.0.Dat$I.t[1], x=x, s=s, h=h, t=t)

## Panel E
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,(qstr*1.25)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(I.t[1:length(r.0.rel.w.I.t)] ~ seq_along(I.t)[1:length(r.0.rel.w.I.t)], lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        lines(I.t[1:length(r.1.rel.w.I.t)] ~ seq_along(I.t)[1:length(r.1.rel.w.I.t)], lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(I.t[1:length(r.2.rel.w.I.t)] ~ seq_along(I.t)[1:length(r.2.rel.w.I.t)], lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.2.Dat)
        lines(I.t[1:length(r.3.rel.w.I.t)] ~ seq_along(I.t)[1:length(r.3.rel.w.I.t)], lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.3.Dat)
        # General solution for r = 0
        lines(I.gensol ~ t, lwd=2, lty=3, col="#252525")
        # Points
        points(I.t[length(r.1.rel.w.I.t)] ~ seq_along(I.t)[length(r.1.rel.w.I.t)], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(I.t[length(r.2.rel.w.I.t)] ~ seq_along(I.t)[length(r.2.rel.w.I.t)], pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.2.Dat)
        points(I.t[length(r.3.rel.w.I.t)] ~ seq_along(I.t)[length(r.3.rel.w.I.t)], pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.3.Dat)
        # benchmark qstar
        abline(h=qstr, lty=2)
        proportionalLabel(0.9, 0.83, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'E', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # reset Colors
     COLS  <-  c('red', '#252525')

## Panel H
# deleterious allele frequencies for r = 1
     # Make plot
     delMax  <-  r.1.Dat$qt.W[1]*1.05
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines((qt.W/(1-I.t))[1:length(r.1.rel.w.I.t)] ~ seq_along(qt.W[1:length(r.1.rel.w.I.t)]), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines((qt.D/(1-I.t))[1:length(r.1.rel.w.I.t)] ~ seq_along(qt.D[1:length(r.1.rel.w.I.t)]), lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines((qt.I.W/I.t)[1:length(r.1.rel.w.I.t)] ~ seq_along(qt.I.W[1:length(r.1.rel.w.I.t)]), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
        # Points indicating inversion extinctions
        points((qt.W/(1-I.t))[length(r.1.rel.w.I.t)] ~ length(r.1.rel.w.I.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points((qt.D/(1-I.t))[length(r.1.rel.w.I.t)] ~ length(r.1.rel.w.I.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points((qt.I.W/I.t)[length(r.1.rel.w.I.t)] ~ length(r.1.rel.w.I.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.1.Dat)
        # axes
        axis(1, las=1, cex=1.3)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'H', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(d)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        




    ###########
    # h = 0.05
    h  <-  0.05
    r.0.Dat   <-  makeAutoDeterministicFigSimData_Nei(U = U, s = s, r = 0,  x = x, h = h, generations = 10000)
    r.5.Dat   <-  makeAutoDeterministicFigSimData_Nei(U = U, s = s, r = 2,  x = x, h = h, generations = 10000)
    r.10.Dat  <-  makeAutoDeterministicFigSimData_Nei(U = U, s = s, r = 4, x = x, h = h, generations = 10000)
    r.15.Dat  <-  makeAutoDeterministicFigSimData_Nei(U = U, s = s, r = 6, x = x, h = h, generations = 10000)

    r.0.rel.w.I.t  <-  c()
    r.5.rel.w.I.t  <-  c()
    r.10.rel.w.I.t  <-  c()
    r.15.rel.w.I.t  <-  c()
    for(i in 2:length(r.0.Dat$I.t)) {
        r.0.rel.w.I.t[i-1]  <-  r.0.Dat$I.t[i]/r.0.Dat$I.t[i-1]
        r.5.rel.w.I.t[i-1]  <-  r.5.Dat$I.t[i]/r.5.Dat$I.t[i-1]
        r.10.rel.w.I.t[i-1]  <-  r.10.Dat$I.t[i]/r.10.Dat$I.t[i-1]
        r.15.rel.w.I.t[i-1]  <-  r.15.Dat$I.t[i]/r.15.Dat$I.t[i-1]
    }
    r.0.rel.w.I.t  <-  r.0.rel.w.I.t[r.0.Dat$I.t < 0.9999]
    r.5.rel.w.I.t  <-  c(r.5.rel.w.I.t,NA)
    r.10.rel.w.I.t  <-  c(r.10.rel.w.I.t,NA)
    r.15.rel.w.I.t  <-  c(r.15.rel.w.I.t,NA)
    r.5.rel.w.I.t  <-  r.5.rel.w.I.t[r.5.Dat$I.t > 5*10^-5]
    r.10.rel.w.I.t  <-  r.10.rel.w.I.t[r.10.Dat$I.t > 5*10^-5]
    r.15.rel.w.I.t  <-  r.15.rel.w.I.t[r.15.Dat$I.t > 5*10^-5]


   # Colors
    COLS     <-  colorRampPalette(wes_palette("Zissou1"))(4)

## Panel C
# Effect of initial del. mut. load on inversion relative fitness 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(relFitMin,relFitMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Inversion relative fitnesses
        lines(r.0.rel.w.I.t[1:length(r.0.rel.w.I.t)] ~ seq_along(r.0.rel.w.I.t[1:length(r.0.rel.w.I.t)]), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1))
        lines(r.5.rel.w.I.t ~ seq_along(r.5.rel.w.I.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r.10.rel.w.I.t ~ seq_along(r.10.rel.w.I.t), lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.15.rel.w.I.t ~ seq_along(r.15.rel.w.I.t), lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1))
        # Benchmakr relative fitness of 1
        abline(h=1, lwd=1, lty=2, col=1)
        # Points depicting (arbitrary) extinction frequency of 0.00001
#        points(r.0.rel.w.I.t[length(r.0.rel.w.I.t)] ~ length(r.0.rel.w.I.t), pch=8, col='#252525')
        points(r.5.rel.w.I.t[length(r.5.rel.w.I.t)] ~ length(r.5.rel.w.I.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r.10.rel.w.I.t[length(r.10.rel.w.I.t)] ~ length(r.10.rel.w.I.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
        points(r.15.rel.w.I.t[length(r.15.rel.w.I.t)] ~ length(r.15.rel.w.I.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        # Annotate 1 + Ux approximation
        points((1+U*x) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, expression(paste(italic(Ux/(hs)==8))), cex=2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(d)==0),
                            expression(italic(d)==2),
                            expression(italic(d)==4),
                            expression(italic(d)==6)),
               lty     =  1,
               lwd     =  2,
               col     =  c(transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[4], opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


qstr      <-  qstar(U=U, q.0=r.0.Dat$I.t[1], x=x, s=s, h=h)
t         <-  seq_along(r.0.rel.w.I.t)
I.gensol  <-  YI.genSol(U=U, YI.init=r.0.Dat$I.t[1], x=x, s=s, h=h, t=t)

## Panel F
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,(qstr*1.25)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        lines(I.t[1:length(r.0.rel.w.I.t)] ~ seq_along(I.t)[1:length(r.0.rel.w.I.t)], lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        lines(I.t[1:length(r.5.rel.w.I.t)] ~ seq_along(I.t)[1:length(r.5.rel.w.I.t)], lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.5.Dat)
        lines(I.t[1:length(r.10.rel.w.I.t)] ~ seq_along(I.t)[1:length(r.10.rel.w.I.t)], lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.10.Dat)
        lines(I.t[1:length(r.15.rel.w.I.t)] ~ seq_along(I.t)[1:length(r.15.rel.w.I.t)], lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.15.Dat)
        # General solution for r = 0
        lines(I.gensol ~ t, lwd=2, lty=3, col="#252525")
        # Points
#        points(I.t[length(r.0.rel.w.I.t)] ~ length(r.0.rel.w.I.t), pch=8, col='#252525', data=r.0.Dat)
        points(I.t[length(r.5.rel.w.I.t)] ~ seq_along(I.t)[length(r.5.rel.w.I.t)], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.5.Dat)
        points(I.t[length(r.10.rel.w.I.t)] ~ seq_along(I.t)[length(r.10.rel.w.I.t)], pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.10.Dat)
        points(I.t[length(r.15.rel.w.I.t)] ~ seq_along(I.t)[length(r.15.rel.w.I.t)], pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.15.Dat)
        # benchmark qstar
        abline(h=qstr, lty=2)
        proportionalLabel(0.9, 0.83, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'F', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # reset Colors
     COLS  <-  c('red', '#252525')

## Panel I
# deleterious allele frequencies for r = 5
     # Make plot
    delMax  <-  r.5.Dat$qt.W[1]*1.05
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines((qt.W/(1-I.t))[1:length(r.5.rel.w.I.t)] ~ seq_along(qt.W[1:length(r.5.rel.w.I.t)]), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.5.Dat)
        lines((qt.D/(1-I.t))[1:length(r.5.rel.w.I.t)] ~ seq_along(qt.D[1:length(r.5.rel.w.I.t)]), lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.5.Dat)
        lines((qt.I.W/I.t)[1:length(r.5.rel.w.I.t)] ~ seq_along(qt.I.W[1:length(r.5.rel.w.I.t)]), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.5.Dat)
        # Points indicating inversion extinctions
        points((qt.W/(1-I.t))[length(r.5.rel.w.I.t)] ~ length(r.5.rel.w.I.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.5.Dat)
        points((qt.D/(1-I.t))[length(r.5.rel.w.I.t)] ~ length(r.5.rel.w.I.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.5.Dat)
        points((qt.I.W/I.t)[length(r.5.rel.w.I.t)] ~ length(r.5.rel.w.I.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.5.Dat)
        # axes
        axis(1, las=1, cex=1.2)
        axis(2, las=1, cex=1.2)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'I', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(d)==5), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        

}





############################################
#' Supplementary figure(s) illustrating 
#' deterministic dynamics for SLR-expanding
#'  inversions under a variety of dominance
#'  and mutation rate values
deterministicSuppFig  <-  function(h = 0.25, U = 0.02) {


    # set plot layout
    layout.mat <- matrix(c(1:20), nrow=4, ncol=5, byrow=FALSE)
    layout     <- layout(layout.mat,respect=TRUE)

    # Make data for plotting (x = 0.1; 1st column of panels)
    s        <-  0.01
    nTot     <-  10^4 
    u        <-  U/nTot
    qHat     <-  (U/(nTot*h*s))
    x        <-  0.05
    rBar     <-  (U*x)/(h*s)
    r.vals   <-  c(0,1)
    # make data for plotting
    r.0.Dat  <-  makeDeterministicFigSimData_SLR_Nei(r = 0, U = U, x = x, h = h, generations = 3500)
    r.1.Dat  <-  makeDeterministicFigSimData_SLR_Nei(r = 1, U = U, x = x, h = h, generations = 3500)

    r.0.rel.w.YI.t   <-  r.0.Dat$rel.w.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.1.rel.w.YI.t   <-  r.1.Dat$rel.w.YI.t[round(r.1.Dat$YI.t,digits=4) > 0]

    r.0.Hap.Dat  <-  makeDeterministicSimHaploidData(r = 0, U = U, x = x, h = h, generations = 3500)
    r.1.Hap.Dat  <-  makeDeterministicSimHaploidData(r = 1, U = U, x = x, h = h, generations = 3500)

    # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  rev(colfunc(5))


## Panel A
# Distribution of captured deleterious alleles (r) given inversion size (x)
# with benchmarks for 5th, 10th, 25th, and 40th quantiles, as well as the mean.
    pDens  <-  dpois(x=c(0:10), lambda=rBar)
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,length(pDens)), ylim = c(0, 0.7), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(pDens~c(0:(length(pDens)-1)), lwd=1.5, lty=1, col=transparentColor("#252525", opacity=1))
        points(pDens~c(0:(length(pDens)-1)), lwd=1.5, lty=1, pch=21, col=transparentColor("#252525", opacity=1), bg=transparentColor("#252525", opacity=0.25))
        abline(v=rBar, col=1, lwd=1)
        abline(v=r.vals, lwd=1, lty=2)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.5, 1.1, substitute(italic(x)==xx, list(xx=x)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.04, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(paste(Pr(italic(d)~"|"~italic(x)))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel( 0.5,  -0.2,  expression(italic(d)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(d)=="{0,1}"),
                            expression(bar(italic(d))==frac(italic(U*x),italic(h*s)))),
               lty     =  c(2,1),
               lwd     =  c(1,1),
               col     =  c(2,1),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

## Panel B
# Effect of initial del. mut. load on inversion relative fitness dynamics 
    # make plot
ymin  <-  min(1-U*x*2, r.1.rel.w.YI.t[length(r.1.rel.w.YI.t)])
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(ymin, 1+U*x*1.7), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(r.0.rel.w.YI.t ~ seq_along(r.0.rel.w.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1))
        lines(r.1.rel.w.YI.t ~ seq_along(r.1.rel.w.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1))
        abline(h=1, lwd=1, lty=2, col=1)
        # Inversion relative fitness -- Haploid model
        lines(rel.w.YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("tomato", opacity=1), data=r.0.Hap.Dat)
        lines(rel.w.YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("tomato", opacity=1), data=r.1.Hap.Dat)
        # Annotate 1 + Ux approximation
        points((1+U*x) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        # Points depicting (arbitrary) extinction frequency of 0.00001
        #points(r.0.rel.w.YI.t[length(r.0.rel.w.YI.t)] ~ length(r.0.rel.w.YI.t), pch=8, col='tomato')
        points(r.1.rel.w.YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        proportionalArrows(0.125, 0.81, 0.06, 0.81, cex=1, adj=c(0.5, 0.5), xpd=NA, log=FALSE, length=0.05)
        proportionalLabel(0.3, 0.81, expression(1+italic(Ux)), cex=1.3, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(paste("Inversion Relative Fitness")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(d)==0),
                            expression(italic(d)==1)),
               lty     =  1,
               lwd     =  1.5,
               col     =  c(transparentColor(COLS[5], opacity=1),
                            transparentColor(COLS[4], opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
        legend(
               x       =  usr[2]*0.7,
               y       =  usr[4],
               legend  =  c(expression("Haploid Model")),
               lty     =  3,
               lwd     =  2,
               col     =  "tomato",
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

## Panel C
# Effect of initial del. mut. load on inversion frequency dynamics 
    qstr  <-  qstar(U=U, q.0=r.0.Dat$YI.t[1], x=x, s=s, h=h)
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,(qstr*1.25)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(YI.t[1:length(r.1.rel.w.YI.t)]), lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.1.Dat)
        # Inversion frequencies -- Haploid model
        lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("tomato", opacity=1), data=r.0.Hap.Dat)
        lines(YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("tomato", opacity=1), data=r.1.Hap.Dat)
        # Points
        points(YI.t[length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.1.rel.w.YI.t)], pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.1.Dat)
        # Benchmark qStar
        abline(h=qstr, lty=2)        
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.9, 0.85, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(paste("Inversion Frequency (", italic(Y[I]), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        legend(
               x       =  usr[2]*0.7,
               y       =  usr[4],
               legend  =  c(expression("Haploid Model")),
               lty     =  3,
               lwd     =  2,
               col     =  "tomato",
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


    # reset Colors
     COLS  <-  c('red', '#252525')

## Panel D
# deleterious allele frequencies for r = 5
     # Make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,max(r.0.Dat$qt.Xf.W)*1.25), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(qt.Xf.W[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.Xf.W[1:length(r.1.rel.w.YI.t)]), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.0.Dat)
        lines(qt.Xf.D[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.Xf.D[1:length(r.1.rel.w.YI.t)]), lwd=1.5, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.0.Dat)
        lines(qt.YI.W[1:length(r.1.rel.w.YI.t)]/YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.YI.W[1:length(r.1.rel.w.YI.t)]), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        # Points
        points(qt.Xf.W[length(r.1.rel.w.YI.t)] ~ seq_along(qt.Xf.W)[length(r.1.rel.w.YI.t)], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.0.Dat)
        points(qt.Xf.D[length(r.1.rel.w.YI.t)] ~ seq_along(qt.Xf.D)[length(r.1.rel.w.YI.t)], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.0.Dat)
        points(qt.YI.W[length(r.1.rel.w.YI.t)]/YI.t[length(r.1.rel.w.YI.t)] ~ seq_along(qt.YI.W)[length(r.1.rel.w.YI.t)], pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.0.Dat)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(italic(d)==0), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(paste("Deleterious Allele Frequency")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4]*0.75,
               legend = c(expression(italic(q)[italic(X[f])]^italic(wt)), 
                          expression(italic(q)[italic(Y[I])]^italic(wt)), 
                          expression(italic(q)[italic(X[f])]^italic(del))),
               lty     =  c(1,1,3),
               lwd     =  1.5,
               col     =  c(transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


## 2nd Column of Panels: x = 0.1
    # Set inversion size, r values
    x        <-  0.1

    # Select r values
    rBar  <-  (U*x)/(h*s)
    r.vals  <-  unique(round(seq(0, rBar, length=3)))

    # make data for plotting
    r.Dat.list      <-  list()
    r.Hap.Dat.list  <-  list()
    r.relW.list     <-  list()
    for(i in 1:length(r.vals)) {
        r.Dat.list[[i]] <-  makeDeterministicFigSimData_SLR_Nei(r = r.vals[i], U = U, x = x, h = h, generations = 3500)
        r.Hap.Dat.list[[i]] <-  makeDeterministicSimHaploidData(r = r.vals[i], U = U, x = x, h = h, generations = 3500)
        if(r.vals[i] == 0) {
            r.relW.list[[i]]  <-  r.Dat.list[[i]]$rel.w.YI.t[round(r.Dat.list[[i]]$YI.t,digits=4) < 1]
        } 
        if(r.vals[i] > 0) {
            r.relW.list[[i]]  <-  r.Dat.list[[i]]$rel.w.YI.t[round(r.Dat.list[[i]]$YI.t,digits=4) > 0]
        }
    }


   # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  rev(colfunc(5))

## Panel E
# Distribution of captured deleterious alleles (r) given inversion size (x)
# with benchmarks for 5th, 10th, 25th, and 40th quantiles, as well as the mean.
    pDens  <-  dpois(x=c(0:10), lambda=rBar)
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,length(pDens)), ylim = c(0, 0.7), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(pDens~c(0:(length(pDens)-1)), lwd=1.5, lty=1, col=transparentColor("#252525", opacity=1))
        points(pDens~c(0:(length(pDens)-1)), lwd=1.5, lty=1, pch=21, col=transparentColor("#252525", opacity=1), bg=transparentColor("#252525", opacity=0.25))
        abline(v=rBar, col=1, lwd=1)
        abline(v=r.vals, col=2, lwd=1, lty=2)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.5, 1.1, substitute(italic(x)==xx, list(xx=x)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.04, 1.05, 'E', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.2,  expression(italic(d)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            substitute(italic(d)=={rs}, list(rs=deparse(r.vals))),
                            expression(bar(italic(d))==frac(italic(U*x),italic(h*s)))),
               lty     =  c(2,1),
               lwd     =  c(1,1),
               col     =  c(2,1),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


## Panel F
# Effect of initial del. mut. load on inversion relative fitness dynamics 
    # make plot
ymin  <-  min(1-U*x, r.relW.list[[i]][length(r.relW.list[[i]])])
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.Dat.list[[1]])), ylim = c(ymin, 1+U*x*1.7), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        for(i in 1:length(r.Dat.list)){
            lines(r.relW.list[[i]] ~ seq_along(r.relW.list[[i]]), lwd=1.5, lty=1, col=transparentColor(COLS[6-i], opacity=1))
            lines(r.Hap.Dat.list[[i]]$rel.w.YI.t[1:length(r.relW.list[[i]])] ~ seq_along(r.relW.list[[i]]), lwd=1.5, lty=3, col='tomato')
        }
        abline(h=1, lwd=1, lty=2, col=1)
        # Annotate 1 + Ux approximation
        points((1+U*x) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        # Points depicting (arbitrary) extinction frequency of 0.00001
        for(i in 2:length(r.relW.list)) {
            points(r.relW.list[[i]][length(r.relW.list[[i]])] ~ length(r.relW.list[[i]]), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[6-i], opacity=0.7))    
        }
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'F', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            paste("r = ", r.vals[1]),
                            paste("r = ", r.vals[2]),
                            paste("r = ", r.vals[3])),
               lty     =  1,
               lwd     =  1.5,
               col     =  c(transparentColor(COLS[5], opacity=1),
                            transparentColor(COLS[4], opacity=1),
                            transparentColor(COLS[3], opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


## Panel G
# Effect of initial del. mut. load on inversion frequency dynamics 
    qstr  <-  qstar(U=U, q.0=r.0.Dat$YI.t[1], x=x, s=s, h=h)
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,(qstr*1.25)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        for(i in 1:length(r.Dat.list)){
            lines(r.Dat.list[[i]]$YI[1:length(r.relW.list[[i]])] ~ seq_along(r.relW.list[[i]]), lwd=1.5, lty=1, col=transparentColor(COLS[6-i], opacity=1))
            lines(r.Hap.Dat.list[[i]]$YI.t[1:length(r.relW.list[[i]])] ~ seq_along(r.relW.list[[i]]), lwd=1.5, lty=3, col='tomato')
        }
        # Points
        for(i in 2:length(r.Dat.list)){
            points(r.Dat.list[[i]]$YI[length(r.relW.list[[i]])] ~ length(r.relW.list[[i]]), pch=21, col=transparentColor(COLS[6-i], opacity=1), bg=transparentColor(COLS[6-i], opacity=0.7))
        }
        # Benchmark qStar
        abline(h=qstr, lty=2)        
        proportionalLabel(0.9, 0.85, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'G', cex=1.2, adj=c(0.5, 0.5), xpd=NA)


    # reset Colors
     COLS  <-  c('red', '#252525')

## Panel H
# deleterious allele frequencies for r = 5
     # Make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.1.Dat)), ylim = c(0,max(r.0.Dat$qt.Xf.W)*1.25), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        lines(qt.Xf.W[1:length(r.relW.list[[2]])] ~ seq_along(qt.Xf.W[1:length(r.relW.list[[2]])]), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.Dat.list[[2]])
        lines(qt.Xf.D[1:length(r.relW.list[[2]])] ~ seq_along(qt.Xf.D[1:length(r.relW.list[[2]])]), lwd=1.5, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.Dat.list[[2]])
        lines(qt.YI.W[1:length(r.relW.list[[2]])]/YI.t[1:length(r.relW.list[[2]])] ~ seq_along(qt.YI.W[1:length(r.relW.list[[2]])]), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.Dat.list[[2]])
        # Points
        points(qt.Xf.W[length(r.relW.list[[2]])] ~ seq_along(qt.Xf.W)[length(r.relW.list[[2]])], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.Dat.list[[2]])
        points(qt.Xf.D[length(r.relW.list[[2]])] ~ seq_along(qt.Xf.D)[length(r.relW.list[[2]])], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.Dat.list[[2]])
        points(qt.YI.W[length(r.relW.list[[2]])]/YI.t[length(r.relW.list[[2]])] ~ seq_along(qt.YI.W)[length(r.relW.list[[2]])], pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.Dat.list[[2]])
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'H', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1, substitute(italic(d)==rs, list(rs=r.vals[2])), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        




## 3rd Column of Panels: x = 0.2
    # Set inversion size, r values
    x        <-  0.2

    # Select r values
    rBar    <-  round((U*x)/(h*s), digits=2)
    r.vals  <-  unique(round(seq(0, rBar, length=4)))

    # make data for plotting
    r.Dat.list      <-  list()
    r.Hap.Dat.list  <-  list()
    r.relW.list     <-  list()
    for(i in 1:length(r.vals)) {
        r.Dat.list[[i]] <-  makeDeterministicFigSimData_SLR_Nei(r = r.vals[i], U = U, x = x, h = h, generations = 3500)
        r.Hap.Dat.list[[i]] <-  makeDeterministicSimHaploidData(r = r.vals[i], U = U, x = x, h = h, generations = 3500)
        if(r.vals[i] == 0) {
            r.relW.list[[i]]  <-  r.Dat.list[[i]]$rel.w.YI.t[round(r.Dat.list[[i]]$YI.t,digits=4) < 1]
        } 
        if(r.vals[i] > 0) {
            r.relW.list[[i]]  <-  r.Dat.list[[i]]$rel.w.YI.t[round(r.Dat.list[[i]]$YI.t,digits=4) > 0]
        }
    }

   # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  rev(colfunc(5))

## Panel I
# Distribution of captured deleterious alleles (r) given inversion size (x)
# with benchmarks for 5th, 10th, 25th, and 40th quantiles, as well as the mean.
    pDens  <-  dpois(x=c(0:max(10,(2.1*rBar))), lambda=rBar)
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,length(pDens)), ylim = c(0, 0.7), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(pDens~c(0:(length(pDens)-1)), lwd=1.5, lty=1, col=transparentColor("#252525", opacity=1))
        points(pDens~c(0:(length(pDens)-1)), lwd=1.5, lty=1, pch=21, col=transparentColor("#252525", opacity=1), bg=transparentColor("#252525", opacity=0.25))
        abline(v=rBar, lwd=1)
        abline(v=r.vals, lwd=1, lty=2)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.5, 1.1, substitute(italic(x)==xx, list(xx=x)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.04, 1.05, 'I', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.2,  expression(italic(d)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            substitute(italic(d)==xx, list(xx=deparse(r.vals))),
                            expression(bar(italic(d))==frac(italic(U*x),italic(h*s)))),
               lty     =  c(2,1),
               lwd     =  c(1,1),
               col     =  c(2,1),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


## Panel J
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.Dat.list[[1]])), ylim = c(1-U*x, 1+U*x*1.7), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        for(i in 1:length(r.Dat.list)){
            lines(r.relW.list[[i]] ~ seq_along(r.relW.list[[i]]), lwd=1.5, lty=1, col=transparentColor(COLS[6-i], opacity=1))
            lines(r.Hap.Dat.list[[i]]$rel.w.YI.t[1:length(r.relW.list[[i]])] ~ seq_along(r.relW.list[[i]]), lwd=1.5, lty=3, col='tomato')
        }
        abline(h=1, lwd=1, lty=2, col=1)
        # Annotate 1 + Ux approximation
        points((1+U*x) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        # Points depicting (arbitrary) extinction frequency of 0.00001
        for(i in 2:length(r.relW.list)) {
            points(r.relW.list[[i]][length(r.relW.list[[i]])] ~ length(r.relW.list[[i]]), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[6-i], opacity=0.7))    
        }
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'J', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            paste("r = ", r.vals[1]),
                            paste("r = ", r.vals[2]),
                            paste("r = ", r.vals[3]),
                            paste("r = ", r.vals[4])),
               lty     =  1,
               lwd     =  1.5,
               col     =  c(transparentColor(COLS[5], opacity=1),
                            transparentColor(COLS[4], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[2], opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )



## Panel K
# Effect of initial del. mut. load on inversion frequency dynamics 
    qstr  <-  qstar(U=U, q.0=r.Dat.list[[1]]$YI.t[1], x=x, s=s, h=h)
    # make plot
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.Dat.list[[1]])), ylim = c(0,(qstr*1.25)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        for(i in 1:length(r.Dat.list)){
            lines(r.Dat.list[[i]]$YI[1:length(r.relW.list[[i]])] ~ seq_along(r.relW.list[[i]]), lwd=1.5, lty=1, col=transparentColor(COLS[6-i], opacity=1))
            lines(r.Hap.Dat.list[[i]]$YI.t[1:length(r.relW.list[[i]])] ~ seq_along(r.relW.list[[i]]), lwd=1.5, lty=3, col='tomato')
        }
        # Points
        for(i in 2:length(r.Dat.list)){
            points(r.Dat.list[[i]]$YI[length(r.relW.list[[i]])] ~ length(r.relW.list[[i]]), pch=21, col=transparentColor(COLS[6-i], opacity=1), bg=transparentColor(COLS[6-i], opacity=0.7))
        }
        # Benchmark qStar
        abline(h=qstr, lty=2)        
        proportionalLabel(0.9, 0.85, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'K', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # reset Colors
     COLS  <-  c('red',
                 '#252525')
## Panel L
# deleterious allele frequencies for r = 1
     # Make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.Dat.list[[1]])), ylim = c(0,max(r.Dat.list[[1]]$qt.Xf.W)*1.25), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        lines(qt.Xf.W[1:length(r.relW.list[[2]])] ~ seq_along(qt.Xf.W[1:length(r.relW.list[[2]])]), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.Dat.list[[2]])
        lines(qt.Xf.D[1:length(r.relW.list[[2]])] ~ seq_along(qt.Xf.D[1:length(r.relW.list[[2]])]), lwd=1.5, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.Dat.list[[2]])
        lines(qt.YI.W[1:length(r.relW.list[[2]])]/YI.t[1:length(r.relW.list[[2]])] ~ seq_along(qt.YI.W[1:length(r.relW.list[[2]])]), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.Dat.list[[2]])
        # Points
        points(qt.Xf.W[length(r.relW.list[[2]])] ~ seq_along(qt.Xf.W)[length(r.relW.list[[2]])], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.Dat.list[[2]])
        points(qt.Xf.D[length(r.relW.list[[2]])] ~ seq_along(qt.Xf.D)[length(r.relW.list[[2]])], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.Dat.list[[2]])
        points(qt.YI.W[length(r.relW.list[[2]])]/YI.t[length(r.relW.list[[2]])] ~ seq_along(qt.YI.W)[length(r.relW.list[[2]])], pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.Dat.list[[2]])
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'L', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1, substitute(italic(d)==rs, list(rs=r.vals[2])), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        



## 4th Column of Panels: x = 0.5
    # Set inversion size, r values
    x        <-  0.5

  # Select r values
    rBar    <-  round((U*x)/(h*s), digits=2)
    r.vals  <-  unique(round(seq(0, rBar, length=5)))

    # make data for plotting
    r.Dat.list      <-  list()
    r.Hap.Dat.list  <-  list()
    r.relW.list     <-  list()
    for(i in 1:length(r.vals)) {
        r.Dat.list[[i]] <-  makeDeterministicFigSimData_SLR_Nei(r = r.vals[i], U = U, x = x, h = h, generations = 3500)
        r.Hap.Dat.list[[i]] <-  makeDeterministicSimHaploidData(r = r.vals[i], U = U, x = x, h = h, generations = 3500)
        if(r.vals[i] == 0) {
            r.relW.list[[i]]  <-  r.Dat.list[[i]]$rel.w.YI.t[round(r.Dat.list[[i]]$YI.t,digits=4) < 1]
        } 
        if(r.vals[i] > 0) {
            r.relW.list[[i]]  <-  r.Dat.list[[i]]$rel.w.YI.t[round(r.Dat.list[[i]]$YI.t,digits=4) > 0]
        }
    }

   # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  rev(colfunc(5))

## Panel M
# Distribution of captured deleterious alleles (r) given inversion size (x)
# with benchmarks for 5th, 10th, 25th, and 40th quantiles, as well as the mean.
    pDens  <-  dpois(x=c(0:max(20, 2.1*rBar)), lambda=rBar)
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,length(pDens)), ylim = c(0, 0.4), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(pDens~c(0:(length(pDens)-1)), lwd=1.5, lty=1, col=transparentColor("#252525", opacity=1))
        points(pDens~c(0:(length(pDens)-1)), lwd=1.5, lty=1, pch=21, col=transparentColor("#252525", opacity=1), bg=transparentColor("#252525", opacity=0.25))
        abline(v=rBar, col=1, lwd=1)
        abline(v=r.vals, col=2, lwd=1, lty=2)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.5, 1.1, substitute(italic(x)==xx, list(xx=x)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.04, 1.05, 'M', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.2,  expression(italic(d)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            substitute(italic(d)==xx, list(xx=deparse(r.vals))),
                            expression(bar(italic(d))==frac(italic(U*x),italic(h*s)))),
               lty     =  c(2,1),
               lwd     =  c(1,1),
               col     =  c(2,1),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


## Panel N
# Effect of initial del. mut. load on inversion relative fitness dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.Dat.list[[1]])), ylim = c(1-U*x, 1+U*x*1.7), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
      plotGrid(lineCol='grey80')
        box()
        # lines
        for(i in 1:length(r.Dat.list)){
            lines(r.relW.list[[i]] ~ seq_along(r.relW.list[[i]]), lwd=1.5, lty=1, col=transparentColor(COLS[6-i], opacity=1))
            lines(r.Hap.Dat.list[[i]]$rel.w.YI.t[1:length(r.relW.list[[i]])] ~ seq_along(r.relW.list[[i]]), lwd=1.5, lty=3, col='tomato')
        }
        abline(h=1, lwd=1, lty=2, col=1)
        # Annotate 1 + Ux approximation
        points((1+U*x) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        # Points depicting (arbitrary) extinction frequency of 0.00001
        for(i in 2:length(r.relW.list)) {
            points(r.relW.list[[i]][length(r.relW.list[[i]])] ~ length(r.relW.list[[i]]), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[6-i], opacity=0.7))    
        }
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'N', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            paste("r = ", r.vals[1]),
                            paste("r = ", r.vals[2]),
                            paste("r = ", r.vals[3]),
                            paste("r = ", r.vals[4])),
               lty     =  1,
               lwd     =  1.5,
               col     =  c(transparentColor(COLS[5], opacity=1),
                            transparentColor(COLS[4], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[2], opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


## Panel O
# Effect of initial del. mut. load on inversion frequency dynamics 
    qstr  <-  qstar(U=U, q.0=r.Dat.list[[1]]$YI.t[1], x=x, s=s, h=h)
    # make plot
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.Dat.list[[1]])), ylim = c(0,(qstr*1.25)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        for(i in 1:length(r.Dat.list)){
            lines(r.Dat.list[[i]]$YI[1:length(r.relW.list[[i]])] ~ seq_along(r.relW.list[[i]]), lwd=1.5, lty=1, col=transparentColor(COLS[6-i], opacity=1))
            lines(r.Hap.Dat.list[[i]]$YI.t[1:length(r.relW.list[[i]])] ~ seq_along(r.relW.list[[i]]), lwd=1.5, lty=3, col='tomato')
        }
        # Points
        for(i in 2:length(r.Dat.list)){
            points(r.Dat.list[[i]]$YI[length(r.relW.list[[i]])] ~ length(r.relW.list[[i]]), pch=21, col=transparentColor(COLS[6-i], opacity=1), bg=transparentColor(COLS[6-i], opacity=0.7))
        }
        # Benchmark qStar
        abline(h=qstr, lty=2)        
        proportionalLabel(0.9, 0.85, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'O', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # reset Colors
     COLS  <-  c('red',
                 '#252525')


## Panel P
# deleterious allele frequencies for r = 1
     # Make plot
plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.Dat.list[[1]])), ylim = c(0,max(r.Dat.list[[1]]$qt.Xf.W)*1.25), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        lines(qt.Xf.W[1:length(r.relW.list[[2]])] ~ seq_along(qt.Xf.W[1:length(r.relW.list[[2]])]), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.Dat.list[[2]])
        lines(qt.Xf.D[1:length(r.relW.list[[2]])] ~ seq_along(qt.Xf.D[1:length(r.relW.list[[2]])]), lwd=1.5, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.Dat.list[[2]])
        lines(qt.YI.W[1:length(r.relW.list[[2]])]/YI.t[1:length(r.relW.list[[2]])] ~ seq_along(qt.YI.W[1:length(r.relW.list[[2]])]), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.Dat.list[[2]])
        # Points
        points(qt.Xf.W[length(r.relW.list[[2]])] ~ seq_along(qt.Xf.W)[length(r.relW.list[[2]])], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.Dat.list[[2]])
        points(qt.Xf.D[length(r.relW.list[[2]])] ~ seq_along(qt.Xf.D)[length(r.relW.list[[2]])], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.Dat.list[[2]])
        points(qt.YI.W[length(r.relW.list[[2]])]/YI.t[length(r.relW.list[[2]])] ~ seq_along(qt.YI.W)[length(r.relW.list[[2]])], pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.Dat.list[[2]])
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'P', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1, substitute(italic(d)==rs, list(rs=r.vals[2])), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        




## 5th Column of Panels: x = 0.8
    # Set inversion size, r values
    x    <-  0.8

  # Select r values
    rBar    <-  round((U*x)/(h*s), digits=2)
    r.vals  <-  unique(round(seq(0, rBar, length=5)))

    # make data for plotting
    r.Dat.list      <-  list()
    r.Hap.Dat.list  <-  list()
    r.relW.list     <-  list()
    for(i in 1:length(r.vals)) {
        r.Dat.list[[i]] <-  makeDeterministicFigSimData_SLR_Nei(r = r.vals[i], U = U, x = x, h = h, generations = 3500)
        r.Hap.Dat.list[[i]] <-  makeDeterministicSimHaploidData(r = r.vals[i], U = U, x = x, h = h, generations = 3500)
        if(r.vals[i] == 0) {
            r.relW.list[[i]]  <-  r.Dat.list[[i]]$rel.w.YI.t[round(r.Dat.list[[i]]$YI.t,digits=4) < 1]
        } 
        if(r.vals[i] > 0) {
            r.relW.list[[i]]  <-  r.Dat.list[[i]]$rel.w.YI.t[round(r.Dat.list[[i]]$YI.t,digits=4) > 0]
        }
    }
   
   # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  rev(colfunc(5))

## Panel Q
# Distribution of captured deleterious alleles (r) given inversion size (x)
# with benchmarks for 5th, 10th, 25th, and 40th quantiles, as well as the mean.
    qHat  <-  r.0.Dat$qt.Y.W[1]
    rBar  <-  round(nTot*x*qHat, digits=2)
    pDens  <-  dpois(x=c(0:max(20, 2.1*rBar)), lambda=U*x/(h*s))
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,length(pDens)), ylim = c(0, 0.7), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(pDens~c(0:(length(pDens)-1)), lwd=1.5, lty=1, col=transparentColor("#252525", opacity=1))
        points(pDens~c(0:(length(pDens)-1)), lwd=1.5, lty=1, pch=21, col=transparentColor("#252525", opacity=1), bg=transparentColor("#252525", opacity=0.25))
        abline(v=rBar, col=1, lwd=1)
        abline(v=r.vals, col=2, lwd=1, lty=2)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.5, 1.1, substitute(italic(x)==xx, list(xx=x)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.04, 1.05, 'Q', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.2,  expression(italic(d)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            substitute(italic(d)==xx, list(xx=deparse(r.vals))),
                            expression(bar(italic(d))==frac(italic(U*x),italic(h*s)))),
               lty     =  c(2,1),
               lwd     =  c(1,1),
               col     =  c(2,1),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


## Panel R
# Effect of initial del. mut. load on inversion relative fitness
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.Dat.list[[1]])), ylim = c(1-U*x, 1+U*x*1.7), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        for(i in 1:length(r.Dat.list)){
            lines(r.relW.list[[i]] ~ seq_along(r.relW.list[[i]]), lwd=1.5, lty=1, col=transparentColor(COLS[6-i], opacity=1))
            lines(r.Hap.Dat.list[[i]]$rel.w.YI.t[1:length(r.relW.list[[i]])] ~ seq_along(r.relW.list[[i]]), lwd=1.5, lty=3, col='tomato')
        }
        abline(h=1, lwd=1, lty=2, col=1)
        # Annotate 1 + Ux approximation
        points((1+U*x) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        # Points depicting (arbitrary) extinction frequency of 0.00001
        for(i in 2:length(r.relW.list)) {
            points(r.relW.list[[i]][length(r.relW.list[[i]])] ~ length(r.relW.list[[i]]), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[6-i], opacity=0.7))    
        }
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'R', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            paste("r = ", r.vals[1]),
                            paste("r = ", r.vals[2]),
                            paste("r = ", r.vals[3]),
                            paste("r = ", r.vals[4])),
               lty     =  1,
               lwd     =  1.5,
               col     =  c(transparentColor(COLS[5], opacity=1),
                            transparentColor(COLS[4], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[2], opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )




## Panel S
# Effect of initial del. mut. load on inversion frequency dynamics 
    qstr  <-  qstar(U=U, q.0=r.Dat.list[[1]]$YI.t[1], x=x, s=s, h=h)
    # make plot
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.Dat.list[[1]])), ylim = c(0,(qstr*1.25)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        for(i in 1:length(r.Dat.list)){
            lines(r.Dat.list[[i]]$YI[1:length(r.relW.list[[i]])] ~ seq_along(r.relW.list[[i]]), lwd=1.5, lty=1, col=transparentColor(COLS[6-i], opacity=1))
            lines(r.Hap.Dat.list[[i]]$YI.t[1:length(r.relW.list[[i]])] ~ seq_along(r.relW.list[[i]]), lwd=1.5, lty=3, col='tomato')
        }
        # Points
        for(i in 2:length(r.Dat.list)){
            points(r.Dat.list[[i]]$YI[length(r.relW.list[[i]])] ~ length(r.relW.list[[i]]), pch=21, col=transparentColor(COLS[6-i], opacity=1), bg=transparentColor(COLS[6-i], opacity=0.7))
        }
        # Benchmark qStar
        abline(h=qstr, lty=2)        
        proportionalLabel(0.9, 0.85, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'S', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # reset Colors
     COLS  <-  c('red',
                 '#252525')


## Panel T
# deleterious allele frequencies for r = 1
     # Make plot
        plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.Dat.list[[1]])), ylim = c(0,max(r.Dat.list[[1]]$qt.Xf.W)*1.25), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        lines(qt.Xf.W[1:length(r.relW.list[[2]])] ~ seq_along(qt.Xf.W[1:length(r.relW.list[[2]])]), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.Dat.list[[2]])
        lines(qt.Xf.D[1:length(r.relW.list[[2]])] ~ seq_along(qt.Xf.D[1:length(r.relW.list[[2]])]), lwd=1.5, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.Dat.list[[2]])
        lines(qt.YI.W[1:length(r.relW.list[[2]])]/YI.t[1:length(r.relW.list[[2]])] ~ seq_along(qt.YI.W[1:length(r.relW.list[[2]])]), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.Dat.list[[2]])
        # Points
        points(qt.Xf.W[length(r.relW.list[[2]])] ~ seq_along(qt.Xf.W)[length(r.relW.list[[2]])], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.Dat.list[[2]])
        points(qt.Xf.D[length(r.relW.list[[2]])] ~ seq_along(qt.Xf.D)[length(r.relW.list[[2]])], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.Dat.list[[2]])
        points(qt.YI.W[length(r.relW.list[[2]])]/YI.t[length(r.relW.list[[2]])] ~ seq_along(qt.YI.W)[length(r.relW.list[[2]])], pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.Dat.list[[2]])
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'T', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1, substitute(italic(d)==rs, list(rs=r.vals[2])), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        

}





#######################################
#######################################
#######################################



#' Fig. 2
#' Fixation probability of different size inversions 
#' expanding the SDR on a Y chromosome: Pr(fix | x) 
#' Illustration of Pr(fix | x) for different 
#' population size and U/s ratios.
PrFixFig4Panel_Approx  <-  function() {

    # Import data for plotting

    ## Autosomal Inversion Data
    ## N = 100k data
    AdatN100kU02  <-  read.csv(file = './data/RECODE/PrFixAutoFig_h0.25_s0.01_N100k_U02.csv', header=TRUE)
    AdatN100kU05  <-  read.csv(file = './data/RECODE/PrFixAutoFig_h0.25_s0.01_N100k_U05.csv', header=TRUE)
    AdatN100kU1   <-  read.csv(file = './data/RECODE/PrFixAutoFig_h0.25_s0.01_N100k_U1.csv', header=TRUE)
    AdatN100kU02$PrFixNe  <-  AdatN100kU02$PrFix*(2*AdatN100kU02$N)
    AdatN100kU05$PrFixNe  <-  AdatN100kU05$PrFix*(2*AdatN100kU05$N)
    AdatN100kU1$PrFixNe  <-  AdatN100kU1$PrFix*(AdatN100kU1$N/2)
    AdatN100kU02$PrFix[AdatN100kU02$PrFix < (1/(2*10^8))]  <-  (1/(2*10^8))
    AdatN100kU05$PrFix[AdatN100kU05$PrFix < (1/(2*10^8))]  <-  (1/(2*10^8))
    AdatN100kU1$PrFix[AdatN100kU1$PrFix   < (1/(2*10^8))]  <-  (1/(2*10^8))

    ## N = 1mil data
    AdatN1milU02  <-  read.csv(file = './data/RECODE/PrFixAutoFig_h0.25_s0.01_N1mil_U02.csv', header=TRUE)
    AdatN1milU05  <-  read.csv(file = './data/RECODE/PrFixAutoFig_h0.25_s0.01_N1mil_U05.csv', header=TRUE)
    AdatN1milU1   <-  read.csv(file = './data/RECODE/PrFixAutoFig_h0.25_s0.01_N1mil_U1.csv', header=TRUE)
    AdatN1milU02$PrFixNe  <-  AdatN1milU02$PrFix*(2*AdatN1milU02$N)
    AdatN1milU05$PrFixNe  <-  AdatN1milU05$PrFix*(2*AdatN1milU05$N)
    AdatN1milU1$PrFixNe  <-  AdatN1milU1$PrFix*(AdatN1milU1$N/2)
    AdatN1milU02$PrFix[AdatN1milU02$PrFix < (1/(2*10^8))]  <-  (1/(2*10^8))
    AdatN1milU05$PrFix[AdatN1milU05$PrFix < (1/(2*10^8))]  <-  (1/(2*10^8))
    AdatN1milU1$PrFix[AdatN1milU1$PrFix < (1/(2*10^8))]    <-  (1/(2*10^8))

    ## SLR Expanding Inversion data
    ## N = 100k data
    datN100kU02  <-  read.csv(file = './data/RECODE/PrFixFig_Haploid_h0.25_s0.01_N100k_U02.csv', header=TRUE)
    datN100kU05  <-  read.csv(file = './data/RECODE/PrFixFig_Haploid_h0.25_s0.01_N100k_U05.csv', header=TRUE)
    datN100kU1   <-  read.csv(file = './data/RECODE/PrFixFig_Haploid_h0.25_s0.01_N100k_U1.csv', header=TRUE)
    datN100kU02$PrFixNe  <-  datN100kU02$PrFix*(datN100kU02$N/2)
    datN100kU05$PrFixNe  <-  datN100kU05$PrFix*(datN100kU05$N/2)
    datN100kU1$PrFixNe   <-  datN100kU1$PrFix*(datN100kU1$N/2)
    datN100kU02$PrFix[datN100kU02$PrFix < (2/(10^8))]  <-  (2/(10^8))
    datN100kU05$PrFix[datN100kU05$PrFix < (2/(10^8))]  <-  (2/(10^8))
    datN100kU1$PrFix[datN100kU1$PrFix < (2/(10^8))]  <-  (2/(10^8))
    rdatN100kU02  <-  read.csv(file = './data/RECODE/PrFixFig_Haploid_rFixedInv_h0.25_s0.01_N100k_U02.csv', header=TRUE)
    rdatN100kU05  <-  read.csv(file = './data/RECODE/PrFixFig_Haploid_rFixedInv_h0.25_s0.01_N100k_U05.csv', header=TRUE)
    rdatN100kU1   <-  read.csv(file = './data/RECODE/PrFixFig_Haploid_rFixedInv_h0.25_s0.01_N100k_U1.csv', header=TRUE)
    rdatN100kU02  <-  subset(rdatN100kU02,rdatN100kU02$rFixedInv == 0)
    rdatN100kU05  <-  subset(rdatN100kU05,rdatN100kU05$rFixedInv == 0)
    rdatN100kU1   <-  subset(rdatN100kU1,rdatN100kU1$rFixedInv == 0)
    # Add 0's for Pr(fix | r=0) for table missing table values
    for(i in 1:length(unique(datN100kU02$x))) {
        if(!any(rdatN100kU02$rInvSize == datN100kU02$x[i])) {
            rdatN100kU02  <-  rbind(rdatN100kU02, rdatN100kU02[i-1,])
            rdatN100kU02$rInvSize[i]  <-  datN100kU02$x[i]
            rdatN100kU02$Freq[i]  <-  (2/(10^8))
        }
        if(!any(rdatN100kU05$rInvSize == datN100kU05$x[i])) {
            rdatN100kU05  <-  rbind(rdatN100kU05, rdatN100kU05[i-1,])
            rdatN100kU05$rInvSize[i]  <-  datN100kU05$x[i]
            rdatN100kU05$Freq[i]  <-  (2/(10^8))
        }
        if(!any(rdatN100kU1$rInvSize == datN100kU1$x[i])) {
            rdatN100kU1  <-  rbind(rdatN100kU1, rdatN100kU1[i-1,])
            rdatN100kU1$rInvSize[i]  <-  datN100kU1$x[i]
            rdatN100kU1$Freq[i]  <-  (2/(10^8))
        }
    }


    ## N = 1 mil data
    datN1milU02  <-  read.csv(file = './data/RECODE/PrFixFig_Haploid_h0.25_s0.01_N1mil_U02.csv', header=TRUE)
    datN1milU05  <-  read.csv(file = './data/RECODE/PrFixFig_Haploid_h0.25_s0.01_N1mil_U05.csv', header=TRUE)
    datN1milU1  <-  read.csv(file = './data/RECODE/PrFixFig_Haploid_h0.25_s0.01_N1mil_U1.csv', header=TRUE)
    datN1milU02$PrFixNe  <-  datN1milU02$PrFix*(datN1milU02$N/2)
    datN1milU05$PrFixNe  <-  datN1milU05$PrFix*(datN1milU05$N/2)
    datN1milU1$PrFixNe  <-  datN1milU1$PrFix*(datN1milU1$N/2)
    datN1milU02$PrFix[datN1milU02$PrFix < (2/(10^8))]  <-  (2/(10^8))
    datN1milU05$PrFix[datN1milU05$PrFix < (2/(10^8))]  <-  (2/(10^8))
    datN1milU1$PrFix[datN1milU1$PrFix < (2/(10^8))]  <-  (2/(10^8))
    rdatN1milU02  <-  read.csv(file = './data/RECODE/PrFixFig_Haploid_rFixedInv_h0.25_s0.01_N1mil_U02.csv', header=TRUE)
    rdatN1milU05  <-  read.csv(file = './data/RECODE/PrFixFig_Haploid_rFixedInv_h0.25_s0.01_N1mil_U05.csv', header=TRUE)
    rdatN1milU1  <-  read.csv(file = './data/RECODE/PrFixFig_Haploid_rFixedInv_h0.25_s0.01_N1mil_U1.csv', header=TRUE)
    rdatN1milU02  <-  subset(rdatN1milU02,rdatN1milU02$rFixedInv == 0)
    rdatN1milU05  <-  subset(rdatN1milU05,rdatN1milU05$rFixedInv == 0)
    rdatN1milU1   <-  subset(rdatN1milU1,rdatN1milU1$rFixedInv == 0)
    # Add 0's for Pr(fix | r=0) for table missing table values
    for(i in 1:length(unique(datN1milU02$x))) {
        if(!any(rdatN1milU02$rInvSize == datN1milU02$x[i])) {
            rdatN1milU02  <-  rbind(rdatN1milU02, rdatN1milU02[i-1,])
            rdatN1milU02$rInvSize[i]  <-  datN1milU02$x[i]
            rdatN1milU02$Freq[i]  <-  (2/(10^8))
        }
        if(!any(rdatN1milU05$rInvSize == datN1milU05$x[i])) {
            rdatN1milU05  <-  rbind(rdatN1milU05, rdatN1milU05[i-1,])
            rdatN1milU05$rInvSize[i]  <-  datN1milU05$x[i]
            rdatN1milU05$Freq[i]  <-  (2/(10^8))
        }
        if(!any(rdatN1milU1$rInvSize == datN1milU1$x[i])) {
            rdatN1milU1  <-  rbind(rdatN1milU1, rdatN1milU1[i-1,])
            rdatN1milU1$rInvSize[i]  <-  datN1milU1$x[i]
            rdatN1milU1$Freq[i]  <-  (2/(10^8))
        }
    }

    # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  colfunc(3)

    # set plot layout
    layout.mat <- matrix(c(1:4), nrow=2, ncol=2, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

## Panel A & B -- Autosomal inversions
## Panel A
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(5,4,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c((1/(2*10^8)), (4/(2*10^5))), log='y', ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80', log='y')
        box()
        # Points
        abline(h=(1/(2*10^5)), lwd=1, lty=2, col=1)
        points(PrFix ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=AdatN100kU02)
        points(PrFix ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=AdatN100kU05)
        points(PrFix ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=AdatN100kU1)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, at=c((1/(2*10^8)), 
                            (1/(2*10^7)), 
                            (1/(2*10^6)), 
                            (1/(2*10^5)), 
                            (1/(2*10^4))), 
                    labels=c(0, 
                            expression(1/(2%*%10^7)), 
                            expression(1/(2%*%10^6)), 
                            expression(1/(2%*%10^5)), 
                            expression(1/(2%*%10^4))))
        proportionalLabel(1.15, 1.2, 'Autosomal', cex=1.5, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(-0.0175, 0.16, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(-0.0175, 0.14, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(0.04, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.5, 1.06, expression(italic(N)==10^5) , cex=1.25, adj=c(0.5, 0.5), xpd=NA, log='y')
        # Plot labels etc.
        proportionalLabel(-0.4, 0.5, expression(Fixation~Probability), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90, log='y')


## Panel B
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c((1/(2*10^8)), (4/(2*10^5))), log='y', ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80', log='y')
        box()
        # Points
        abline(h=(1/(2*10^6)), lwd=1, lty=2, col=1)
        abline(h=1, lwd=1, lty=1, col=1)
        points(PrFix ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=AdatN1milU02)
        points(PrFix ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=AdatN1milU05)
        points(PrFix ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=AdatN1milU1)
        # axes
        axis(1, las=1, labels=NA)
         axis(2, las=1, at=c((1/(2*10^8)), 
                             (1/(2*10^7)), 
                             (1/(2*10^6)), 
                             (1/(2*10^5)), 
                             (1/(2*10^4))), 
                labels=NA)
        proportionalLabel(-0.0175, 0.16, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(-0.0175, 0.14, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(0.04, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.5, 1.06, expression(italic(N)==10^6) , cex=1.25, adj=c(0.5, 0.5), xpd=NA, log='y')
        # Plot labels etc.
        # Legend
        legend(
               x       =  1.02,
               y       =  0.00003,
               legend  =  c(
                            expression(paste(italic(U/(hs)), " = ", 8)),
                            expression(paste(italic(U/(hs)), " = ", 20)),
                            expression(paste(italic(U/(hs)), " = ", 40))),
               pch     =  c(21,22,24),
               col     =  transparentColor(COLS[1], opacity=1),
               pt.bg     =  transparentColor(COLS[2], opacity=0.6),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


## Panels C & D -- SLR Expanding Inversions
## Panel C
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c((2/10^8), 8/(10^5)), log='y', ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80', log='y')
        box()
        # Benchmark 2/N
        abline(h=(2/100000), lwd=1, lty=2, col=1)
        # PrFix Points
        points(PrFix ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN100kU02)
        points(PrFix ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN100kU05)
        points(PrFix ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN100kU1)
        # r0 PrFix Points
        points(Freq ~ rInvSize, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rdatN100kU02)
        points(Freq ~ rInvSize, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rdatN100kU05)
        points(Freq ~ rInvSize, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rdatN100kU1)
        # f0*qStar
qStars  <-  qstar(q.0=(2/datN100kU02$N[1]), U=datN100kU02$U[1], x=datN100kU02$x, s=datN100kU02$s[1], h=datN100kU02$h[1])        
p0  <-  exp(-(datN100kU02$U[1] * datN100kU02$x)/(datN100kU02$s[1]*datN100kU02$h[1]))
detExp  <-  p0*qStars
detExp[detExp < 2/(10^8)]  <-  2/(10^8)
points(detExp ~ datN100kU02$x, pch=21)
qStars  <-  qstar(q.0=(2/datN100kU05$N[1]), U=datN100kU05$U[1], x=datN100kU05$x, s=datN100kU05$s[1], h=datN100kU05$h[1])        
p0  <-  exp(-(datN100kU05$U[1] * datN100kU05$x)/(datN100kU05$s[1]*datN100kU05$h[1]))
detExp  <-  p0*qStars
detExp[detExp < 2/(10^8)]  <-  2/(10^8)
points(detExp ~ datN100kU05$x, pch=22)
qStars  <-  qstar(q.0=(2/datN100kU1$N[1]), U=datN100kU1$U[1], x=datN100kU1$x, s=datN100kU1$s[1], h=datN100kU1$h[1])        
p0  <-  exp(-(datN100kU1$U[1] * datN100kU1$x)/(datN100kU1$s[1]*datN100kU1$h[1]))
detExp  <-  p0*qStars
detExp[detExp < 2/(10^8)]  <-  2/(10^8)
points(detExp ~ datN100kU1$x, pch=24)
        # axes
        axis(1, las=1)
        axis(2, las=1, at=c(2/10^8, 
                            2/10^7, 
                            2/10^6, 
                            2/10^5, 
                            2/10^4), 
                    labels=c(0, 
                            expression(2/10^7), 
                            expression(2/10^6), 
                            expression(2/10^5), 
                            expression(2/10^4)))
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(1.15, 1.2, 'Expand SLR', cex=1.5, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.5, 1.06, expression(italic(N)==10^5) , cex=1.25, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(-0.4,  0.5,   expression(paste("Fixation Probability")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90, log='y')        
        proportionalLabel(-0.0175, 0.16, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(-0.0175, 0.14, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel( 0.5,  -0.225,  expression(paste("Inversion size (", italic(x),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')

## Panel D
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c((2/10^8), 8/(10^5)), log='y', ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80', log='y')
        box()
        # Benchmark 2/N
        abline(h=(2/10^6), lwd=1, lty=2, col=1)
        # PrFix Points
        points(PrFix ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN1milU02)
        points(PrFix ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN1milU05)
        points(PrFix ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN1milU1)
        # r0 PrFix Points
        points(Freq ~ rInvSize, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rdatN1milU02)
        points(Freq ~ rInvSize, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rdatN1milU05)
        points(Freq ~ rInvSize, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rdatN1milU1)
        # f0*qStar
qStars  <-  qstar(q.0=(2/datN1milU02$N[1]), U=datN1milU02$U[1], x=datN1milU02$x, s=datN1milU02$s[1], h=datN1milU02$h[1])        
p0  <-  exp(-(datN1milU02$U[1] * datN1milU02$x)/(datN1milU02$s[1]*datN1milU02$h[1]))
detExp  <-  p0*qStars
detExp[detExp < 2/(10^8)]  <-  2/(10^8)
points(detExp ~ datN1milU02$x, pch=21)
qStars  <-  qstar(q.0=(2/datN1milU05$N[1]), U=datN1milU05$U[1], x=datN1milU05$x, s=datN1milU05$s[1], h=datN1milU05$h[1])        
p0  <-  exp(-(datN1milU05$U[1] * datN1milU05$x)/(datN1milU05$s[1]*datN1milU05$h[1]))
detExp  <-  p0*qStars
detExp[detExp < 2/(10^8)]  <-  2/(10^8)
points(detExp ~ datN1milU05$x, pch=22)
qStars  <-  qstar(q.0=(2/datN1milU1$N[1]), U=datN1milU1$U[1], x=datN1milU1$x, s=datN1milU1$s[1], h=datN1milU1$h[1])        
p0  <-  exp(-(datN1milU1$U[1] * datN1milU1$x)/(datN1milU1$s[1]*datN1milU1$h[1]))
detExp  <-  p0*qStars
detExp[detExp < 2/(10^8)]  <-  2/(10^8)
points(detExp ~ datN1milU1$x, pch=24)
        # axes
        axis(1, las=1)
        axis(2, las=1, at=c(2/10^7, 2/10^6, 2/10^5, 2/10^4, 2/10^3, 2/10^2, 2/10^1), labels=NA)
        # Plot labels etc.
        proportionalLabel(0.5, 1.06, expression(italic(N)==10^6) , cex=1.25, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.04, 1.05, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(-0.0175, 0.16, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(-0.0175, 0.14, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel( 0.5,  -0.225,  expression(paste("Inversion size (", italic(x),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        # Legend
        legend(
               x       =  0.88,
               y       =  1.075/(10^4),
               legend  =  c(
                            expression(paste("Large ", italic(N), " approx. (", italic(f)[0],italic(q)^"*", ")")),
                            expression(paste("All fixed inversions (all ", italic(d), ")")),
                            expression(paste("Initially unloaded (",italic(d), " = ", 0, ")"))),
               pch     =  c(24,24,24),
               col     =  transparentColor(COLS[1], opacity=1),
               pt.bg   =  c(NA,
                            transparentColor(COLS[2], opacity=0.6), 
                            transparentColor('tomato', opacity=0.6)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
        legend(
               x       =  0.16,
               y       =  1.1/(10^4),
               legend  =  c(NA,NA,NA),
               pch     =  c(22,22,22),
               col     =  transparentColor(COLS[1], opacity=1),
               pt.bg   =  c(NA,
                            transparentColor(COLS[2], opacity=0.6), 
                            transparentColor('tomato', opacity=0.6)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
        legend(
               x       =  0.095,
               y       =  1.1/(10^4),
               legend  =  c(NA,NA,NA),
               pch     =  c(21,21,21),
               col     =  transparentColor(COLS[1], opacity=1),
               pt.bg   =  c(NA,
                            transparentColor(COLS[2], opacity=0.6), 
                            transparentColor('tomato', opacity=0.6)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

}





#' Fig. D1
#' Fixation probability of different size inversions 
#' expanding the SDR on a Y chromosome: Pr(fix | x) 
#' Illustration of Pr(fix | x) for different 
#' population size and U/s ratios.
PrFixFig4Panel_Exact  <-  function(suppVersion=FALSE) {

    # Import data for plotting

    ## Autosomal Inversion Data
    AdatN25kU02  <-  read.csv(file = './data/RECODE/PrFixFig_AutoExact_h0.25_s0.01_N250k_U02_Nei.csv', header=TRUE)
    AdatN25kU05  <-  read.csv(file = './data/RECODE/PrFixFig_AutoExact_h0.25_s0.01_N250k_U05_Nei.csv', header=TRUE)
    AdatN25kU1   <-  read.csv(file = './data/RECODE/PrFixFig_AutoExact_h0.25_s0.01_N250k_U1_Nei.csv', header=TRUE)
    AdatN25kU02$PrFixNe  <-  AdatN25kU02$PrFix*(AdatN25kU02$N*2)
    AdatN25kU05$PrFixNe  <-  AdatN25kU05$PrFix*(AdatN25kU05$N*2)
    AdatN25kU1$PrFixNe   <-  AdatN25kU1$PrFix*(AdatN25kU1$N*2)
    AdatN25kU02$PrFix[AdatN25kU02$PrFix < (1/(2*10^7))]  <-  (1/(2*10^7))
    AdatN25kU05$PrFix[AdatN25kU05$PrFix < (1/(2*10^7))]  <-  (1/(2*10^7))
    AdatN25kU1$PrFix[AdatN25kU1$PrFix < (1/(2*10^7))]    <-  (1/(2*10^7))
    rAdatN25kU02  <-  read.csv(file = './data/RECODE/PrFixFig_AutoExact_rFixedInv_h0.25_s0.01_N250k_U02_Nei.csv', header=TRUE)
    rAdatN25kU05  <-  read.csv(file = './data/RECODE/PrFixFig_AutoExact_rFixedInv_h0.25_s0.01_N250k_U05_Nei.csv', header=TRUE)
    rAdatN25kU1   <-  read.csv(file = './data/RECODE/PrFixFig_AutoExact_rFixedInv_h0.25_s0.01_N250k_U1_Nei.csv', header=TRUE)
    rAdatN25kU02  <-  subset(rAdatN25kU02,rAdatN25kU02$rFixedInv == 0)
    rAdatN25kU05  <-  subset(rAdatN25kU05,rAdatN25kU05$rFixedInv == 0)
    rAdatN25kU1   <-  subset(rAdatN25kU1,rAdatN25kU1$rFixedInv == 0)
    # Add 0's for Pr(fix | r=0) for table missing table values
    for(i in 1:length(unique(AdatN25kU02$x))) {
        if(!any(rAdatN25kU02$rInvSize == AdatN25kU02$x[i])) {
            rAdatN25kU02  <-  rbind(rAdatN25kU02, rAdatN25kU02[i-1,])
            rAdatN25kU02$rInvSize[i]  <-  AdatN25kU02$x[i]
            rAdatN25kU02$Freq[i]  <-  (1/(2*10^7))
        }
        if(!any(rAdatN25kU05$rInvSize == AdatN25kU05$x[i])) {
            rAdatN25kU05  <-  rbind(rAdatN25kU05, rAdatN25kU05[i-1,])
            rAdatN25kU05$rInvSize[i]  <-  AdatN25kU05$x[i]
            rAdatN25kU05$Freq[i]  <-  (1/(2*10^7))
        }
        if(!any(rAdatN25kU1$rInvSize == AdatN25kU1$x[i])) {
            rAdatN25kU1  <-  rbind(rAdatN25kU1, rAdatN25kU1[i-1,])
            rAdatN25kU1$rInvSize[i]  <-  AdatN25kU1$x[i]
            rAdatN25kU1$Freq[i]  <-  (1/(2*10^7))
        }
    }

    ## N = 1 mil data
    AdatN250kU02  <-  read.csv(file = './data/RECODE/PrFixFig_AutoExact_h0.25_s0.01_N250k_U02_Nei.csv', header=TRUE)
    AdatN250kU05  <-  read.csv(file = './data/RECODE/PrFixFig_AutoExact_h0.25_s0.01_N250k_U05_Nei.csv', header=TRUE)
    AdatN250kU1   <-  read.csv(file = './data/RECODE/PrFixFig_AutoExact_h0.25_s0.01_N250k_U1_Nei.csv', header=TRUE)
    AdatN250kU02$PrFixNe  <-  AdatN250kU02$PrFix*(AdatN250kU02$N*2)
    AdatN250kU05$PrFixNe  <-  AdatN250kU05$PrFix*(AdatN250kU05$N*2)
    AdatN250kU1$PrFixNe   <-  AdatN250kU1$PrFix*(AdatN250kU1$N*2)
    AdatN250kU02$PrFix[AdatN250kU02$PrFix < (1/(2*10^7))]  <-  (1/(2*10^7))
    AdatN250kU05$PrFix[AdatN250kU05$PrFix < (1/(2*10^7))]  <-  (1/(2*10^7))
    AdatN250kU1$PrFix[AdatN250kU1$PrFix < (1/(2*10^7))]    <-  (1/(2*10^7))
    rAdatN250kU02  <-  read.csv(file = './data/RECODE/PrFixFig_AutoExact_rFixedInv_h0.25_s0.01_N250k_U02_Nei.csv', header=TRUE)
    rAdatN250kU05  <-  read.csv(file = './data/RECODE/PrFixFig_AutoExact_rFixedInv_h0.25_s0.01_N250k_U05_Nei.csv', header=TRUE)
    rAdatN250kU1   <-  read.csv(file = './data/RECODE/PrFixFig_AutoExact_rFixedInv_h0.25_s0.01_N250k_U1_Nei.csv', header=TRUE)    
    rAdatN250kU02  <-  subset(rAdatN250kU02,rAdatN250kU02$rFixedInv == 0)
    rAdatN250kU05  <-  subset(rAdatN250kU05,rAdatN250kU05$rFixedInv == 0)
    rAdatN250kU1   <-  subset(rAdatN250kU1,rAdatN250kU1$rFixedInv == 0)
    # Add 0's for Pr(fix | r=0) for table missing table values
    for(i in 1:length(unique(AdatN250kU02$x))) {
        if(!any(rAdatN250kU02$rInvSize == AdatN250kU02$x[i])) {
            rAdatN250kU02  <-  rbind(rAdatN250kU02, rAdatN250kU02[i-1,])
            rAdatN250kU02$rInvSize[i]  <-  AdatN250kU02$x[i]
            rAdatN250kU02$Freq[i]  <-  (1/(2*10^7))
        }
        if(!any(rAdatN250kU05$rInvSize == AdatN250kU05$x[i])) {
            rAdatN250kU05  <-  rbind(rAdatN250kU05, rAdatN250kU05[i-1,])
            rAdatN250kU05$rInvSize[i]  <-  AdatN250kU05$x[i]
            rAdatN250kU05$Freq[i]  <-  (1/(2*10^7))
        }
        if(!any(rAdatN250kU1$rInvSize == AdatN250kU1$x[i])) {
            rAdatN250kU1  <-  rbind(rAdatN250kU1, rAdatN250kU1[i-1,])
            rAdatN250kU1$rInvSize[i]  <-  AdatN250kU1$x[i]
            rAdatN250kU1$Freq[i]  <-  (1/(2*10^7))
        }
    }


    ## SLR Expanding Inversion data
    ## N = 100k data
    datN100kU02  <-  read.csv(file = './data/RECODE/PrFixFig_SLR-Nei_h0.25_s0.01_N100k_U02_Nei.csv', header=TRUE)
    datN100kU05  <-  read.csv(file = './data/RECODE/PrFixFig_SLR-Nei_h0.25_s0.01_N100k_U05_Nei.csv', header=TRUE)
    datN100kU1   <-  read.csv(file = './data/RECODE/PrFixFig_SLR-Nei_h0.25_s0.01_N100k_U1_Nei.csv', header=TRUE)
    datN100kU02$PrFixNe  <-  datN100kU02$PrFix*(datN100kU02$N/2)
    datN100kU05$PrFixNe  <-  datN100kU05$PrFix*(datN100kU05$N/2)
    datN100kU1$PrFixNe   <-  datN100kU1$PrFix*(datN100kU1$N/2)
    datN100kU02$PrFix[datN100kU02$PrFix < (2/(10^8))]  <-  (2/(10^8))
    datN100kU05$PrFix[datN100kU05$PrFix < (2/(10^8))]  <-  (2/(10^8))
    datN100kU1$PrFix[datN100kU1$PrFix < (2/(10^8))]    <-  (2/(10^8))
    rdatN100kU02  <-  read.csv(file = './data/RECODE/PrFixFig_SLR-Nei_rFixedInv_h0.25_s0.01_N100k_U02_Nei.csv', header=TRUE)
    rdatN100kU05  <-  read.csv(file = './data/RECODE/PrFixFig_SLR-Nei_rFixedInv_h0.25_s0.01_N100k_U05_Nei.csv', header=TRUE)
    rdatN100kU1   <-  read.csv(file = './data/RECODE/PrFixFig_SLR-Nei_rFixedInv_h0.25_s0.01_N100k_U1_Nei.csv', header=TRUE)
    rdatN100kU02  <-  subset(rdatN100kU02,rdatN100kU02$rFixedInv == 0)
    rdatN100kU05  <-  subset(rdatN100kU05,rdatN100kU05$rFixedInv == 0)
    rdatN100kU1   <-  subset(rdatN100kU1,rdatN100kU1$rFixedInv == 0)
    # Add 0's for Pr(fix | r=0) for table missing table values
    for(i in 1:length(unique(datN100kU02$x))) {
        if(!any(rdatN100kU02$rInvSize == datN100kU02$x[i])) {
            rdatN100kU02  <-  rbind(rdatN100kU02, rdatN100kU02[i-1,])
            rdatN100kU02$rInvSize[i]  <-  datN100kU02$x[i]
            rdatN100kU02$Freq[i]  <-  (2/(10^8))
        }
        if(!any(rdatN100kU05$rInvSize == datN100kU05$x[i])) {
            rdatN100kU05  <-  rbind(rdatN100kU05, rdatN100kU05[i-1,])
            rdatN100kU05$rInvSize[i]  <-  datN100kU05$x[i]
            rdatN100kU05$Freq[i]  <-  (2/(10^8))
        }
        if(!any(rdatN100kU1$rInvSize == datN100kU1$x[i])) {
            rdatN100kU1  <-  rbind(rdatN100kU1, rdatN100kU1[i-1,])
            rdatN100kU1$rInvSize[i]  <-  datN100kU1$x[i]
            rdatN100kU1$Freq[i]  <-  (2/(10^8))
        }
    }

    ## N = 1 mil data
    datN250kU02  <-  read.csv(file = './data/RECODE/PrFixFig_SLR-Nei_h0.25_s0.01_N250k_U02_Nei.csv', header=TRUE)
    datN250kU05  <-  read.csv(file = './data/RECODE/PrFixFig_SLR-Nei_h0.25_s0.01_N250k_U05_Nei.csv', header=TRUE)
    datN250kU1   <-  read.csv(file = './data/RECODE/PrFixFig_SLR-Nei_h0.25_s0.01_N250k_U1_Nei.csv', header=TRUE)
    datN250kU02$PrFixNe  <-  datN250kU02$PrFix*(datN250kU02$N/2)
    datN250kU05$PrFixNe  <-  datN250kU05$PrFix*(datN250kU05$N/2)
    datN250kU1$PrFixNe   <-  datN250kU1$PrFix*(datN250kU1$N/2)
    datN250kU02$PrFix[datN250kU02$PrFix < (2/(10^8))]  <-  (2/(10^8))
    datN250kU05$PrFix[datN250kU05$PrFix < (2/(10^8))]  <-  (2/(10^8))
    datN250kU1$PrFix[datN250kU1$PrFix   < (2/(10^8))]  <-  (2/(10^8))
    rdatN250kU02  <-  read.csv(file = './data/RECODE/PrFixFig_SLR-Nei_rFixedInv_h0.25_s0.01_N250k_U02_Nei.csv', header=TRUE)
    rdatN250kU05  <-  read.csv(file = './data/RECODE/PrFixFig_SLR-Nei_rFixedInv_h0.25_s0.01_N250k_U05_Nei.csv', header=TRUE)
    rdatN250kU1   <-  read.csv(file = './data/RECODE/PrFixFig_SLR-Nei_rFixedInv_h0.25_s0.01_N250k_U1_Nei.csv', header=TRUE)
    rdatN250kU02  <-  subset(rdatN250kU02,rdatN250kU02$rFixedInv == 0)
    rdatN250kU05  <-  subset(rdatN250kU05,rdatN250kU05$rFixedInv == 0)
    rdatN250kU1   <-  subset(rdatN250kU1,rdatN250kU1$rFixedInv == 0)
    # Add 0's for Pr(fix | r=0) for table missing table values
    for(i in 1:length(unique(datN250kU02$x))) {
        if(!any(rdatN250kU02$rInvSize == datN250kU02$x[i])) {
            rdatN250kU02  <-  rbind(rdatN250kU02, rdatN250kU02[i-1,])
            rdatN250kU02$rInvSize[i]  <-  datN250kU02$x[i]
            rdatN250kU02$Freq[i]  <-  (2/(10^8))
        }
        if(!any(rdatN250kU05$rInvSize == datN250kU05$x[i])) {
            rdatN250kU05  <-  rbind(rdatN250kU05, rdatN250kU05[i-1,])
            rdatN250kU05$rInvSize[i]  <-  datN250kU05$x[i]
            rdatN250kU05$Freq[i]  <-  (2/(10^8))
        }
        if(!any(rdatN250kU1$rInvSize == datN250kU1$x[i])) {
            rdatN250kU1  <-  rbind(rdatN250kU1, rdatN250kU1[i-1,])
            rdatN250kU1$rInvSize[i]  <-  datN250kU1$x[i]
            rdatN250kU1$Freq[i]  <-  (2/(10^8))
        }
    }


    # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  colfunc(3)

    # set plot layout
    layout.mat <- matrix(c(1:4), nrow=2, ncol=2, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

## Panel A & B -- Autosomal inversions
## Panel A
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(5,4,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c((1/(2*10^7)), (4/(2*10^4))), log='y', ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80', log='y')
        box()
        # Points
        abline(h=(1/(2*2.5*10^4)), lwd=1, lty=2, col=1)
        points(PrFix ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=AdatN25kU02)
        points(PrFix ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=AdatN25kU05)
        points(PrFix ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=AdatN25kU1)
        # r0 PrFix Points
        points(Freq ~ rInvSize, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rAdatN25kU02)
        points(Freq ~ rInvSize, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rAdatN25kU05)
        points(Freq ~ rInvSize, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rAdatN25kU1)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, at=c((1/(2*10^7)), 
                            (1/(2*10^6)), 
                            (1/(2*10^5)), 
                            (1/(2*10^4)), 
                            (1/(2*10^3))), 
                    labels=c(0, 
                            expression(1/(2%*%10^7)), 
                            expression(1/(2%*%10^6)), 
                            expression(1/(2%*%10^5)), 
                            expression(1/(2%*%10^4))))
        proportionalLabel(1.15, 1.15, 'Autosomal', cex=1.5, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(-0.0175, 0.16, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(-0.0175, 0.14, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(0.04, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.5, 1.06, expression(italic(N)==2.5%*%10^4) , cex=1.25, adj=c(0.5, 0.5), xpd=NA, log='y')
        # Plot labels etc.
        proportionalLabel(-0.4, 0.5, expression(Fixation~Probability), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90, log='y')


## Panel B
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c((1/(2*10^7)), (4/(2*10^4))), log='y', ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80', log='y')
        box()
        # Points
        abline(h=(1/(2*2.5*10^5)), lwd=1, lty=2, col=1)
        abline(h=1, lwd=1, lty=1, col=1)
        points(PrFix ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=AdatN250kU02)
        points(PrFix ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=AdatN250kU05)
        points(PrFix ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=AdatN250kU1)
        # r0 PrFix Points
        points(Freq ~ rInvSize, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rAdatN250kU02)
        points(Freq ~ rInvSize, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rAdatN250kU05)
        points(Freq ~ rInvSize, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rAdatN250kU1)
        # axes
        axis(1, las=1, labels=NA)
         axis(2, las=1, at=c((1/(2*10^7)), (1/(2*10^6)), (1/(2*10^5)), (1/(2*10^4)), (1/(2*10^3))), 
                labels=NA)
        proportionalLabel(-0.0175, 0.16, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(-0.0175, 0.14, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(0.04, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.5, 1.06, expression(italic(N)==2.5%*%10^5) , cex=1.25, adj=c(0.5, 0.5), xpd=NA, log='y')
        # Plot labels etc.
        # Legend
        legend(
               x       =  1.02,
               y       =  0.0003,
               legend  =  c(
                            expression(paste(italic(U/(hs)), " = ", 8)),
                            expression(paste(italic(U/(hs)), " = ", 20)),
                            expression(paste(italic(U/(hs)), " = ", 40))),
               pch     =  c(21,22,24),
               col     =  transparentColor(COLS[1], opacity=1),
               pt.bg     =  transparentColor(COLS[2], opacity=0.6),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


## Panels C & D -- SLR Expanding Inversions
## Panel C
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c((2/10^8), 4/(10^4)), log='y', ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80', log='y')
        box()
        # Benchmark 2/N
        abline(h=(2/10^5), lwd=1, lty=2, col=1)
        # PrFix Points
        points(PrFix ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN100kU02)
        points(PrFix ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN100kU05)
        points(PrFix ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN100kU1)
        # r0 PrFix Points
        points(Freq ~ rInvSize, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rdatN100kU02)
        points(Freq ~ rInvSize, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rdatN100kU05)
        points(Freq ~ rInvSize, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rdatN100kU1)
        # f0*qStar
        if(suppVersion){
            qStars  <-  qstar(q.0=(2/datN100kU02$N[1]), U=datN100kU02$U[1], x=datN100kU02$x, s=datN100kU02$s[1], h=datN100kU02$h[1])        
            p0  <-  exp(-(datN100kU02$U[1] * datN100kU02$x)/(datN100kU02$s[1]*datN100kU02$h[1]))
            points((p0*qStars) ~ datN100kU02$x, pch=21)
            qStars  <-  qstar(q.0=(2/datN100kU05$N[1]), U=datN100kU05$U[1], x=datN100kU05$x, s=datN100kU05$s[1], h=datN100kU05$h[1])        
            p0  <-  exp(-(datN100kU05$U[1] * datN100kU05$x)/(datN100kU05$s[1]*datN100kU05$h[1]))
            points((p0*qStars) ~ datN100kU05$x, pch=22)
            qStars  <-  qstar(q.0=(2/datN100kU1$N[1]), U=datN100kU1$U[1], x=datN100kU1$x, s=datN100kU1$s[1], h=datN100kU1$h[1])        
            p0  <-  exp(-(datN100kU1$U[1] * datN100kU1$x)/(datN100kU1$s[1]*datN100kU1$h[1]))
            points((p0*qStars) ~ datN100kU1$x, pch=24)
        }
        # axes
        axis(1, las=1)
        axis(2, las=1, at=c(2/10^8, 
                            2/10^7, 
                            2/10^6, 
                            2/10^5, 
                            2/10^4), 
                    labels=c(0, 
                            expression(2/10^7), 
                            expression(2/10^6), 
                            expression(2/10^5), 
                            expression(2/10^4)))
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(1.15, 1.15, 'Expand SLR', cex=1.5, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.5, 1.06, expression(italic(N)==10^5) , cex=1.25, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(-0.4,  0.5,   expression(paste("Fixation Probability")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90, log='y')        
        proportionalLabel(-0.0175, 0.16, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(-0.0175, 0.14, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel( 0.5,  -0.225,  expression(paste("Inversion size (", italic(x),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')

## Panel D
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c((2/10^8), 4/(10^4)), log='y', ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80', log='y')
        box()
        # Benchmark 2/N
        abline(h=(2/(2.5*10^5)), lwd=1, lty=2, col=1)
        # PrFix Points
        points(PrFix ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN250kU02)
        points(PrFix ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN250kU05)
        points(PrFix ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN250kU1)
        # r0 PrFix Points
        points(Freq ~ rInvSize, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rdatN250kU02)
        points(Freq ~ rInvSize, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rdatN250kU05)
        points(Freq ~ rInvSize, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rdatN250kU1)
        if(suppVersion){
        # f0*qStar
            qStars  <-  qstar(q.0=(2/datN250kU02$N[1]), U=datN250kU02$U[1], x=datN250kU02$x, s=datN250kU02$s[1], h=datN250kU02$h[1])        
            p0  <-  exp(-(datN250kU02$U[1] * datN250kU02$x)/(datN250kU02$s[1]*datN250kU02$h[1]))
            points((p0*qStars) ~ datN250kU02$x, pch=21)
            qStars  <-  qstar(q.0=(2/datN250kU05$N[1]), U=datN250kU05$U[1], x=datN250kU05$x, s=datN250kU05$s[1], h=datN250kU05$h[1])        
            p0  <-  exp(-(datN250kU05$U[1] * datN250kU05$x)/(datN250kU05$s[1]*datN250kU05$h[1]))
            points((p0*qStars) ~ datN250kU05$x, pch=22)
            qStars  <-  qstar(q.0=(2/datN250kU1$N[1]), U=datN250kU1$U[1], x=datN250kU1$x, s=datN250kU1$s[1], h=datN250kU1$h[1])        
            p0  <-  exp(-(datN250kU1$U[1] * datN250kU1$x)/(datN250kU1$s[1]*datN250kU1$h[1]))
            points((p0*qStars) ~ datN250kU1$x, pch=24)
        }
        # axes
        axis(1, las=1)
        axis(2, las=1, at=c(2/10^8, 
                            2/10^7, 
                            2/10^6, 
                            2/10^5, 
                            2/10^4), 
                    labels=NA)
        # Plot labels etc.
        proportionalLabel(0.5, 1.06, expression(italic(N)==2.5%*%10^5) , cex=1.25, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.04, 1.05, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(-0.0175, 0.16, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(-0.0175, 0.14, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel( 0.5,  -0.225,  expression(paste("Inversion size (", italic(x),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        if(!suppVersion){
                # Legend
                legend(
                       x       =  1.02,
                       y       =  6/(10^4),
                       legend  =  c(
                                    expression(paste("All ", italic(d))),
                                    expression(paste(italic(d), " = ", 0))),
                       pch     =  c(24,24),
                       col     =  transparentColor(COLS[1], opacity=1),
                       pt.bg   =  c(transparentColor(COLS[2], opacity=0.6), 
                                    transparentColor('tomato', opacity=0.6)),
                       cex     =  1,
                       xjust   =  1,
                       yjust   =  1,
                       bty     =  'n',
                       border  =  NA
                       )
                legend(
                       x       =  0.81,
                       y       =  6/(10^4),
                       legend  =  c(NA,NA),
                       pch     =  c(22,22),
                       col     =  transparentColor(COLS[1], opacity=1),
                       pt.bg   =  c(transparentColor(COLS[2], opacity=0.6), 
                                    transparentColor('tomato', opacity=0.6)),
                       cex     =  1,
                       xjust   =  1,
                       yjust   =  1,
                       bty     =  'n',
                       border  =  NA
                       )
                legend(
                       x       =  0.745,
                       y       =  6/(10^4),
                       legend  =  c(NA,NA,NA),
                       pch     =  c(21,21, NA),
                       col     =  transparentColor(COLS[1], opacity=1),
                       pt.bg   =  c(transparentColor(COLS[2], opacity=0.6), 
                                    transparentColor('tomato', opacity=0.6)),
                       cex     =  1,
                       xjust   =  1,
                       yjust   =  1,
                       bty     =  'n',
                       border  =  NA
                       )
        }
        if(suppVersion){
                # Legend
                legend(
                       x       =  1.02,
                       y       =  6/(10^4),
                       legend  =  c(
                                    expression(paste("All ", italic(d))),
                                    expression(paste(italic(d), " = ", 0)),
                                    expression(paste(italic(f)[0],italic(q)^"*"))),
                       pch     =  c(24,24,24),
                       col     =  transparentColor(COLS[1], opacity=1),
                       pt.bg   =  c(transparentColor(COLS[2], opacity=0.6), 
                                    transparentColor('tomato', opacity=0.6),
                                    NA),
                       cex     =  1,
                       xjust   =  1,
                       yjust   =  1,
                       bty     =  'n',
                       border  =  NA
                       )
                legend(
                       x       =  0.81,
                       y       =  6/(10^4),
                       legend  =  c(NA,NA,NA),
                       pch     =  c(22,22,22),
                       col     =  transparentColor(COLS[1], opacity=1),
                       pt.bg   =  c(transparentColor(COLS[2], opacity=0.6), 
                                    transparentColor('tomato', opacity=0.6),
                                    NA),
                       cex     =  1,
                       xjust   =  1,
                       yjust   =  1,
                       bty     =  'n',
                       border  =  NA
                       )
                legend(
                       x       =  0.745,
                       y       =  6/(10^4),
                       legend  =  c(NA,NA,NA),
                       pch     =  c(21,21,21),
                       col     =  transparentColor(COLS[1], opacity=1),
                       pt.bg   =  c(transparentColor(COLS[2], opacity=0.6), 
                                    transparentColor('tomato', opacity=0.6),
                                    NA),
                       cex     =  1,
                       xjust   =  1,
                       yjust   =  1,
                       bty     =  'n',
                       border  =  NA
                       )
        }
}





####################################
#' Pr(fix | x) Supplement Figure
#' 


#' Fixation probability of different size inversions 
#' expanding the SDR on a Y chromosome: Pr(fix | x) 
#' Illustration of Pr(fix | x) for different 
#' population size and U/s ratios.
PrFixFig4Panel_SuppFig <-  function() {

    # Import data for plotting
    
    ## Autosomal Haploid Model
    AdatHapN100kU02  <-  read.csv(file = './data/RECODE/PrFixAutoFig_h0.1_s0.01_N100k_U02.csv', header=TRUE)
    AdatHapN100kU05  <-  read.csv(file = './data/RECODE/PrFixAutoFig_h0.1_s0.01_N100k_U05.csv', header=TRUE)
    AdatHapN100kU1   <-  read.csv(file = './data/RECODE/PrFixAutoFig_h0.1_s0.01_N100k_U1.csv', header=TRUE)
    AdatHapN100kU02$PrFixNe  <-  AdatHapN100kU02$PrFix*(AdatHapN100kU02$N*2)
    AdatHapN100kU05$PrFixNe  <-  AdatHapN100kU05$PrFix*(AdatHapN100kU05$N*2)
    AdatHapN100kU1$PrFixNe   <-  AdatHapN100kU1$PrFix*(AdatHapN100kU1$N*2)
    AdatHapN100kU02$PrFix[AdatHapN100kU02$PrFix < (1/(2*10^7))]  <-  (1/(2*10^7))
    AdatHapN100kU05$PrFix[AdatHapN100kU05$PrFix < (1/(2*10^7))]  <-  (1/(2*10^7))
    AdatHapN100kU1$PrFix[AdatHapN100kU1$PrFix   < (1/(2*10^7))]    <-  (1/(2*10^7))

    ## Autosomal Expanded Model
    ## N = 1 mil data
    AdatN250kU02  <-  read.csv(file = './data/RECODE/PrFixFig_AutoExact_h0.1_s0.01_N250k_U02_Nei.csv', header=TRUE)
    AdatN250kU05  <-  read.csv(file = './data/RECODE/PrFixFig_AutoExact_h0.1_s0.01_N250k_U05_Nei.csv', header=TRUE)
    AdatN250kU1   <-  read.csv(file = './data/RECODE/PrFixFig_AutoExact_h0.1_s0.01_N250k_U1_Nei.csv', header=TRUE)
    AdatN250kU02$PrFixNe  <-  AdatN250kU02$PrFix*(AdatN250kU02$N*2)
    AdatN250kU05$PrFixNe  <-  AdatN250kU05$PrFix*(AdatN250kU05$N*2)
    AdatN250kU1$PrFixNe   <-  AdatN250kU1$PrFix*(AdatN250kU1$N*2)
    AdatN250kU02$PrFix[AdatN250kU02$PrFix < (1/(2*10^7))]  <-  (1/(2*10^7))
    AdatN250kU05$PrFix[AdatN250kU05$PrFix < (1/(2*10^7))]  <-  (1/(2*10^7))
    AdatN250kU1$PrFix[AdatN250kU1$PrFix < (1/(2*10^7))]    <-  (1/(2*10^7))
    rAdatN250kU02  <-  read.csv(file = './data/RECODE/PrFixFig_AutoExact_rFixedInv_h0.1_s0.01_N250k_U02_Nei.csv', header=TRUE)
    rAdatN250kU05  <-  read.csv(file = './data/RECODE/PrFixFig_AutoExact_rFixedInv_h0.1_s0.01_N250k_U05_Nei.csv', header=TRUE)
    rAdatN250kU1   <-  read.csv(file = './data/RECODE/PrFixFig_AutoExact_rFixedInv_h0.1_s0.01_N250k_U1_Nei.csv', header=TRUE)    
    rAdatN250kU02  <-  subset(rAdatN250kU02,rAdatN250kU02$rFixedInv == 0)
    rAdatN250kU05  <-  subset(rAdatN250kU05,rAdatN250kU05$rFixedInv == 0)
    rAdatN250kU1   <-  subset(rAdatN250kU1,rAdatN250kU1$rFixedInv == 0)
    # Add 0's for Pr(fix | r=0) for table missing table values
    for(i in 1:length(unique(AdatN250kU02$x))) {
        if(!any(rAdatN250kU02$rInvSize == AdatN250kU02$x[i])) {
            rAdatN250kU02  <-  rbind(rAdatN250kU02, rAdatN250kU02[i-1,])
            rAdatN250kU02$rInvSize[i]  <-  AdatN250kU02$x[i]
            rAdatN250kU02$Freq[i]  <-  (1/(2*10^7))
        }
        if(!any(rAdatN250kU05$rInvSize == AdatN250kU05$x[i])) {
            rAdatN250kU05  <-  rbind(rAdatN250kU05, rAdatN250kU05[i-1,])
            rAdatN250kU05$rInvSize[i]  <-  AdatN250kU05$x[i]
            rAdatN250kU05$Freq[i]  <-  (1/(2*10^7))
        }
        if(!any(rAdatN250kU1$rInvSize == AdatN250kU1$x[i])) {
            rAdatN250kU1  <-  rbind(rAdatN250kU1, rAdatN250kU1[i-1,])
            rAdatN250kU1$rInvSize[i]  <-  AdatN250kU1$x[i]
            rAdatN250kU1$Freq[i]  <-  (1/(2*10^7))
        }
    }

    ## SLR Expanding Inversions Haploid Model
    ## N = 1 mil data
    datHapN1milU02  <-  read.csv(file = './data/RECODE/PrFixFig_Haploid_h0.1_s0.01_N1mil_U02.csv', header=TRUE)
    datHapN1milU05  <-  read.csv(file = './data/RECODE/PrFixFig_Haploid_h0.1_s0.01_N1mil_U05.csv', header=TRUE)
    datHapN1milU1  <-  read.csv(file = './data/RECODE/PrFixFig_Haploid_h0.1_s0.01_N1mil_U1.csv', header=TRUE)
    datHapN1milU02$PrFixNe  <-  datHapN1milU02$PrFix*(datHapN1milU02$N/2)
    datHapN1milU05$PrFixNe  <-  datHapN1milU05$PrFix*(datHapN1milU05$N/2)
    datHapN1milU1$PrFixNe  <-  datHapN1milU1$PrFix*(datHapN1milU1$N/2)
    datHapN1milU02$PrFix[datHapN1milU02$PrFix < (2/(10^8))]  <-  (2/(10^8))
    datHapN1milU05$PrFix[datHapN1milU05$PrFix < (2/(10^8))]  <-  (2/(10^8))
    datHapN1milU1$PrFix[datHapN1milU1$PrFix < (2/(10^8))]  <-  (2/(10^8))
    rdatHapN1milU02  <-  read.csv(file = './data/RECODE/PrFixFig_Haploid_rFixedInv_h0.1_s0.01_N1mil_U02.csv', header=TRUE)
    rdatHapN1milU05  <-  read.csv(file = './data/RECODE/PrFixFig_Haploid_rFixedInv_h0.1_s0.01_N1mil_U05.csv', header=TRUE)
    rdatHapN1milU1  <-  read.csv(file = './data/RECODE/PrFixFig_Haploid_rFixedInv_h0.1_s0.01_N1mil_U1.csv', header=TRUE)
    rdatHapN1milU02  <-  subset(rdatHapN1milU02,rdatHapN1milU02$rFixedInv == 0)
    rdatHapN1milU05  <-  subset(rdatHapN1milU05,rdatHapN1milU05$rFixedInv == 0)
    rdatHapN1milU1   <-  subset(rdatHapN1milU1,rdatHapN1milU1$rFixedInv == 0)
    # Add 0's for Pr(fix | r=0) for table missing table values
    for(i in 1:length(unique(datHapN1milU02$x))) {
        if(!any(rdatHapN1milU02$rInvSize == datHapN1milU02$x[i])) {
            rdatHapN1milU02  <-  rbind(rdatHapN1milU02, rdatHapN1milU02[i-1,])
            rdatHapN1milU02$rInvSize[i]  <-  datHapN1milU02$x[i]
            rdatHapN1milU02$Freq[i]  <-  (2/(10^8))
        }
        if(!any(rdatHapN1milU05$rInvSize == datHapN1milU05$x[i])) {
            rdatHapN1milU05  <-  rbind(rdatHapN1milU05, rdatHapN1milU05[i-1,])
            rdatHapN1milU05$rInvSize[i]  <-  datHapN1milU05$x[i]
            rdatHapN1milU05$Freq[i]  <-  (2/(10^8))
        }
        if(!any(rdatHapN1milU1$rInvSize == datHapN1milU1$x[i])) {
            rdatHapN1milU1  <-  rbind(rdatHapN1milU1, rdatHapN1milU1[i-1,])
            rdatHapN1milU1$rInvSize[i]  <-  datHapN1milU1$x[i]
            rdatHapN1milU1$Freq[i]  <-  (2/(10^8))
        }
    }

    ## SLR Expanding Inversions Expanded Model
    ## N = 100k data
    datN250kU02  <-  read.csv(file = './data/RECODE/PrFixFig_SLR-Nei_h0.1_s0.01_N250k_U02_Nei.csv', header=TRUE)
    datN250kU05  <-  read.csv(file = './data/RECODE/PrFixFig_SLR-Nei_h0.1_s0.01_N250k_U05_Nei.csv', header=TRUE)
    datN250kU1   <-  read.csv(file = './data/RECODE/PrFixFig_SLR-Nei_h0.1_s0.01_N250k_U1_Nei.csv', header=TRUE)
    datN250kU02$PrFixNe  <-  datN250kU02$PrFix*(datN250kU02$N/2)
    datN250kU05$PrFixNe  <-  datN250kU05$PrFix*(datN250kU05$N/2)
    datN250kU1$PrFixNe   <-  datN250kU1$PrFix*(datN250kU1$N/2)
    datN250kU02$PrFix[datN250kU02$PrFix < (2/(10^8))]  <-  (2/(10^8))
    datN250kU05$PrFix[datN250kU05$PrFix < (2/(10^8))]  <-  (2/(10^8))
    datN250kU1$PrFix[datN250kU1$PrFix < (2/(10^8))]  <-  (2/(10^8))
    rdatN250kU02  <-  read.csv(file = './data/RECODE/PrFixFig_SLR-Nei_rFixedInv_h0.1_s0.01_N250k_U02_Nei.csv', header=TRUE)
    rdatN250kU05  <-  read.csv(file = './data/RECODE/PrFixFig_SLR-Nei_rFixedInv_h0.1_s0.01_N250k_U05_Nei.csv', header=TRUE)
    rdatN250kU1   <-  read.csv(file = './data/RECODE/PrFixFig_SLR-Nei_rFixedInv_h0.1_s0.01_N250k_U1_Nei.csv', header=TRUE)
    rdatN250kU02  <-  subset(rdatN250kU02,rdatN250kU02$rFixedInv == 0)
    rdatN250kU05  <-  subset(rdatN250kU05,rdatN250kU05$rFixedInv == 0)
    rdatN250kU1   <-  subset(rdatN250kU1,rdatN250kU1$rFixedInv == 0)
    # Add 0's for Pr(fix | r=0) for table missing table values
    for(i in 1:length(unique(datN250kU02$x))) {
        if(!any(rdatN250kU02$rInvSize == datN250kU02$x[i])) {
            rdatN250kU02  <-  rbind(rdatN250kU02, rdatN250kU02[i-1,])
            rdatN250kU02$rInvSize[i]  <-  datN250kU02$x[i]
            rdatN250kU02$Freq[i]  <-  (2/(10^8))
        }
        if(!any(rdatN250kU05$rInvSize == datN250kU05$x[i])) {
            rdatN250kU05  <-  rbind(rdatN250kU05, rdatN250kU05[i-1,])
            rdatN250kU05$rInvSize[i]  <-  datN250kU05$x[i]
            rdatN250kU05$Freq[i]  <-  (2/(10^8))
        }
        if(!any(rdatN250kU1$rInvSize == datN250kU1$x[i])) {
            rdatN250kU1  <-  rbind(rdatN250kU1, rdatN250kU1[i-1,])
            rdatN250kU1$rInvSize[i]  <-  datN250kU1$x[i]
            rdatN250kU1$Freq[i]  <-  (2/(10^8))
        }
    }


    # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  colfunc(3)

    # set plot layout
    layout.mat <- matrix(c(1:4), nrow=2, ncol=2, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

## Panel A & B -- Autosomal inversions
## Panel A
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(5,4,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c((1/(2*10^7)), (4/(2*10^4))), log='y', ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80', log='y')
        box()
        # Points
        abline(h=(1/(2*10^5)), lwd=1, lty=2, col=1)
        points(PrFix ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=AdatHapN100kU02)
        points(PrFix ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=AdatHapN100kU05)
        points(PrFix ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=AdatHapN100kU1)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, at=c((1/(2*10^7)), 
                            (1/(2*10^6)), 
                            (1/(2*10^5)), 
                            (1/(2*10^4)), 
                            (1/(2*10^3))), 
                    labels=c(0, 
                            expression(1/(2%*%10^6)), 
                            expression(1/(2%*%10^5)), 
                            expression(1/(2%*%10^4)), 
                            expression(1/(2%*%10^3))))
        proportionalLabel(1.15, 1.3, 'Autosomal', cex=1.75, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.5, 1.15, 'Equil. Approx. Model', cex=1.5, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(-0.0175, 0.16, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(-0.0175, 0.14, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(0.04, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.5, 1.06, expression(italic(N)==10^5) , cex=1.25, adj=c(0.5, 0.5), xpd=NA, log='y')
        # Plot labels etc.
        proportionalLabel(-0.4, 0.5, expression(Fixation~Probability), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90, log='y')


## Panel B
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c((1/(2*10^7)), (4/(2*10^4))), log='y', ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80', log='y')
        box()
        # Points
        abline(h=(1/(2*2.5*10^5)), lwd=1, lty=2, col=1)
        points(PrFix ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=AdatN250kU02)
        points(PrFix ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=AdatN250kU05)
        points(PrFix ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=AdatN250kU1)
        # r0 PrFix Points
        points(Freq ~ rInvSize, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rAdatN250kU02)
        points(Freq ~ rInvSize, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rAdatN250kU05)
        points(Freq ~ rInvSize, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rAdatN250kU1)
        # axes
        axis(1, las=1, labels=NA)
#        axis(2, las=1)
         axis(2, las=1, at=c((1/(2*10^7)), (1/(2*10^6)), (1/(2*10^5)), (1/(2*10^4)), (1/(2*10^3))), 
                labels=NA)
        proportionalLabel(-0.0175, 0.16, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(-0.0175, 0.14, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(0.5, 1.15, 'Expanded Model', cex=1.5, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.04, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.5, 1.06, expression(italic(N)==2.5%*%10^5) , cex=1.25, adj=c(0.5, 0.5), xpd=NA, log='y')
        # Plot labels etc.
        # Legend
        legend(
               x       =  1.02,
               y       =  0.0003,
               legend  =  c(
                            expression(paste(italic(U/(hs)), " = ", 8)),
                            expression(paste(italic(U/(hs)), " = ", 20)),
                            expression(paste(italic(U/(hs)), " = ", 40))),
               pch     =  c(21,22,24),
               col     =  transparentColor(COLS[1], opacity=1),
               pt.bg     =  transparentColor(COLS[2], opacity=0.6),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


## Panels C & D -- SLR Expanding Inversions
## Panel C
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c((2/10^8), 8/(10^5)), log='y', ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80', log='y')
        box()
        # Benchmark 2/N
        abline(h=(2/1000000), lwd=1, lty=2, col=1)
        # PrFix Points
        points(PrFix ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datHapN1milU02)
        points(PrFix ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datHapN1milU05)
        points(PrFix ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datHapN1milU1)
        # r0 PrFix Points
#        points(Freq ~ rInvSize, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rdatHapN1milU02)
#        points(Freq ~ rInvSize, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rdatHapN1milU05)
#        points(Freq ~ rInvSize, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rdatHapN1milU1)
        # f0*qStar
qStars  <-  qstar(q.0=(2/datHapN1milU02$N[1]), U=datHapN1milU02$U[1], x=datHapN1milU02$x, s=datHapN1milU02$s[1], h=datHapN1milU02$h[1])        
p0  <-  exp(-(datHapN1milU02$U[1] * datHapN1milU02$x)/(datHapN1milU02$s[1]*datHapN1milU02$h[1]))
detExp  <-  p0*qStars
detExp[detExp < 2/(10^8)]  <-  2/(10^8)
points(detExp ~ datHapN1milU02$x, pch=21)
qStars  <-  qstar(q.0=(2/datHapN1milU05$N[1]), U=datHapN1milU05$U[1], x=datHapN1milU05$x, s=datHapN1milU05$s[1], h=datHapN1milU05$h[1])        
p0  <-  exp(-(datHapN1milU05$U[1] * datHapN1milU05$x)/(datHapN1milU05$s[1]*datHapN1milU05$h[1]))
detExp  <-  p0*qStars
detExp[detExp < 2/(10^8)]  <-  2/(10^8)
points(detExp ~ datHapN1milU05$x, pch=22)
qStars  <-  qstar(q.0=(2/datHapN1milU1$N[1]), U=datHapN1milU1$U[1], x=datHapN1milU1$x, s=datHapN1milU1$s[1], h=datHapN1milU1$h[1])        
p0  <-  exp(-(datHapN1milU1$U[1] * datHapN1milU1$x)/(datHapN1milU1$s[1]*datHapN1milU1$h[1]))
detExp  <-  p0*qStars
detExp[detExp < 2/(10^8)]  <-  2/(10^8)
points(detExp ~ datHapN1milU1$x, pch=24)
        # axes
        axis(1, las=1)
        axis(2, las=1, at=c(2/10^8, 
                            2/10^7, 
                            2/10^6, 
                            2/10^5, 
                            2/10^4), 
                    labels=c(0, 
                            expression(2/10^7), 
                            expression(2/10^6), 
                            expression(2/10^5), 
                            expression(2/10^4)))
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(1.15, 1.2, 'Expand SLR', cex=1.75, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.5, 1.06, expression(italic(N)==10^6) , cex=1.25, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(-0.4,  0.5,   expression(paste("Fixation Probability")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90, log='y')        
        proportionalLabel(-0.0175, 0.16, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(-0.0175, 0.14, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel( 0.5,  -0.225,  expression(paste("Inversion size (", italic(x),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')

## Panel D
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
#head(datN1milU05)
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c((2/10^8), 8/(10^5)), log='y', ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80', log='y')
        box()
        # Benchmark 2/N
        abline(h=(2/(2.5*10^5)), lwd=1, lty=2, col=1)
        # PrFix Points
        points(PrFix ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN250kU02)
        points(PrFix ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN250kU05)
        points(PrFix ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN250kU1)
        # r0 PrFix Points
        points(Freq ~ rInvSize, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rdatN250kU02)
        points(Freq ~ rInvSize, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rdatN250kU05)
        points(Freq ~ rInvSize, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rdatN250kU1)
        # f0*qStar
 qStars  <-  qstar(q.0=(2/datN250kU02$N[1]), U=datN250kU02$U[1], x=datN250kU02$x, s=datN250kU02$s[1], h=datN250kU02$h[1])        
 p0  <-  exp(-(datN250kU02$U[1] * datN250kU02$x)/(datN250kU02$s[1]*datN250kU02$h[1]))
 detExp  <-  p0*qStars
 detExp[detExp < 2/(10^8)]  <-  2/(10^8)
 points(detExp ~ datN250kU02$x, pch=21)
 qStars  <-  qstar(q.0=(2/datN250kU05$N[1]), U=datN250kU05$U[1], x=datN250kU05$x, s=datN250kU05$s[1], h=datN250kU05$h[1])        
 p0  <-  exp(-(datN250kU05$U[1] * datN250kU05$x)/(datN250kU05$s[1]*datN250kU05$h[1]))
 detExp  <-  p0*qStars
 detExp[detExp < 2/(10^8)]  <-  2/(10^8)
 points(detExp ~ datN250kU05$x, pch=22)
 qStars  <-  qstar(q.0=(2/datN250kU1$N[1]), U=datN250kU1$U[1], x=datN250kU1$x, s=datN250kU1$s[1], h=datN250kU1$h[1])        
 p0  <-  exp(-(datN250kU1$U[1] * datN250kU1$x)/(datN250kU1$s[1]*datN250kU1$h[1]))
 detExp  <-  p0*qStars
 detExp[detExp < 2/(10^8)]  <-  2/(10^8)
 points(detExp ~ datN250kU1$x, pch=24)
        # axes
        axis(1, las=1)
#        axis(2, las=1, at=c(2/10^6, 2/10^5, 2/10^4, 2/10^3, 2/10^2, 2/10^1), labels=NA)
        axis(2, las=1, at=c(2/10^7, 2/10^6, 2/10^5, 2/10^4, 2/10^3, 2/10^2, 2/10^1), labels=NA)
        # Plot labels etc.
        proportionalLabel(0.5, 1.06, expression(italic(N)==2.5%*%10^5) , cex=1.25, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.04, 1.05, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(-0.0175, 0.16, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(-0.0175, 0.14, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel( 0.5,  -0.225,  expression(paste("Inversion size (", italic(x),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        # Legend
        legend(
               x       =  1.02,
               y       =  1.075/(10^4),
               legend  =  c(
                            expression(paste("All ", italic(d))),
                            expression(paste(italic(d), " = ", 0)),
                            expression(paste(italic(f)[0],italic(q)^"*"))),
               pch     =  c(24,24,24),
               col     =  transparentColor(COLS[1], opacity=1),
               pt.bg   =  c(transparentColor(COLS[2], opacity=0.6), 
                            transparentColor('tomato', opacity=0.6),
                            NA),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
        legend(
               x       =  0.81,
               y       =  1.1/(10^4),
               legend  =  c(NA,NA,NA),
               pch     =  c(22,22,22),
               col     =  transparentColor(COLS[1], opacity=1),
               pt.bg   =  c(transparentColor(COLS[2], opacity=0.6), 
                            transparentColor('tomato', opacity=0.6),
                            NA),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
        legend(
               x       =  0.745,
               y       =  1.1/(10^4),
               legend  =  c(NA,NA,NA),
               pch     =  c(21,21,21),
               col     =  transparentColor(COLS[1], opacity=1),
               pt.bg   =  c(transparentColor(COLS[2], opacity=0.6), 
                            transparentColor('tomato', opacity=0.6),
                            NA),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

}






PrFixFig2Panel_SuppFig <-  function() {

    # Import data for plotting
    
    ## Autosomal Haploid Model
    AdatHapN100kU02  <-  read.csv(file = './data/RECODE/PrFixAutoFig_h0.1_s0.01_N100k_U02.csv', header=TRUE)
    AdatHapN100kU05  <-  read.csv(file = './data/RECODE/PrFixAutoFig_h0.1_s0.01_N100k_U05.csv', header=TRUE)
    AdatHapN100kU1   <-  read.csv(file = './data/RECODE/PrFixAutoFig_h0.1_s0.01_N100k_U1.csv', header=TRUE)
    AdatHapN100kU02$PrFixNe  <-  AdatHapN100kU02$PrFix*(AdatHapN100kU02$N*2)
    AdatHapN100kU05$PrFixNe  <-  AdatHapN100kU05$PrFix*(AdatHapN100kU05$N*2)
    AdatHapN100kU1$PrFixNe   <-  AdatHapN100kU1$PrFix*(AdatHapN100kU1$N*2)
    AdatHapN100kU02$PrFix[AdatHapN100kU02$PrFix < (1/(2*10^8))]  <-  (1/(2*10^8))
    AdatHapN100kU05$PrFix[AdatHapN100kU05$PrFix < (1/(2*10^8))]  <-  (1/(2*10^8))
    AdatHapN100kU1$PrFix[AdatHapN100kU1$PrFix   < (1/(2*10^8))]    <-  (1/(2*10^8))


    ## SLR Expanding Inversions Haploid Model
    ## N = 1 mil data
    datHapN100kU02  <-  read.csv(file = './data/RECODE/PrFixFig_Haploid_h0.1_s0.01_N1mil_U02.csv', header=TRUE)
    datHapN100kU05  <-  read.csv(file = './data/RECODE/PrFixFig_Haploid_h0.1_s0.01_N1mil_U05.csv', header=TRUE)
    datHapN100kU1  <-  read.csv(file = './data/RECODE/PrFixFig_Haploid_h0.1_s0.01_N1mil_U1.csv', header=TRUE)
    datHapN100kU02$PrFixNe  <-  datHapN100kU02$PrFix*(datHapN100kU02$N/2)
    datHapN100kU05$PrFixNe  <-  datHapN100kU05$PrFix*(datHapN100kU05$N/2)
    datHapN100kU1$PrFixNe  <-  datHapN100kU1$PrFix*(datHapN100kU1$N/2)
    datHapN100kU02$PrFix[datHapN100kU02$PrFix < (2/(10^8))]  <-  (2/(10^8))
    datHapN100kU05$PrFix[datHapN100kU05$PrFix < (2/(10^8))]  <-  (2/(10^8))
    datHapN100kU1$PrFix[datHapN100kU1$PrFix < (2/(10^8))]  <-  (2/(10^8))
    rdatHapN100kU02  <-  read.csv(file = './data/RECODE/PrFixFig_Haploid_rFixedInv_h0.1_s0.01_N1mil_U02.csv', header=TRUE)
    rdatHapN100kU05  <-  read.csv(file = './data/RECODE/PrFixFig_Haploid_rFixedInv_h0.1_s0.01_N1mil_U05.csv', header=TRUE)
    rdatHapN100kU1  <-  read.csv(file = './data/RECODE/PrFixFig_Haploid_rFixedInv_h0.1_s0.01_N1mil_U1.csv', header=TRUE)
    rdatHapN100kU02  <-  subset(rdatHapN100kU02,rdatHapN100kU02$rFixedInv == 0)
    rdatHapN100kU05  <-  subset(rdatHapN100kU05,rdatHapN100kU05$rFixedInv == 0)
    rdatHapN100kU1   <-  subset(rdatHapN100kU1,rdatHapN100kU1$rFixedInv == 0)
    # Add 0's for Pr(fix | r=0) for table missing table values
    for(i in 1:length(unique(datHapN100kU02$x))) {
        if(!any(rdatHapN100kU02$rInvSize == datHapN100kU02$x[i])) {
            rdatHapN100kU02  <-  rbind(rdatHapN100kU02, rdatHapN100kU02[i-1,])
            rdatHapN100kU02$rInvSize[i]  <-  datHapN100kU02$x[i]
            rdatHapN100kU02$Freq[i]  <-  (2/(10^8))
        }
        if(!any(rdatHapN100kU05$rInvSize == datHapN100kU05$x[i])) {
            rdatHapN100kU05  <-  rbind(rdatHapN100kU05, rdatHapN100kU05[i-1,])
            rdatHapN100kU05$rInvSize[i]  <-  datHapN100kU05$x[i]
            rdatHapN100kU05$Freq[i]  <-  (2/(10^8))
        }
        if(!any(rdatHapN100kU1$rInvSize == datHapN100kU1$x[i])) {
            rdatHapN100kU1  <-  rbind(rdatHapN100kU1, rdatHapN100kU1[i-1,])
            rdatHapN100kU1$rInvSize[i]  <-  datHapN100kU1$x[i]
            rdatHapN100kU1$Freq[i]  <-  (2/(10^8))
        }
    }



    # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  colfunc(3)

    # set plot layout
    layout.mat <- matrix(c(1:2), nrow=1, ncol=2, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

## Panel A & B -- Autosomal inversions
## Panel A
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(5,4,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c((1/(2*10^8)), (6/(2*10^5))), log='y', ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80', log='y')
        box()
        # Points
        abline(h=(1/(2*10^5)), lwd=1, lty=2, col=1)
        points(PrFix ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=AdatHapN100kU02)
        points(PrFix ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=AdatHapN100kU05)
        points(PrFix ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=AdatHapN100kU1)
        # axes
        axis(1, las=1)
        axis(2, las=1, at=c((1/(2*10^8)), 
                            (1/(2*10^7)), 
                            (1/(2*10^6)), 
                            (1/(2*10^5)), 
                            (1/(2*10^4))), 
                    labels=c(0, 
                            expression(1/(2%*%10^7)), 
                            expression(1/(2%*%10^6)), 
                            expression(1/(2%*%10^5)), 
                            expression(1/(2%*%10^4))))
        proportionalLabel(0.5, 1.25, 'Autosomal', cex=1.75, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.04, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(-0.0175, 0.16, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(-0.0175, 0.14, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(0.5, 1.06, expression(italic(N)==10^5) , cex=1.25, adj=c(0.5, 0.5), xpd=NA, log='y')
        # Plot labels etc.
        proportionalLabel(-0.475, 0.5, expression(Fixation~Probability), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90, log='y')
        proportionalLabel( 0.5,  -0.25,  expression(paste("Inversion size (", italic(x),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')


## Panel B

# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c((2/10^8), 12/(10^5)), log='y', ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80', log='y')
        box()
        # Benchmark 2/N
        abline(h=(2/10^6), lwd=1, lty=2, col=1)
        # PrFix Points
        points(PrFix ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datHapN100kU02)
        points(PrFix ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datHapN100kU05)
        points(PrFix ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datHapN100kU1)
        # r0 PrFix Points
#        points(Freq ~ rInvSize, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rdatHapN100kU02)
#        points(Freq ~ rInvSize, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rdatHapN100kU05)
#        points(Freq ~ rInvSize, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rdatHapN100kU1)
        # f0*qStar
qStars  <-  qstar(q.0=(2/datHapN100kU02$N[1]), U=datHapN100kU02$U[1], x=datHapN100kU02$x, s=datHapN100kU02$s[1], h=datHapN100kU02$h[1])        
p0  <-  exp(-(datHapN100kU02$U[1] * datHapN100kU02$x)/(datHapN100kU02$s[1]*datHapN100kU02$h[1]))
detExp  <-  p0*qStars
detExp[detExp < 2/(10^8)]  <-  2/(10^8)
points(detExp ~ datHapN100kU02$x, pch=21)
qStars  <-  qstar(q.0=(2/datHapN100kU05$N[1]), U=datHapN100kU05$U[1], x=datHapN100kU05$x, s=datHapN100kU05$s[1], h=datHapN100kU05$h[1])        
p0  <-  exp(-(datHapN100kU05$U[1] * datHapN100kU05$x)/(datHapN100kU05$s[1]*datHapN100kU05$h[1]))
detExp  <-  p0*qStars
detExp[detExp < 2/(10^8)]  <-  2/(10^8)
points(detExp ~ datHapN100kU05$x, pch=22)
qStars  <-  qstar(q.0=(2/datHapN100kU1$N[1]), U=datHapN100kU1$U[1], x=datHapN100kU1$x, s=datHapN100kU1$s[1], h=datHapN100kU1$h[1])        
p0  <-  exp(-(datHapN100kU1$U[1] * datHapN100kU1$x)/(datHapN100kU1$s[1]*datHapN100kU1$h[1]))
detExp  <-  p0*qStars
detExp[detExp < 2/(10^8)]  <-  2/(10^8)
points(detExp ~ datHapN100kU1$x, pch=24)
        # axes
        axis(1, las=1)
        axis(2, las=1, at=c(2/10^8, 
                            2/10^7, 
                            2/10^6, 
                            2/10^5, 
                            2/10^4), 
                    labels=c(0, 
                            expression(2/10^7), 
                            expression(2/10^6), 
                            expression(2/10^5), 
                            expression(2/10^4)))
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.5, 1.25, 'Expand SLR', cex=1.75, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.5, 1.06, expression(italic(N)==10^6) , cex=1.25, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(-0.0175, 0.16, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(-0.0175, 0.14, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Inversion size (", italic(x),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')

        # Legend
        legend(
               x       =  1.02,
               y       =  2/(10^4),
               legend  =  c(
                            expression(paste("All ", italic(d))),
                            expression(paste(italic(f)[0],italic(q)^"*"))),
               pch     =  c(24,24,24),
               col     =  transparentColor(COLS[1], opacity=1),
               pt.bg   =  c(transparentColor(COLS[2], opacity=0.6), 
                            NA),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
        legend(
               x       =  0.81,
               y       =  2/(10^4),
               legend  =  c(NA,NA),
               pch     =  c(22,22),
               col     =  transparentColor(COLS[1], opacity=1),
               pt.bg   =  c(transparentColor(COLS[2], opacity=0.6), 
                            NA),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
        legend(
               x       =  0.745,
               y       =  2/(10^4),
               legend  =  c(NA,NA),
               pch     =  c(21,21),
               col     =  transparentColor(COLS[1], opacity=1),
               pt.bg   =  c(transparentColor(COLS[2], opacity=0.6), 
                            NA),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

}




####################################

deterministicSuppFig_AutoExact  <-  function(h = 0.25, U = 0.02) {


    # set plot layout
    layout.mat <- matrix(c(1:20), nrow=4, ncol=5, byrow=FALSE)
    layout     <- layout(layout.mat,respect=TRUE)

    # Make data for plotting (x = 0.1; 1st column of panels)
    s        <-  0.01
    nTot     <-  10^4 
    u        <-  U/nTot
    qHat     <-  (U/(nTot*h*s))
    x        <-  0.05
    rBar     <-  (U*x)/(h*s)
    r.vals   <-  c(0,1)
    # make data for plotting
    r.0.Dat  <-  makeAutoDeterministicFigSimData_Nei(r = 0, U = U, x = x, h = h, generations = 3500)
    r.1.Dat  <-  makeAutoDeterministicFigSimData_Nei(r = 1, U = U, x = x, h = h, generations = 3500)

    r.0.rel.w.I.t  <-  c()
    r.1.rel.w.I.t  <-  c()
    for(i in 2:length(r.0.Dat$I.t)) {
        r.0.rel.w.I.t[i-1]  <-  r.0.Dat$I.t[i]/r.0.Dat$I.t[i-1]
        r.1.rel.w.I.t[i-1]  <-  r.1.Dat$I.t[i]/r.1.Dat$I.t[i-1]
    }
    r.0.rel.w.I.t  <-  c(r.0.rel.w.I.t,NA)
    r.1.rel.w.I.t  <-  c(r.1.rel.w.I.t,NA)
    r.1.rel.w.I.t  <-  r.1.rel.w.I.t[r.1.Dat$I.t > 5*10^-5]

    # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  rev(colfunc(5))


## Panel A
# Distribution of captured deleterious alleles (r) given inversion size (x)
# with benchmarks for 5th, 10th, 25th, and 40th quantiles, as well as the mean.
    pDens  <-  dpois(x=c(0:10), lambda=rBar)
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,length(pDens)), ylim = c(0, 0.7), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(pDens~c(0:(length(pDens)-1)), lwd=1.5, lty=1, col=transparentColor("#252525", opacity=1))
        points(pDens~c(0:(length(pDens)-1)), lwd=1.5, lty=1, pch=21, col=transparentColor("#252525", opacity=1), bg=transparentColor("#252525", opacity=0.25))
        abline(v=rBar, col=1, lwd=1)
        abline(v=r.vals, lwd=1, lty=2)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.5, 1.1, substitute(italic(x)==xx, list(xx=x)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.04, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(paste(Pr(italic(d)~"|"~italic(x)))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel( 0.5,  -0.2,  expression(italic(d)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(d)=="{0,1}"),
                            expression(bar(italic(d))==frac(italic(U*x),italic(h*s)))),
               lty     =  c(2,1),
               lwd     =  c(1,1),
               col     =  c(2,1),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

## Panel B
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(min(na.omit(r.1.rel.w.I.t)), 1+U*x*1.7), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(r.0.rel.w.I.t ~ seq_along(r.0.rel.w.I.t), lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1))
        lines(r.1.rel.w.I.t ~ seq_along(r.1.rel.w.I.t), lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1))
        abline(h=1, lwd=1, lty=2, col=1)
        # Inversion relative fitness -- Haploid model
#        lines(rel.w.I.t[1:length(r.0.rel.w.I.t)] ~ seq_along(rel.w.I.t)[1:length(r.0.rel.w.I.t)], lwd=2, lty=3, col=transparentColor("tomato", opacity=1), data=r.0.Hap.Dat)
#        lines(rel.w.I.t[1:length(r.1.rel.w.I.t)] ~ seq_along(rel.w.I.t)[1:length(r.1.rel.w.I.t)], lwd=2, lty=3, col=transparentColor("tomato", opacity=1), data=r.1.Hap.Dat)
        # Annotate 1 + Ux approximation
        points((1+U*x) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        # Points depicting (arbitrary) extinction frequency of 0.00001
        #points(r.0.rel.w.I.t[length(r.0.rel.w.I.t)] ~ length(r.0.rel.w.I.t), pch=8, col='tomato')
        points(r.1.rel.w.I.t[length(r.1.rel.w.I.t)] ~ length(r.1.rel.w.I.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        proportionalArrows(0.125, 0.81, 0.06, 0.81, cex=1, adj=c(0.5, 0.5), xpd=NA, log=FALSE, length=0.05)
        proportionalLabel(0.3, 0.81, expression(1+italic(Ux)), cex=1.3, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(paste("Inversion Relative Fitness")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(d)==0),
                            expression(italic(d)==1)),
               lty     =  1,
               lwd     =  1.5,
               col     =  c(transparentColor(COLS[5], opacity=1),
                            transparentColor(COLS[4], opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

## Panel C
# Effect of initial del. mut. load on inversion frequency dynamics 
    qstr  <-  qstar(U=U, q.0=r.0.Dat$I.t[1], x=x, s=s, h=h)
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,(qstr*1.25)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        lines(I.t[1:length(r.0.rel.w.I.t)] ~ seq_along(I.t)[1:length(r.0.rel.w.I.t)], lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1), data=r.0.Dat)
        lines(I.t[1:length(r.1.rel.w.I.t)] ~ seq_along(I.t[1:length(r.1.rel.w.I.t)]), lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.1.Dat)
        # Inversion frequencies -- Haploid model
#        lines(I.t[1:length(r.0.rel.w.I.t)] ~ seq_along(I.t)[1:length(r.0.rel.w.I.t)], lwd=2, lty=3, col=transparentColor("tomato", opacity=1), data=r.0.Hap.Dat)
#        lines(I.t[1:length(r.1.rel.w.I.t)] ~ seq_along(I.t)[1:length(r.1.rel.w.I.t)], lwd=2, lty=3, col=transparentColor("tomato", opacity=1), data=r.1.Hap.Dat)
        # Points
#        points(I.t[length(r.0.rel.w.I.t)] ~ length(r.0.rel.w.I.t), pch=8, col='tomato', data=r.0.Dat)
        points(I.t[length(r.1.rel.w.I.t)] ~ seq_along(I.t)[length(r.1.rel.w.I.t)], pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.1.Dat)
        # Benchmark qStar
        abline(h=qstr, lty=2)        
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.9, 0.85, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(paste("Inversion Frequency (", italic(Y[I]), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
#        legend(
#               x       =  usr[2]*0.7,
#               y       =  usr[4],
#               legend  =  c(expression("Haploid Model")),
#               lty     =  3,
#               lwd     =  2,
#               col     =  "tomato",
#               cex     =  1,
#               xjust   =  1,
#               yjust   =  1,
#               bty     =  'n',
#               border  =  NA
#               )


    # reset Colors
     COLS  <-  c('red', '#252525')

## Panel D
# deleterious allele frequencies for r = 5
     # Make plot
    delMax  <-  r.1.Dat$qt.W[1]*1.25
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines((qt.W/(1-I.t))[1:length(r.1.rel.w.I.t)] ~ seq_along(qt.W[1:length(r.1.rel.w.I.t)]), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.0.Dat)
        lines((qt.D/(1-I.t))[1:length(r.1.rel.w.I.t)] ~ seq_along(qt.D[1:length(r.1.rel.w.I.t)]), lwd=1.5, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.0.Dat)
        lines((qt.I.W/I.t)[1:length(r.1.rel.w.I.t)] ~ seq_along(qt.I.W[1:length(r.1.rel.w.I.t)]), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        # Points indicating inversion extinctions
        points((qt.W/(1-I.t))[length(r.1.rel.w.I.t)] ~ length(r.1.rel.w.I.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.0.Dat)
        points((qt.D/(1-I.t))[length(r.1.rel.w.I.t)] ~ length(r.1.rel.w.I.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.0.Dat)
        points((qt.I.W/I.t)[length(r.1.rel.w.I.t)] ~ length(r.1.rel.w.I.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.0.Dat)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(italic(d)==0), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(paste("Deleterious Allele Frequency")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4]*0.35,
               legend = c(expression(italic(Q)[italic(N)]^italic(W)), 
                          expression(italic(Q)[italic(I)]^italic(W)), 
                          expression(italic(Q)[italic(N)]^italic(D))),
               lty     =  c(1,1,3),
               lwd     =  1.5,
               col     =  c(transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


## 2nd Column of Panels: x = 0.1
    # Set inversion size, r values
    x        <-  0.1

    # Select r values
    rBar  <-  (U*x)/(h*s)
    r.vals  <-  unique(round(seq(0, rBar, length=3)))

    # make data for plotting
    r.Dat.list      <-  list()
    r.relW.list     <-  list()
    for(i in 1:length(r.vals)) {
        r.Dat.list[[i]] <-  makeAutoDeterministicFigSimData_Nei(r = r.vals[i], U = U, x = x, h = h, generations = 3500)
        r.relW.list[[i]]  <-  c(NA)
        for(j in 2:length(r.Dat.list[[i]]$I.t)) {
                r.relW.list[[i]][j-1]  <-  r.Dat.list[[i]]$I.t[j]/r.Dat.list[[i]]$I.t[j-1]
            }
        r.relW.list[[i]]  <-  c(r.relW.list[[i]],NA)
        r.relW.list[[i]]  <-  r.relW.list[[i]][r.Dat.list[[i]]$I.t > 5*10^-5]
        r.relW.list[[i]]  <-  r.relW.list[[i]][r.Dat.list[[i]]$I.t < 0.9996]
        r.relW.list[[i]]  <-  na.omit(r.relW.list[[i]])
    } 
    

   # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  rev(colfunc(5))

## Panel E
# Distribution of captured deleterious alleles (r) given inversion size (x)
# with benchmarks for 5th, 10th, 25th, and 40th quantiles, as well as the mean.
    pDens  <-  dpois(x=c(0:10), lambda=rBar)
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,length(pDens)), ylim = c(0, 0.7), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(pDens~c(0:(length(pDens)-1)), lwd=1.5, lty=1, col=transparentColor("#252525", opacity=1))
        points(pDens~c(0:(length(pDens)-1)), lwd=1.5, lty=1, pch=21, col=transparentColor("#252525", opacity=1), bg=transparentColor("#252525", opacity=0.25))
        abline(v=rBar, col=1, lwd=1)
        abline(v=r.vals, col=2, lwd=1, lty=2)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.5, 1.1, substitute(italic(x)==xx, list(xx=x)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.04, 1.05, 'E', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.2,  expression(italic(d)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            substitute(italic(d)=={rs}, list(rs=deparse(r.vals))),
#                            substitute(italic(d)==xx, list(xx=paste(r.vals, sep=", "))),
                            expression(bar(italic(d))==frac(italic(U*x),italic(h*s)))),
               lty     =  c(2,1),
               lwd     =  c(1,1),
               col     =  c(2,1),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


## Panel F
# Effect of initial del. mut. load on inversion relative fitness 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.Dat.list[[1]])), ylim = c(min(na.omit(r.relW.list[[length(r.relW.list)]])), 1+U*x*1.7), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        for(i in 1:length(r.Dat.list)){
            lines(r.relW.list[[i]] ~ seq_along(r.relW.list[[i]]), lwd=1.5, lty=1, col=transparentColor(COLS[6-i], opacity=1))
        }
        abline(h=1, lwd=1, lty=2, col=1)
        # Annotate 1 + Ux approximation
        points((1+U*x) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        # Points depicting (arbitrary) extinction frequency of 0.00001
#        points(r.0.rel.w.I.t[length(r.0.rel.w.I.t)] ~ length(r.0.rel.w.I.t), pch=8, col='tomato')
        for(i in 2:length(r.relW.list)) {
            points(r.relW.list[[i]][length(r.relW.list[[i]])] ~ length(r.relW.list[[i]]), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[6-i], opacity=0.7))    
        }
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'F', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            paste("r = ", r.vals[1]),
                            paste("r = ", r.vals[2]),
                            paste("r = ", r.vals[3])),
               lty     =  1,
               lwd     =  1.5,
               col     =  c(transparentColor(COLS[5], opacity=1),
                            transparentColor(COLS[4], opacity=1),
                            transparentColor(COLS[3], opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


## Panel G
# Effect of initial del. mut. load on inversion frequency dynamics 
    qstr  <-  qstar(U=U, q.0=r.0.Dat$I.t[1], x=x, s=s, h=h)
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,(qstr*1.25)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        for(i in 1:length(r.Dat.list)){
            lines(r.Dat.list[[i]]$I.t[1:length(r.relW.list[[i]])] ~ seq_along(r.relW.list[[i]]), lwd=1.5, lty=1, col=transparentColor(COLS[6-i], opacity=1))
        }
        # Points
        for(i in 2:length(r.Dat.list)){
            points(r.Dat.list[[i]]$I.t[length(r.relW.list[[i]])] ~ length(r.relW.list[[i]]), pch=21, col=transparentColor(COLS[6-i], opacity=1), bg=transparentColor(COLS[6-i], opacity=0.7))
        }
        # Benchmark qStar
        abline(h=qstr, lty=2)        
        proportionalLabel(0.9, 0.85, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'G', cex=1.2, adj=c(0.5, 0.5), xpd=NA)


    # reset Colors
     COLS  <-  c('red', '#252525')

## Panel H
# deleterious allele frequencies for r = 5
     # Make plot
    delMax  <-  r.Dat.list[[2]]$qt.W[1]*1.25
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.Dat.list[[2]])), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines((qt.W/(1-I.t))[1:length(r.relW.list[[2]])] ~ seq_along(qt.W[1:length(r.relW.list[[2]])]), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.Dat.list[[2]])
        lines((qt.D/(1-I.t))[1:length(r.relW.list[[2]])] ~ seq_along(qt.D[1:length(r.relW.list[[2]])]), lwd=1.5, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.Dat.list[[2]])
        lines((qt.I.W/I.t)[1:length(r.relW.list[[2]])] ~ seq_along(qt.I.W[1:length(r.relW.list[[2]])]), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.Dat.list[[2]])
        # Points indicating inversion extinctions
        points((qt.W/(1-I.t))[length(r.relW.list[[2]])] ~ length(r.relW.list[[2]]), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.Dat.list[[2]])
        points((qt.D/(1-I.t))[length(r.relW.list[[2]])] ~ length(r.relW.list[[2]]), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.Dat.list[[2]])
        points((qt.I.W/I.t)[length(r.relW.list[[2]])] ~ length(r.relW.list[[2]]), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.Dat.list[[2]])
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'H', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1, substitute(italic(d)==rs, list(rs=r.vals[2])), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        




## 3rd Column of Panels: x = 0.2
    # Set inversion size, r values
    x        <-  0.2

    # Select r values
    rBar    <-  round((U*x)/(h*s), digits=2)
    r.vals  <-  unique(round(seq(0, rBar, length=4)))

    # make data for plotting
    r.Dat.list      <-  list()
    r.relW.list     <-  list()
    for(i in 1:length(r.vals)) {
        r.Dat.list[[i]] <-  makeAutoDeterministicFigSimData_Nei(r = r.vals[i], U = U, x = x, h = h, generations = 3500)
        r.relW.list[[i]]  <-  c(NA)
        for(j in 2:length(r.Dat.list[[i]]$I.t)) {
                r.relW.list[[i]][j-1]  <-  r.Dat.list[[i]]$I.t[j]/r.Dat.list[[i]]$I.t[j-1]
            }
        r.relW.list[[i]]  <-  c(r.relW.list[[i]],NA)
        r.relW.list[[i]]  <-  r.relW.list[[i]][r.Dat.list[[i]]$I.t > 5*10^-5]
        r.relW.list[[i]]  <-  r.relW.list[[i]][r.Dat.list[[i]]$I.t < 0.9996]
        r.relW.list[[i]]  <-  na.omit(r.relW.list[[i]])
    } 
 
   # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  rev(colfunc(5))

## Panel I
# Distribution of captured deleterious alleles (r) given inversion size (x)
# with benchmarks for 5th, 10th, 25th, and 40th quantiles, as well as the mean.
    pDens  <-  dpois(x=c(0:max(10,(2.1*rBar))), lambda=rBar)
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,length(pDens)), ylim = c(0, 0.7), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(pDens~c(0:(length(pDens)-1)), lwd=1.5, lty=1, col=transparentColor("#252525", opacity=1))
        points(pDens~c(0:(length(pDens)-1)), lwd=1.5, lty=1, pch=21, col=transparentColor("#252525", opacity=1), bg=transparentColor("#252525", opacity=0.25))
        abline(v=rBar, lwd=1)
        abline(v=r.vals, lwd=1, lty=2)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.5, 1.1, substitute(italic(x)==xx, list(xx=x)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.04, 1.05, 'I', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.2,  expression(italic(d)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            substitute(italic(d)==xx, list(xx=deparse(r.vals))),
                            expression(bar(italic(d))==frac(italic(U*x),italic(h*s)))),
               lty     =  c(2,1),
               lwd     =  c(1,1),
               col     =  c(2,1),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


## Panel J
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.Dat.list[[1]])), ylim = c(min(na.omit(r.relW.list[[length(r.relW.list)]])), 1+U*x*1.7), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        for(i in 1:length(r.Dat.list)){
            lines(r.relW.list[[i]] ~ seq_along(r.relW.list[[i]]), lwd=1.5, lty=1, col=transparentColor(COLS[6-i], opacity=1))
        }
        abline(h=1, lwd=1, lty=2, col=1)
        # Annotate 1 + Ux approximation
        points((1+U*x) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        # Points depicting (arbitrary) extinction frequency of 0.00001
#        points(r.0.rel.w.I.t[length(r.0.rel.w.I.t)] ~ length(r.0.rel.w.I.t), pch=8, col='tomato')
        for(i in 2:length(r.relW.list)) {
            points(r.relW.list[[i]][length(r.relW.list[[i]])] ~ length(r.relW.list[[i]]), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[6-i], opacity=0.7))    
        }
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'J', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            paste("r = ", r.vals[1]),
                            paste("r = ", r.vals[2]),
                            paste("r = ", r.vals[3]),
                            paste("r = ", r.vals[4])),
               lty     =  1,
               lwd     =  1.5,
               col     =  c(transparentColor(COLS[5], opacity=1),
                            transparentColor(COLS[4], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[2], opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )



## Panel K
# Effect of initial del. mut. load on inversion frequency dynamics 
    qstr  <-  qstar(U=U, q.0=r.Dat.list[[1]]$I.t[1], x=x, s=s, h=h)
    # make plot
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.Dat.list[[1]])), ylim = c(0,(qstr*1.25)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        for(i in 1:length(r.Dat.list)){
            lines(r.Dat.list[[i]]$I.t[1:length(r.relW.list[[i]])] ~ seq_along(r.relW.list[[i]]), lwd=1.5, lty=1, col=transparentColor(COLS[6-i], opacity=1))
        }
        # Points
        for(i in 2:length(r.Dat.list)){
            points(r.Dat.list[[i]]$I.t[length(r.relW.list[[i]])] ~ length(r.relW.list[[i]]), pch=21, col=transparentColor(COLS[6-i], opacity=1), bg=transparentColor(COLS[6-i], opacity=0.7))
        }
        # Benchmark qStar
        abline(h=qstr, lty=2)        
        proportionalLabel(0.9, 0.85, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'K', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # reset Colors
     COLS  <-  c('red', '#252525')

## Panel L
# deleterious allele frequencies for r = 1
     # Make plot
    delMax  <-  r.Dat.list[[2]]$qt.W[1]*1.25
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.Dat.list[[2]])), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines((qt.W/(1-I.t))[1:length(r.relW.list[[2]])] ~ seq_along(qt.W[1:length(r.relW.list[[2]])]), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.Dat.list[[2]])
        lines((qt.D/(1-I.t))[1:length(r.relW.list[[2]])] ~ seq_along(qt.D[1:length(r.relW.list[[2]])]), lwd=1.5, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.Dat.list[[2]])
        lines((qt.I.W/I.t)[1:length(r.relW.list[[2]])] ~ seq_along(qt.I.W[1:length(r.relW.list[[2]])]), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.Dat.list[[2]])
        # Points indicating inversion extinctions
        points((qt.W/(1-I.t))[length(r.relW.list[[2]])] ~ length(r.relW.list[[2]]), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.Dat.list[[2]])
        points((qt.D/(1-I.t))[length(r.relW.list[[2]])] ~ length(r.relW.list[[2]]), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.Dat.list[[2]])
        points((qt.I.W/I.t)[length(r.relW.list[[2]])] ~ length(r.relW.list[[2]]), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.Dat.list[[2]])
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'L', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(italic(d)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        



## 4th Column of Panels: x = 0.1
    # Set inversion size, r values
    x        <-  0.5

  # Select r values
    rBar    <-  round((U*x)/(h*s), digits=2)
    r.vals  <-  unique(round(seq(0, rBar, length=5)))

    # make data for plotting
    r.Dat.list      <-  list()
    r.relW.list     <-  list()
    for(i in 1:length(r.vals)) {
        r.Dat.list[[i]] <-  makeAutoDeterministicFigSimData_Nei(r = r.vals[i], U = U, x = x, h = h, generations = 3500)
        r.relW.list[[i]]  <-  c(NA)
        for(j in 2:length(r.Dat.list[[i]]$I.t)) {
                r.relW.list[[i]][j-1]  <-  r.Dat.list[[i]]$I.t[j]/r.Dat.list[[i]]$I.t[j-1]
            }
        r.relW.list[[i]]  <-  c(r.relW.list[[i]],NA)
        r.relW.list[[i]]  <-  r.relW.list[[i]][r.Dat.list[[i]]$I.t > 5*10^-5]
        r.relW.list[[i]]  <-  r.relW.list[[i]][r.Dat.list[[i]]$I.t < 0.9996]
        r.relW.list[[i]]  <-  na.omit(r.relW.list[[i]])
    } 

   # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  rev(colfunc(5))

## Panel M
# Distribution of captured deleterious alleles (r) given inversion size (x)
# with benchmarks for 5th, 10th, 25th, and 40th quantiles, as well as the mean.
    pDens  <-  dpois(x=c(0:max(20, 2.1*rBar)), lambda=rBar)
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,length(pDens)), ylim = c(0, 0.4), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(pDens~c(0:(length(pDens)-1)), lwd=1.5, lty=1, col=transparentColor("#252525", opacity=1))
        points(pDens~c(0:(length(pDens)-1)), lwd=1.5, lty=1, pch=21, col=transparentColor("#252525", opacity=1), bg=transparentColor("#252525", opacity=0.25))
        abline(v=rBar, col=1, lwd=1)
        abline(v=r.vals, col=2, lwd=1, lty=2)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.5, 1.1, substitute(italic(x)==xx, list(xx=x)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.04, 1.05, 'M', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.2,  expression(italic(d)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            substitute(italic(d)==xx, list(xx=deparse(r.vals))),
                            expression(bar(italic(d))==frac(italic(U*x),italic(h*s)))),
               lty     =  c(2,1),
               lwd     =  c(1,1),
               col     =  c(2,1),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


## Panel N
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.Dat.list[[1]])), ylim = c(min(na.omit(r.relW.list[[length(r.relW.list)]])), 1+U*x*1.7), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        for(i in 1:length(r.Dat.list)){
            lines(r.relW.list[[i]] ~ seq_along(r.relW.list[[i]]), lwd=1.5, lty=1, col=transparentColor(COLS[6-i], opacity=1))
        }
        abline(h=1, lwd=1, lty=2, col=1)
        # Annotate 1 + Ux approximation
        points((1+U*x) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        # Points depicting (arbitrary) extinction frequency of 0.00001
#        points(r.0.rel.w.I.t[length(r.0.rel.w.I.t)] ~ length(r.0.rel.w.I.t), pch=8, col='tomato')
        for(i in 2:length(r.relW.list)) {
            points(r.relW.list[[i]][length(r.relW.list[[i]])] ~ length(r.relW.list[[i]]), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[6-i], opacity=0.7))    
        }
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'N', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            paste("r = ", r.vals[1]),
                            paste("r = ", r.vals[2]),
                            paste("r = ", r.vals[3]),
                            paste("r = ", r.vals[4])),
               lty     =  1,
               lwd     =  1.5,
               col     =  c(transparentColor(COLS[5], opacity=1),
                            transparentColor(COLS[4], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[2], opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


## Panel O
# Effect of initial del. mut. load on inversion frequency dynamics 
    qstr  <-  qstar(U=U, q.0=r.Dat.list[[1]]$I.t[1], x=x, s=s, h=h)
    # make plot
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.Dat.list[[1]])), ylim = c(0,(qstr*1.25)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        for(i in 1:length(r.Dat.list)){
            lines(r.Dat.list[[i]]$I.t[1:length(r.relW.list[[i]])] ~ seq_along(r.relW.list[[i]]), lwd=1.5, lty=1, col=transparentColor(COLS[6-i], opacity=1))
        }
        # Points
        for(i in 2:length(r.Dat.list)){
            points(r.Dat.list[[i]]$I.t[length(r.relW.list[[i]])] ~ length(r.relW.list[[i]]), pch=21, col=transparentColor(COLS[6-i], opacity=1), bg=transparentColor(COLS[6-i], opacity=0.7))
        }
        # Benchmark qStar
        abline(h=qstr, lty=2)        
        proportionalLabel(0.9, 0.85, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'O', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # reset Colors
     COLS  <-  c('red','#252525')

## Panel P
# deleterious allele frequencies for r = 1
     # Make plot
    delMax  <-  r.Dat.list[[2]]$qt.W[1]*1.25
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.Dat.list[[2]])), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines((qt.W/(1-I.t))[1:length(r.relW.list[[2]])] ~ seq_along(qt.W[1:length(r.relW.list[[2]])]), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.Dat.list[[2]])
        lines((qt.D/(1-I.t))[1:length(r.relW.list[[2]])] ~ seq_along(qt.D[1:length(r.relW.list[[2]])]), lwd=1.5, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.Dat.list[[2]])
        lines((qt.I.W/I.t)[1:length(r.relW.list[[2]])] ~ seq_along(qt.I.W[1:length(r.relW.list[[2]])]), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.Dat.list[[2]])
        # Points indicating inversion extinctions
        points((qt.W/(1-I.t))[length(r.relW.list[[2]])] ~ length(r.relW.list[[2]]), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.Dat.list[[2]])
        points((qt.D/(1-I.t))[length(r.relW.list[[2]])] ~ length(r.relW.list[[2]]), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.Dat.list[[2]])
        points((qt.I.W/I.t)[length(r.relW.list[[2]])] ~ length(r.relW.list[[2]]), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.Dat.list[[2]])
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'P', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(italic(d)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        




## 5th Column of Panels: x = 0.8
    # Set inversion size, r values
    x    <-  0.8

  # Select r values
    rBar    <-  round((U*x)/(h*s), digits=2)
    r.vals  <-  unique(round(seq(0, rBar, length=5)))

    # make data for plotting
    r.Dat.list      <-  list()
    r.relW.list     <-  list()
    for(i in 1:length(r.vals)) {
        r.Dat.list[[i]] <-  makeAutoDeterministicFigSimData_Nei(r = r.vals[i], U = U, x = x, h = h, generations = 3500)
        r.relW.list[[i]]  <-  c(NA)
        for(j in 2:length(r.Dat.list[[i]]$I.t)) {
                r.relW.list[[i]][j-1]  <-  r.Dat.list[[i]]$I.t[j]/r.Dat.list[[i]]$I.t[j-1]
            }
        r.relW.list[[i]]  <-  c(r.relW.list[[i]],NA)
        r.relW.list[[i]]  <-  r.relW.list[[i]][r.Dat.list[[i]]$I.t > 5*10^-5]
        r.relW.list[[i]]  <-  r.relW.list[[i]][r.Dat.list[[i]]$I.t < 0.9996]
        r.relW.list[[i]]  <-  na.omit(r.relW.list[[i]])
    } 

   # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  rev(colfunc(5))

## Panel Q
# Distribution of captured deleterious alleles (r) given inversion size (x)
# with benchmarks for 5th, 10th, 25th, and 40th quantiles, as well as the mean.
    qHat  <-  r.0.Dat$qt.Y.W[1]
#    rBar  <-  round(nTot*x*qHat, digits=2)
    pDens  <-  dpois(x=c(0:max(20, 2.1*rBar)), lambda=U*x/(h*s))
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,length(pDens)), ylim = c(0, 0.7), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(pDens~c(0:(length(pDens)-1)), lwd=1.5, lty=1, col=transparentColor("#252525", opacity=1))
        points(pDens~c(0:(length(pDens)-1)), lwd=1.5, lty=1, pch=21, col=transparentColor("#252525", opacity=1), bg=transparentColor("#252525", opacity=0.25))
        abline(v=rBar, col=1, lwd=1)
        abline(v=r.vals, col=2, lwd=1, lty=2)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.5, 1.1, substitute(italic(x)==xx, list(xx=x)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.04, 1.05, 'Q', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.2,  expression(italic(d)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            substitute(italic(d)==xx, list(xx=deparse(r.vals))),
                            expression(bar(italic(d))==frac(italic(U*x),italic(h*s)))),
               lty     =  c(2,1),
               lwd     =  c(1,1),
               col     =  c(2,1),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


## Panel R
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.Dat.list[[1]])), ylim = c(min(na.omit(r.relW.list[[length(r.relW.list)]])), 1+U*x*1.7), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        for(i in 1:length(r.Dat.list)){
            lines(r.relW.list[[i]] ~ seq_along(r.relW.list[[i]]), lwd=1.5, lty=1, col=transparentColor(COLS[6-i], opacity=1))
        }
        abline(h=1, lwd=1, lty=2, col=1)
        # Annotate 1 + Ux approximation
        points((1+U*x) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        # Points depicting (arbitrary) extinction frequency of 0.00001
#        points(r.0.rel.w.I.t[length(r.0.rel.w.I.t)] ~ length(r.0.rel.w.I.t), pch=8, col='tomato')
        for(i in 2:length(r.relW.list)) {
            points(r.relW.list[[i]][length(r.relW.list[[i]])] ~ length(r.relW.list[[i]]), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[6-i], opacity=0.7))    
        }
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'R', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            paste("r = ", r.vals[1]),
                            paste("r = ", r.vals[2]),
                            paste("r = ", r.vals[3]),
                            paste("r = ", r.vals[4])),
               lty     =  1,
               lwd     =  1.5,
               col     =  c(transparentColor(COLS[5], opacity=1),
                            transparentColor(COLS[4], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[2], opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )




## Panel S
# Effect of initial del. mut. load on inversion frequency dynamics 
    qstr  <-  qstar(U=U, q.0=r.Dat.list[[1]]$I.t[1], x=x, s=s, h=h)
    # make plot
    plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.Dat.list[[1]])), ylim = c(0,(qstr*1.25)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        for(i in 1:length(r.Dat.list)){
            lines(r.Dat.list[[i]]$I.t[1:length(r.relW.list[[i]])] ~ seq_along(r.relW.list[[i]]), lwd=1.5, lty=1, col=transparentColor(COLS[6-i], opacity=1))
        }
        # Points
        for(i in 2:length(r.Dat.list)){
            points(r.Dat.list[[i]]$I.t[length(r.relW.list[[i]])] ~ length(r.relW.list[[i]]), pch=21, col=transparentColor(COLS[6-i], opacity=1), bg=transparentColor(COLS[6-i], opacity=0.7))
        }
        # Benchmark qStar
        abline(h=qstr, lty=2)        
        proportionalLabel(0.9, 0.85, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'S', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # reset Colors
     COLS  <-  c('red', '#252525')


## Panel T
# deleterious allele frequencies for r = 1
     # Make plot
    delMax  <-  r.Dat.list[[2]]$qt.W[1]*1.25
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.Dat.list[[2]])), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines((qt.W/(1-I.t))[1:length(r.relW.list[[2]])] ~ seq_along(qt.W[1:length(r.relW.list[[2]])]), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.Dat.list[[2]])
        lines((qt.D/(1-I.t))[1:length(r.relW.list[[2]])] ~ seq_along(qt.D[1:length(r.relW.list[[2]])]), lwd=1.5, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.Dat.list[[2]])
        lines((qt.I.W/I.t)[1:length(r.relW.list[[2]])] ~ seq_along(qt.I.W[1:length(r.relW.list[[2]])]), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.Dat.list[[2]])
        # Points indicating inversion extinctions
        points((qt.W/(1-I.t))[length(r.relW.list[[2]])] ~ length(r.relW.list[[2]]), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.Dat.list[[2]])
        points((qt.D/(1-I.t))[length(r.relW.list[[2]])] ~ length(r.relW.list[[2]]), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.Dat.list[[2]])
        points((qt.I.W/I.t)[length(r.relW.list[[2]])] ~ length(r.relW.list[[2]]), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.Dat.list[[2]])
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'T', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(italic(d)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        

}










####################################
####################################
####################################
#' WF-simulation allele & inversion
#' frequency dynamics figure
SLR_WF_DynamicsFig  <-  function(invSize) {

   # Import data for plotting
    dat  <-  read.csv(file = './data/RECODE/SLR-Nei-WF-Dynamics_h0.25_s0.01_U0.02_N500k.csv', header=TRUE)
#    head(dat)

    dat$invRelFit  <-  invRelFit2(n = (10^4*dat$x - dat$r), 
                                 r = dat$r, 
                                 s = dat$s, 
                                 h = dat$h, 
                                 q.YI.W = dat$qt.YI.W, 
                                 q.Y.W = dat$qt.Y.W, 
                                 q.Y.D = dat$qt.Y.D, 
                                 q.Xf.W = dat$qt.Xf.W, 
                                 q.Xf.D = dat$qt.Xf.D)


   # Colors
    COLS  <- c("#252525", colorRampPalette(wes_palette("Zissou1"))(4))

    # set plot layout
    layout.mat <- matrix(c(1:3), nrow=3, ncol=1, byrow=FALSE)
    layout     <- layout(layout.mat,respect=TRUE)

    # subset data
    subDat  <-  subset(dat, dat$x == invSize)


## Panel A
# Overview of inversion frequency and relative fitness dynamics
    # make plot
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(subDat)), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Inversion frequency
        lines(YI.t ~ seq_along(YI.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=subDat)
        axis(1, las=1)
        axis(4, las=1, cex=1.3)
        # Inversion relative fitness
        par(new=TRUE)
        plot(invRelFit ~ seq_along(invRelFit), type='l', ylim=c((1-(max(invRelFit) - 1)), (1 + subDat$U[1]*subDat$x[1]*1.25)), axes=FALSE, lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=subDat)
        # Benchmark for relative fitness of 1
        abline(h=1, lty=2, lwd=1.5, col=COLS[1])
        # Annotate 1 + Ux approximation
        points((1 + subDat$U[1]*subDat$x[1]) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        proportionalArrows(0.125, 0.86, 0.06, 0.86, cex=1, adj=c(0.5, 0.5), xpd=NA, log=FALSE, length=0.05)
        proportionalLabel(0.25, 0.86, expression(1+italic(Ux)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.5, 1.15, expression(paste(italic(Ux/(hs)==1.6))), cex=2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(1.3,  0.5,   expression(paste("Inversion Frequency (", italic(Y[I]), ")")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=270)
        proportionalLabel(-0.3,  0.5,  expression(paste("Inversion Relative Fitness")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)
        # Legend


    # Deterministic expectations
    t             <-  seq_along(subDat$YI.t)
    YI.gensol     <-  YI.genSol(U=subDat$U[1], YI.init=subDat$YI.t[1], x=subDat$x[1], s=subDat$s[1], h=subDat$h[1], t=t)
    qstr          <-  qstar(U=subDat$U[1], q.0=subDat$YI.t[1], x=subDat$x[1], s=subDat$s[1], h=subDat$h[1])
    subDat2       <-  subDat[(qstr-YI.gensol)/qstr > 0.01,]
    r.0.Dat       <-  makeDeterministicFigSimData(U = subDat2$U[1], s = subDat2$s[1], r = subDat2$r[1], x = subDat2$x[1], h = subDat2$h[1], generations = nrow(subDat2))
    HapInvRelFit  <-  w.YI.X(n=(10^4*subDat2$x[1]), r=subDat2$r[1], t=c(1:nrow(subDat2)), u=subDat2$U[1]/10^4, s=subDat2$s[1], h=subDat2$h[1], qHat=subDat2$qt.Xf.W[1])/ w.Y.X(n=(10^4*subDat2$x[1]), u=subDat2$U[1]/10^4) 
    relGeomS      <-  round((exp(mean(log(subDat2$invRelFit)))-1)/(exp(mean(log(HapInvRelFit)))-1), digits=3)


## Panel B
# Zoom in on initial phase where inversion relative fitness is greater than 1
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(subDat2)), ylim = c(0,max(subDat2$YI.t)*1.15), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        axis(1, las=1)
        axis(4, las=1, cex=1.3)
        # Inversion Frequency
        lines(YI.t ~ seq_along(YI.t), type='l', lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=subDat2)
        lines(YI.gensol ~ seq_along(YI.gensol), type='l', lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=subDat2)
        # Benchmark qstar
        abline(h=qstr, lty=2, lwd=1.5, col='tomato')
        # Inversion relative fitness
        par(new=TRUE)
        plot(invRelFit ~ seq_along(invRelFit), type='l', ylim=c((1-(max(invRelFit) - 1)), (1 + subDat$U[1]*subDat$x[1]*1.25)), axes=FALSE, lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=subDat2)
        lines(rel.w.YI.t ~ seq_along(rel.w.YI.t), lwd=2, lty=3, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
#        lines((r.0.Dat$rel.w.YI.t[1:nrow(subDat2)] - subDat2$invRelFit) ~ seq_along(r.0.Dat$rel.w.YI.t[1:nrow(subDat2)]), lwd=2, lty=3, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        lines(HapInvRelFit ~ seq_along(HapInvRelFit), lty=3, col='tomato')
        # Benchmark for relative fitness of 1
        abline(h=1, lty=2, lwd=1.5, col=COLS[1])
        # axes
        axis(2, las=1, cex=1.3)
        # Annotate 1 + Ux approximation
        points((1 + subDat$U[1]*subDat$x[1]) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(1.3,  0.5,   expression(paste("Inversion Frequency (", italic(Y[I]), ")")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=270)
        proportionalLabel(-0.3,  0.5,  expression(paste("Inversion Relative Fitness")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.15,  0.95,  expression(paste(italic(tilde(s)[WF])/italic(tilde(s)[Det.])," = ")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(0.35, 0.95, substitute(xx, list(xx=relGeomS)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        usr  <-  par('usr')
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(paste("Inv. Rel. Fit.")),
                            expression(paste("Determ. Inv. Rel. Fit.")),
                            expression(italic(I[t]))),
               lty     =  c(1,3,1),
               lwd     =  2,
               col     =  c(transparentColor(COLS[2], opacity=1),
                            transparentColor('tomato', opacity=1),
                            transparentColor("#252525", opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


## Panel C
# Deleterious mutation frequency dynamics
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(subDat2)), ylim = c(min(subDat2$qt.Xf.W),max(subDat2$qt.Xf.W)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        axis(1, las=1)
        axis(2, las=1, cex=1.3)
        # Frequency of del. mutations on Xf
        lines(qt.Xf.W ~ seq_along(qt.Xf.W), type='l', lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=subDat2)
        # Benchmark initial equilibrium frequency
        abline(h=subDat2$qt.Xf.W[1], lty=2, lwd=1.5, col=COLS[1])
        # Frequency of del. mutations on YI
        par(new=TRUE)
        plot(qt.YI.W ~ seq_along(qt.YI.W), type='l', ylim = c(0, max(subDat2$qt.Xf.W)), lwd=2, lty=1, col=transparentColor('red', opacity=1), data=subDat2, axes=FALSE)
        # axes
        axis(4, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5, expression(paste(italic(q[X[m]]))), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel(1.3,  0.5,   expression(paste(italic(q[Y[I]]))), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=270)        
        proportionalLabel( 0.5,  -0.15,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        usr  <-  par('usr')
        legend(
               x       =  usr[2]*0.25,
               y       =  usr[4],
               legend  =  c(
                            expression(italic(q[Y[I]])),
                            expression(italic(q[X[f]]))),
               lty     =  1,
               lwd     =  2,
               col     =  c(transparentColor("red", opacity=1),
                            transparentColor("#252525", opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

}



##########################
#' WF-simulation allele & inversion
#' frequency dynamics Supplementary figure
SLR_WF_DynamicsSuppFig  <-  function() {

   # Import data for plotting
    dat  <-  read.csv(file = './data/RECODE/SLR-Nei-WF-Dynamics_h0.25_s0.01_U0.02_N100k.csv', header=TRUE)

   # Colors
    COLS  <- c("#252525", colorRampPalette(wes_palette("Zissou1"))(4))

    # set plot layout
    layout.mat <- matrix(c(1:30), nrow=3, ncol=10, byrow=FALSE)
    layout     <- layout(layout.mat,respect=TRUE)

    invSizes  <-  unique(dat$x)
    for(i in 1:length(invSizes)) {
        # subset data
        subDat  <-  subset(dat, dat$x == invSizes[i])

## Panel A
# Overview of inversion frequency and relative fitness dynamics
    # make plot
    if(i == 1) {
        par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(4,5,2,2), bty='o', xaxt='s', yaxt='s')
    }
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(subDat)), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2, )
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Inversion frequency
        lines(YI.t ~ seq_along(YI.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=subDat)
        axis(1, las=1)
        axis(4, las=1, cex=1.3)
        # Inversion relative fitness
        par(new=TRUE)
head(subDat)
subDat$invRelFit  
        plot(invRelFit ~ seq_along(invRelFit), type='l', ylim=c((1-(max(invRelFit) - 1)), (1 + subDat$U[1]*subDat$x[1]*1.25)), axes=FALSE, lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), ylab='', xlab='', data=subDat)
        lines(invRelFit ~ seq_along(invRelFit), lwd=2, lty=2, col=transparentColor(COLS[3], opacity=1), data=subDat)
        # Benchmark for relative fitness of 1
        abline(h=1, lty=2, lwd=1.5, col=COLS[1])
        # Annotate 1 + Ux approximation
        points((1 + subDat$U[1]*subDat$x[1]) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        proportionalArrows(0.125, 0.86, 0.06, 0.86, cex=1, adj=c(0.5, 0.5), xpd=NA, log=FALSE, length=0.05)
        proportionalLabel(0.25, 0.86, expression(1+italic(Ux)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # Annotate Inversion Size
                proportionalLabel(0.5, 1.1, substitute(italic(x)==xx, list(xx=subDat$x[1])), cex=2.5, adj=c(0.5, 0.5), xpd=NA)

        # axes
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
#        proportionalLabel(0.04, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.5, 1.15, expression(paste(italic(Ux/(hs)==1.6))), cex=2, adj=c(0.5, 0.5), xpd=NA)
        if(i == 10){        
                    proportionalLabel(1.3,  0.5,   expression(paste("Inversion Frequency (", italic(Y[I]), ")")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=270)
                }
        if(i == 1){
                    proportionalLabel(-0.3,  0.5,  expression(paste("Inversion Relative Fitness")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)
                }
        # Legend


    # Deterministic expectations
    t             <-  seq_along(subDat$YI.t)
    YI.gensol     <-  YI.genSol(U=subDat$U[1], YI.init=subDat$YI.t[1], x=subDat$x[1], s=subDat$s[1], h=subDat$h[1], t=t)
    qstr          <-  qstar(U=subDat$U[1], q.0=subDat$YI.t[1], x=subDat$x[1], s=subDat$s[1], h=subDat$h[1])
    subDat2       <-  subDat[(qstr-YI.gensol)/qstr > 0.01,]
    r.0.Dat       <-  makeDeterministicFigSimData_SLR_Nei(U = subDat2$U[1], s = subDat2$s[1], r = subDat2$r[1], x = subDat2$x[1], h = subDat2$h[1], generations = nrow(subDat2))
    HapInvRelFit  <-  w.YI.X(n=(10^4*subDat2$x[1]), r=subDat2$r[1], t=c(1:nrow(subDat2)), u=subDat2$U[1]/10^4, s=subDat2$s[1], h=subDat2$h[1], qHat=subDat2$qt.Xf.W[1])/ w.Y.X(n=(10^4*subDat2$x[1]), u=subDat2$U[1]/10^4) 
    tildeSWF      <-  round((exp(mean(log(subDat2$invRelFit)))-1), digits=3)
    tildeSDet     <-  round((exp(mean(log(HapInvRelFit)))-1), digits=3)
    relGeomS      <-  round((exp(mean(log(subDat2$invRelFit)))-1)/(exp(mean(log(HapInvRelFit)))-1), digits=3)


## Panel B
# Zoom in on initial phase where inversion relative fitness is greater than 1
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(subDat2)), ylim = c(0,max(subDat2$YI.t)*1.15), ylab='', xlab='', cex.lab=1.2, yaxt="n")
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        axis(1, las=1)
        axis(4, las=1, cex=1.3)
        # Inversion Frequency
        lines(YI.t ~ seq_along(YI.t), type='l', lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=subDat2)
        lines(YI.gensol ~ seq_along(YI.gensol), type='l', lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=subDat2)
        # Benchmark qstar
        abline(h=qstr, lty=2, lwd=1.5, col='tomato')
        # Inversion relative fitness
        par(new=TRUE)
        plot(invRelFit ~ seq_along(invRelFit), type='l', ylim=c((1-(max(invRelFit) - 1)), (1 + subDat$U[1]*subDat$x[1]*1.25)), axes=FALSE, lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), ylab='', xlab='', data=subDat2)
        lines(rel.w.YI.t ~ seq_along(rel.w.YI.t), lwd=2, lty=3, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
#        lines((r.0.Dat$rel.w.YI.t[1:nrow(subDat2)] - subDat2$invRelFit) ~ seq_along(r.0.Dat$rel.w.YI.t[1:nrow(subDat2)]), lwd=2, lty=3, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        lines(HapInvRelFit ~ seq_along(HapInvRelFit), lty=3, col='tomato')
#        lines(invRelFit2 ~ seq_along(invRelFit2), lwd=2, lty=2, col=transparentColor(COLS[3], opacity=1), data=subDat2)
        # Benchmark for relative fitness of 1
        abline(h=1, lty=2, lwd=1.5, col=COLS[1])
        # axes
        axis(2, las=1, cex=1.3)
        # Annotate 1 + Ux approximation
        points((1 + subDat$U[1]*subDat$x[1]) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        # Plot labels etc.
#        proportionalLabel(0.04, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        if(i == 10){        
                    proportionalLabel(1.3,  0.5,   expression(paste("Inversion Frequency (", italic(Y[I]), ")")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=270)
                }
        if(i == 1){
                    proportionalLabel(-0.3,  0.5,  expression(paste("Inversion Relative Fitness")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)
                }
        proportionalLabel(0.15,  0.95,  expression(paste(italic(tilde(s)[WF])/italic(tilde(s)[Det.])," = ")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(0.35, 0.95, substitute(xx, list(xx=relGeomS)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.7, 0.2,  expression(paste(italic(tilde(s)[WF])," = ")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(0.85, 0.2, substitute(xx, list(xx=tildeSWF)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.7, 0.1,  expression(paste(italic(tilde(s)[Det.])," = ")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(0.85, 0.1, substitute(xx, list(xx=tildeSDet)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        usr  <-  par('usr')
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(paste("Inv. Rel. Fit.")),
                            expression(paste("Determ. Inv. Rel. Fit.")),
                            expression(italic(I[t]))),
               lty     =  c(1,3,1),
               lwd     =  2,
               col     =  c(transparentColor(COLS[2], opacity=1),
                            transparentColor('tomato', opacity=1),
                            transparentColor("#252525", opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


## Panel C
# Deleterious mutation frequency dynamics
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(subDat2)), ylim = c(min(subDat2$qt.Xf.W),max(subDat2$qt.Xf.W)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        axis(1, las=1)
        axis(2, las=1, cex=1.3)
        # Frequency of del. mutations on Xf
        lines(qt.Xf.W ~ seq_along(qt.Xf.W), type='l', lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=subDat2)
        # Benchmark initial equilibrium frequency
        abline(h=subDat2$qt.Xf.W[1], lty=2, lwd=1.5, col=COLS[1])
        # Frequency of del. mutations on YI
        par(new=TRUE)
        plot(qt.YI.W ~ seq_along(qt.YI.W), type='l', ylim = c(0, max(subDat2$qt.Xf.W)), lwd=2, lty=1, col=transparentColor('red', opacity=1), data=subDat2, axes=FALSE, ylab='', xlab='')
        # axes
        axis(4, las=1, cex=1.3)
        # Plot labels etc.
        if(i == 10){        
        proportionalLabel(1.3,  0.5,   expression(paste(italic(q[Y[I]]))), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=270)
                }
        if(i == 1){
            proportionalLabel(-0.3,  0.5, expression(paste(italic(q[X[m]]))), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)
                }    
        proportionalLabel( 0.5,  -0.15,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        usr  <-  par('usr')
        legend(
               x       =  usr[2]*0.25,
               y       =  usr[4],
               legend  =  c(
                            expression(italic(q[Y[I]])),
                            expression(italic(q[X[f]]))),
               lty     =  1,
               lwd     =  2,
               col     =  c(transparentColor("red", opacity=1),
                            transparentColor("#252525", opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

    }

}







####################################
####################################
####################################
#' WF-simulation allele & inversion
#' frequency dynamics figure
AutoExact_WF_DynamicsFig  <-  function(invSize) {

   # Import data for plotting
    dat  <-  read.csv(file = './data/RECODE/AutoExact-WF-Dynamics_h0.25_s0.01_U0.02_N500k.csv', header=TRUE)

    dat$invRelFit  <-  autoInvRelFit(qN = dat$qt.W, 
                                     qI = dat$qt.I.W, 
                                     qD = dat$qt.D,
                                     I  = dat$I.t,
                                     h  = dat$h,
                                     s  = dat$s,
                                     n  = (10^4*dat$x - dat$r),
                                     r  = dat$r)
    dat$invRelFit[is.nan(dat$invRelFit)]  <-  NA

   # Colors
    COLS  <- c("#252525", colorRampPalette(wes_palette("Zissou1"))(4))

    # set plot layout
    layout.mat <- matrix(c(1:3), nrow=3, ncol=1, byrow=FALSE)
    layout     <- layout(layout.mat,respect=TRUE)

    # subset data
    subDat  <-  subset(dat, dat$x == invSize)


## Panel A
# Overview of inversion frequency and relative fitness dynamics
    # make plot
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(subDat)), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Inversion frequency
        lines(I.t ~ seq_along(I.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=subDat)
        axis(1, las=1)
        axis(4, las=1, cex=1.3)
        # Inversion relative fitness
        par(new=TRUE)
        plot(invRelFit ~ seq_along(invRelFit), type='l', ylim=c((1-(max(na.omit(invRelFit)) - 1)), (1 + subDat$U[1]*subDat$x[1]*1.25)), axes=FALSE, lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=subDat)
        # Benchmark for relative fitness of 1
        abline(h=1, lty=2, lwd=1.5, col=COLS[1])
        # Annotate 1 + Ux approximation
        points((1 + subDat$U[1]*subDat$x[1]) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        proportionalArrows(0.125, 0.86, 0.06, 0.86, cex=1, adj=c(0.5, 0.5), xpd=NA, log=FALSE, length=0.05)
        proportionalLabel(0.25, 0.86, expression(1+italic(Ux)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.5, 1.15, expression(paste(italic(Ux/(hs)==1.6))), cex=2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(1.3,  0.5,   expression(paste("Inversion Frequency (", italic(Y[I]), ")")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=270)
        proportionalLabel(-0.3,  0.5,   expression(paste("Inversion Relative Fitness")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)
        # Legend

        # deterministic expectations
        t         <-  seq_along(subDat$I.t)
        I.gensol  <-  I.genSol(U=subDat$U[1], q.0=subDat$I.t[1], x=subDat$x[1], s=subDat$s[1], h=subDat$h[1], t=t)
        qstr      <-  qstar(U=subDat$U[1], q.0=subDat$I.t[1], x=subDat$x[1], s=subDat$s[1], h=subDat$h[1])

        subDat2   <-  subDat[(qstr-I.gensol)/qstr > 0.01,]
        r.0.Dat  <-  makeAutoDeterministicFigSimData_Nei(U = subDat2$U[1], s = subDat2$s[1], r = 0, x = subDat2$x[1], h = subDat2$h[1], generations = nrow(subDat2))
        r.0.rel.w.I.t  <-  c()
        for(i in 2:length(r.0.Dat$I.t)) {
            r.0.rel.w.I.t[i-1]  <-  r.0.Dat$I.t[i]/r.0.Dat$I.t[i-1]
        }
        r.0.rel.w.I.t  <-  c(r.0.rel.w.I.t,NA)
        relGeomS  <-  round((exp(mean(log(subDat2$invRelFit)))-1)/(exp(mean(log(na.omit(r.0.rel.w.I.t))))-1), digits=3)


## Panel B
# Zoom in on initial phase where inversion relative fitness is greater than 1
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(subDat2)), ylim = c(0,max(subDat2$I.t)*1.15), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        axis(1, las=1)
        axis(4, las=1, cex=1.3)
        # Inversion Frequency
        lines(I.t ~ seq_along(I.t), type='l', lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=subDat2)
        lines(I.gensol ~ seq_along(I.gensol), type='l', lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=subDat2)
        # Benchmark qstar
        abline(h=qstr, lty=2, lwd=1.5, col='tomato')
        # Inversion relative fitness
        par(new=TRUE)
        plot(invRelFit ~ seq_along(invRelFit), type='l', ylim=c((1-(max(invRelFit) - 1)), (1 + subDat$U[1]*subDat$x[1]*1.25)), axes=FALSE, lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=subDat2)
        lines(r.0.rel.w.I.t ~ seq_along(r.0.rel.w.I.t), type='l', lwd=2, lty=3, col=transparentColor('tomato', opacity=1))
        # Benchmark for relative fitness of 1
        abline(h=1, lty=2, lwd=1.5, col=COLS[1])
        # axes
        axis(2, las=1, cex=1.3)
        # Annotate 1 + Ux approximation
        points((1 + subDat$U[1]*subDat$x[1]) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(1.3,  0.5,   expression(paste("Inversion Frequency (", italic(Y[I]), ")")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=270)
        proportionalLabel(-0.3,  0.5,   expression(paste("Inversion Relative Fitness")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel(0.15,  0.95,  expression(paste(italic(tilde(s)[WF])/italic(tilde(s)[Det.])," = ")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(0.35, 0.95, substitute(xx, list(xx=relGeomS)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        usr  <-  par('usr')
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(paste("Inv. Rel. Fit.")),
                            expression(paste("Determ. Inv. Rel. Fit.")),
                            expression(italic(I[t]))),
               lty     =  c(1,3,1),
               lwd     =  2,
               col     =  c(transparentColor(COLS[2], opacity=1),
                            transparentColor('tomato', opacity=1),
                            transparentColor("#252525", opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


## Panel C
# Deleterious mutation frequency dynamics
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(subDat2)), ylim = c(min(subDat2$qt.W/(1 - subDat2$I.t)), max(subDat2$qt.W/(1 - subDat2$I.t))), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        axis(1, las=1)
        axis(2, las=1, cex=1.3)
        # Frequency of del. mutations on Xf
        lines(qt.W/(1 - I.t) ~ seq_along(qt.W), type='l', lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=subDat2)
        # Benchmark initial equilibrium frequency
        abline(h=(subDat2$qt.W[1]/(1 - subDat2$I.t[1])), lty=2, lwd=1.5, col=COLS[1])
        # Frequency of del. mutations on YI
        par(new=TRUE)
        plot(qt.I.W/I.t ~ seq_along(qt.I.W), type='l', ylim = c(0, max(subDat2$qt.I.W/subDat2$I.t)), lwd=2, lty=1, col=transparentColor('red', opacity=1), data=subDat2, axes=FALSE)
        # axes
        axis(4, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5, expression(paste(italic(q[N]^W))), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel(1.3,  0.5,   expression(paste(italic(q[I]))), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=270)        
        proportionalLabel( 0.5,  -0.15,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        usr  <-  par('usr')
        legend(
               x       =  usr[2]*0.25,
               y       =  usr[4],
               legend  =  c(
                            expression(paste(italic(q[I]))),
                            expression(paste(italic(q[N]^W)))),
               lty     =  1,
               lwd     =  2,
               col     =  c(transparentColor("red", opacity=1),
                            transparentColor("#252525", opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

}




##########################
#' WF-simulation allele & inversion
#' frequency dynamics Supplementary figure

AutoExact_WF_DynamicsFigSuppFig  <-  function() {

   # Import data for plotting
    dat  <-  read.csv(file = './data/RECODE/AutoExact-WF-Dynamics_h0.25_s0.01_U0.02_N500k.csv', header=TRUE)

    dat$invRelFit  <-  autoInvRelFit(qN = dat$qt.W, 
                                     qI = dat$qt.I.W, 
                                     qD = dat$qt.D,
                                     I  = dat$I.t,
                                     h  = dat$h,
                                     s  = dat$s,
                                     n  = (10^4*dat$x - dat$r),
                                     r  = dat$r)
    dat$invRelFit[is.nan(dat$invRelFit)]  <-  NA

   # Colors
    COLS  <- c("#252525", colorRampPalette(wes_palette("Zissou1"))(4))

    # set plot layout
    layout.mat <- matrix(c(1:30), nrow=3, ncol=10, byrow=FALSE)
    layout     <- layout(layout.mat,respect=TRUE)

    invSizes  <-  unique(dat$x)
    for(i in 1:length(invSizes)) {
        # subset data
        subDat  <-  subset(dat, dat$x == invSizes[i])

## Panel A
# Overview of inversion frequency and relative fitness dynamics
    # make plot
    if(i == 1) {
        par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(4,5,2,2), bty='o', xaxt='s', yaxt='s')
    }
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(subDat)), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Inversion frequency
        lines(I.t ~ seq_along(I.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=subDat)
        axis(1, las=1)
        axis(4, las=1, cex=1.3)
        # Inversion relative fitness
        par(new=TRUE)
        plot(invRelFit ~ seq_along(invRelFit), type='l', ylim=c((1-(max(na.omit(invRelFit)) - 1)), (1 + subDat$U[1]*subDat$x[1]*1.25)), axes=FALSE, ylab="", xlab='', lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=subDat)
        # Benchmark for relative fitness of 1
        abline(h=1, lty=2, lwd=1.5, col=COLS[1])
        # Annotate 1 + Ux approximation
        points((1 + subDat$U[1]*subDat$x[1]) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        proportionalArrows(0.125, 0.86, 0.06, 0.86, cex=1, adj=c(0.5, 0.5), xpd=NA, log=FALSE, length=0.05)
        proportionalLabel(0.25, 0.86, expression(1+italic(Ux)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.5, 1.15, expression(paste(italic(Ux/(hs)==1.6))), cex=2, adj=c(0.5, 0.5), xpd=NA)
         # Annotate Inversion Size
                proportionalLabel(0.5, 1.1, substitute(italic(x)==xx, list(xx=subDat$x[1])), cex=2.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
#        proportionalLabel(0.04, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.5, 1.15, expression(paste(italic(Ux/(hs)==1.6))), cex=2, adj=c(0.5, 0.5), xpd=NA)
        if(i == 10){        
                    proportionalLabel(1.3,  0.5,   expression(paste("Inversion Frequency (", italic(Y[I]), ")")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=270)
                }
        if(i == 1){
                    proportionalLabel(-0.3,  0.5,  expression(paste("Inversion Relative Fitness")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)
                }
        # Legend


        # deterministic expectations
        t         <-  seq_along(subDat$I.t)
        I.gensol  <-  I.genSol(U=subDat$U[1], q.0=subDat$I.t[1], x=subDat$x[1], s=subDat$s[1], h=subDat$h[1], t=t)
        qstr      <-  qstar(U=subDat$U[1], q.0=subDat$I.t[1], x=subDat$x[1], s=subDat$s[1], h=subDat$h[1])

        subDat2   <-  subDat[(qstr-I.gensol)/qstr > 0.01,]
        r.0.Dat  <-  makeAutoDeterministicFigSimData_Nei(U = subDat2$U[1], s = subDat2$s[1], r = 0, x = subDat2$x[1], h = subDat2$h[1], generations = nrow(subDat2))
        r.0.rel.w.I.t  <-  c()
        for(i in 2:length(r.0.Dat$I.t)) {
            r.0.rel.w.I.t[i-1]  <-  r.0.Dat$I.t[i]/r.0.Dat$I.t[i-1]
        }
        r.0.rel.w.I.t  <-  c(r.0.rel.w.I.t,NA)
        tildeSWF  <-  round((exp(mean(log(subDat2$invRelFit)))-1), digits=4)
        tildeSDet <-  round((exp(mean(log(na.omit(r.0.rel.w.I.t))))-1), digits=4)
        relGeomS  <-  round((exp(mean(log(subDat2$invRelFit)))-1)/(exp(mean(log(na.omit(r.0.rel.w.I.t))))-1), digits=3)


## Panel B
# Zoom in on initial phase where inversion relative fitness is greater than 1
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(subDat2)), ylim = c(0,max(subDat2$I.t)*1.15), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        axis(1, las=1)
        axis(4, las=1, cex=1.3)
        # Inversion Frequency
        lines(I.t ~ seq_along(I.t), type='l', lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=subDat2)
        lines(I.gensol ~ seq_along(I.gensol), type='l', lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=subDat2)
        # Benchmark qstar
        abline(h=qstr, lty=2, lwd=1.5, col='tomato')
        # Inversion relative fitness
        par(new=TRUE)
        plot(invRelFit ~ seq_along(invRelFit), type='l', ylim=c((1-(max(invRelFit) - 1)), (1 + subDat$U[1]*subDat$x[1]*1.25)), axes=FALSE, lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=subDat2, ylab="", xlab='')
        lines(r.0.rel.w.I.t ~ seq_along(r.0.rel.w.I.t), type='l', lwd=2, lty=3, col=transparentColor('tomato', opacity=1))
        # Benchmark for relative fitness of 1
        abline(h=1, lty=2, lwd=1.5, col=COLS[1])
        # axes
        axis(2, las=1, cex=1.3)
        # Annotate 1 + Ux approximation
        points((1 + subDat$U[1]*subDat$x[1]) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        # Plot labels etc.
        if(i == 10){        
                    proportionalLabel(1.3,  0.5,   expression(paste("Inversion Frequency (", italic(Y[I]), ")")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=270)
                }
        if(i == 1){
                    proportionalLabel(-0.3,  0.5,  expression(paste("Inversion Relative Fitness")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)
                }
        proportionalLabel(0.15,  0.95,  expression(paste(italic(tilde(s)[WF])/italic(tilde(s)[Det.])," = ")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(0.35, 0.95, substitute(xx, list(xx=relGeomS)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.7, 0.2,  expression(paste(italic(tilde(s)[WF])," = ")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(0.85, 0.2, substitute(xx, list(xx=tildeSWF)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.7, 0.1,  expression(paste(italic(tilde(s)[Det.])," = ")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=0)
        proportionalLabel(0.85, 0.1, substitute(xx, list(xx=tildeSDet)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        usr  <-  par('usr')
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(paste("Inv. Rel. Fit.")),
                            expression(paste("Determ. Inv. Rel. Fit.")),
                            expression(italic(I[t]))),
               lty     =  c(1,3,1),
               lwd     =  2,
               col     =  c(transparentColor(COLS[2], opacity=1),
                            transparentColor('tomato', opacity=1),
                            transparentColor("#252525", opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


## Panel C
# Deleterious mutation frequency dynamics
    # make plot
## Panel C
# Deleterious mutation frequency dynamics
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(subDat2)), ylim = c(min(subDat2$qt.W/(1 - subDat2$I.t)), max(subDat2$qt.W/(1 - subDat2$I.t))), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        axis(1, las=1)
        axis(2, las=1, cex=1.3)
        # Frequency of del. mutations on Xf
        lines(qt.W/(1 - I.t) ~ seq_along(qt.W), type='l', lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=subDat2)
        # Benchmark initial equilibrium frequency
        abline(h=(subDat2$qt.W[1]/(1 - subDat2$I.t[1])), lty=2, lwd=1.5, col=COLS[1])
        # Frequency of del. mutations on YI
        par(new=TRUE)
        plot(qt.I.W/I.t ~ seq_along(qt.I.W), type='l', ylim = c(0, max(subDat2$qt.I.W/subDat2$I.t)), lwd=2, lty=1, col=transparentColor('red', opacity=1), data=subDat2, axes=FALSE, ylab="", xlab='')
        # axes
        axis(4, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        if(i == 10){        
        proportionalLabel(1.3,  0.5,   expression(paste(italic(q[I]))), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=270)
                }
        if(i == 1){
            proportionalLabel(-0.3,  0.5, expression(paste(italic(q[N]^W))), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)
                }    
        proportionalLabel( 0.5,  -0.15,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        usr  <-  par('usr')
        legend(
               x       =  usr[2]*0.25,
               y       =  usr[4],
               legend  =  c(
                            expression(paste(italic(q[I]))),
                            expression(paste(italic(q[N]^W)))),
               lty     =  1,
               lwd     =  2,
               col     =  c(transparentColor("red", opacity=1),
                            transparentColor("#252525", opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

    }

}



##################################


deterministicMutRateIllustration_SLR_Nei  <-  function(wHap=FALSE) {

    # set plot layout
    layout.mat <- matrix(c(1:9), nrow=3, ncol=3, byrow=FALSE)
    layout     <- layout(layout.mat,respect=TRUE)

    # Make data for plotting 
    s        <-  0.01
    h        <-  0.25
    nTot     <-  10^4 

    ###########
    # U = 0.02
    U        <-  0.02
    u        <-  U/nTot
    x        <-  0.2

    # make data for plotting
    r.0.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 0, x = x, h = h, generations = 2500)
    r.1.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 1, x = x, h = h, generations = 2500)
    r.0.rel.w.YI.t  <-  r.0.Dat$rel.w.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.1.rel.w.YI.t  <-  r.1.Dat$rel.w.YI.t[round(r.1.Dat$YI.t,digits=4) > 0]

    r.0.Hap.Dat  <-  makeDeterministicSimHaploidData(U = U, s = s, r = 0, x = x, h = h, generations = 2500)
    r.1.Hap.Dat  <-  makeDeterministicSimHaploidData(U = U, s = s, r = 1, x = x, h = h, generations = 2500)

   # Colors
    COLS  <- colorRampPalette(wes_palette("Zissou1"))(4)

    # Axes limits
    relFitMin  <-  0.997
    relFitMax  <-  1 + U*x*1.25


## Panel A
# Effect of initial del. mut. load on inversion relative fitness 
    # make plot
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(relFitMin, relFitMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Inversion relative frequencies
        lines(r.0.rel.w.YI.t ~ seq_along(r.0.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1))
        lines(r.1.rel.w.YI.t ~ seq_along(r.1.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1))
        # Inversion relative fitness -- Haploid model
        if(wHap) {
            lines(rel.w.YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.0.Hap.Dat)
            lines(rel.w.YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.1.Hap.Dat)
        }
        # Benchmark inversion relative fitness of 1
        abline(h=1, lwd=1, lty=2, col=1)
        # Points depicting (arbitrary) extinction frequency of 0.00001
        points(r.1.rel.w.YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        # Annotate 1 + Ux approximation
        points((1 + U*x) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        proportionalArrows(0.125, 0.85, 0.06, 0.85, cex=1, adj=c(0.5, 0.5), xpd=NA, log=FALSE, length=0.05)
        proportionalLabel(0.25, 0.85, expression(1+italic(Ux)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, expression(paste(italic(Ux/(hs)==1.6))), cex=2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,   expression(paste("Inversion Relative Fitness")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(d)==0),
                            expression(italic(d)==1)),
               lty     =  1,
               lwd     =  2,
               col     =  c(transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
        legend(
               x       =  usr[2],
               y       =  usr[4]*0.998,
               legend  =  c(expression("Haploid Model")),
               lty     =  3,
               lwd     =  2,
               col     =  1,
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


qstr       <-  qstar(U=U, q.0=r.0.Dat$YI.t[1], x=x, s=s, h=h)
t          <-  seq_along(r.0.rel.w.YI.t)
YI.gensol  <-  YI.genSol(U=U, YI.init=r.0.Dat$YI.t[1], x=x, s=s, h=h, t=t)

## Panel D
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,(qstr*1.25)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Inversion Frequencies
        lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        # General solution for r = 0
        lines(YI.gensol ~ t, lwd=2, lty=3, col="tomato")
        if(wHap) {
            # Inversion frequencies -- Haploid model
            lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.0.Hap.Dat)
            lines(YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.1.Hap.Dat)
        }
        # Points
        points(YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        # benchmark qstar
        abline(h=qstr, lty=2)
        proportionalLabel(0.9, 0.85, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,   expression(paste("Inversion Frequency (", italic(Y[I]), ")")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        # Legend
        legend(
               x       =  usr[2]*0.95,
               y       =  usr[4],
               legend  =  c(
                            expression(paste("General solution for ", italic(Y[I]), " when ", italic(d)==0))),
               lty     =  3,
               lwd     =  2,
               col     =  "tomato",
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
        if(wHap) {
        legend(
               x       =  usr[2]*0.6,
               y       =  usr[4]*0.9,
               legend  =  c(
                            expression(paste("Haploid Model for ", italic(Y[I])))),
               lty     =  3,
               lwd     =  2,
               col     =  1,
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
        }

    # reset Colors
     COLS  <-  c('red', '#252525')

## Panel G
# deleterious allele frequencies for r = 1
     # Make plot
    delMax  <-  r.1.Dat$qt.Xf.W[1]*1.05
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(qt.Xf.W[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.Xf.W[1:length(r.1.rel.w.YI.t)]), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qt.Xf.D[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.Xf.D[1:length(r.1.rel.w.YI.t)]), lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qt.YI.W[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.YI.W[1:length(r.1.rel.w.YI.t)]), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
        # Points indicating inversion extinctions
        points(qt.Xf.W[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(qt.Xf.D[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(qt.YI.W[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.1.Dat)
        # axes
        axis(1, las=1, cex=1.3)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'G', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(d)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,   expression(paste("Deleterious Allele Frequency")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend = c(expression(italic(q)[italic(X[f])]^italic(W)), 
                          expression(italic(q)[italic(Y[I])]^italic(W)), 
                          expression(italic(q)[italic(X[f])]^italic(D))),
               lty     =  c(1,1,3),
               lwd     =  2,
               col     =  c(transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )



    ###########
    # U = 0.05
    U  <-  0.05
    r.0.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 0, x = x, h = h, generations = 5000)
    r.1.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 1, x = x, h = h, generations = 5000)
    r.2.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 2, x = x, h = h, generations = 5000)
    r.3.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 3, x = x, h = h, generations = 5000)

    r.0.rel.w.YI.t  <-  r.0.Dat$rel.w.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.1.rel.w.YI.t  <-  r.1.Dat$rel.w.YI.t[round(r.1.Dat$YI.t,digits=4) > 0]
    r.2.rel.w.YI.t  <-  r.2.Dat$rel.w.YI.t[round(r.2.Dat$YI.t,digits=4) > 0]
    r.3.rel.w.YI.t  <-  r.3.Dat$rel.w.YI.t[round(r.3.Dat$YI.t,digits=4) > 0]

    r.0.Hap.Dat  <-  makeDeterministicSimHaploidData(U = U, s = s, r = 0, x = x, h = h, generations = 5000)
    r.1.Hap.Dat  <-  makeDeterministicSimHaploidData(U = U, s = s, r = 1, x = x, h = h, generations = 5000)
    r.2.Hap.Dat  <-  makeDeterministicSimHaploidData(U = U, s = s, r = 2, x = x, h = h, generations = 5000)
    r.3.Hap.Dat  <-  makeDeterministicSimHaploidData(U = U, s = s, r = 3, x = x, h = h, generations = 5000)

   # Colors
    COLS     <-  colorRampPalette(wes_palette("Zissou1"))(4)

    relFitMin  <-  min(r.3.Dat$rel.w.YI.t)
    relFitMax  <-  1 + U*x*1.25

## Panel B
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(relFitMin,relFitMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(r.3.rel.w.YI.t ~ seq_along(r.3.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1))
        lines(r.2.rel.w.YI.t ~ seq_along(r.2.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.1.rel.w.YI.t ~ seq_along(r.1.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r.0.rel.w.YI.t ~ seq_along(r.0.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1))
        # Inversion relative fitness -- Haploid model
        if(wHap) {
            lines(rel.w.YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.0.Hap.Dat)
            lines(rel.w.YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.1.Hap.Dat)
            lines(rel.w.YI.t[1:length(r.2.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.2.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.2.Hap.Dat)
            lines(rel.w.YI.t[1:length(r.3.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.3.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.3.Hap.Dat)
        }
        # Benchmark inversion relative fitness of 1
        abline(h=1, lwd=1, lty=2, col=1)
        # Points depicting (arbitrary) extinction frequency of 0.00001
        points(r.1.rel.w.YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r.2.rel.w.YI.t[length(r.2.rel.w.YI.t)] ~ length(r.2.rel.w.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
        points(r.3.rel.w.YI.t[length(r.3.rel.w.YI.t)] ~ length(r.3.rel.w.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        # Annotate 1 + Ux approximation
        points((1 + U*x) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, expression(paste(italic(Ux/(hs)==4))), cex=2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(d)==0),
                            expression(italic(d)==1),
                            expression(italic(d)==2),
                            expression(italic(d)==3)),
               lty     =  1,
               lwd     =  2,
               col     =  c(transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[4], opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )



qstr       <-  qstar(U=U, q.0=r.0.Dat$YI.t[1], x=x, s=s, h=h)
t          <-  seq_along(r.0.rel.w.YI.t)
YI.gensol  <-  YI.genSol(U=U, YI.init=r.0.Dat$YI.t[1], x=x, s=s, h=h, t=t)

## Panel E
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
#    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,(qstr*1.25)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(YI.t[1:length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.2.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.2.Dat)
        lines(YI.t[1:length(r.3.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.3.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.3.Dat)
        # General solution for r = 0
        lines(YI.gensol ~ t, lwd=2, lty=3, col="tomato")
        if(wHap) {
            # Inversion frequencies -- Haploid model
            lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.0.Hap.Dat)
            lines(YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.1.Hap.Dat)
            lines(YI.t[1:length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.2.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.2.Hap.Dat)
            lines(YI.t[1:length(r.3.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.3.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.3.Hap.Dat)
        }
        # Points
        points(YI.t[length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.1.rel.w.YI.t)], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(YI.t[length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.2.rel.w.YI.t)], pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.2.Dat)
        points(YI.t[length(r.3.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.3.rel.w.YI.t)], pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.3.Dat)
        # benchmark qstar
        abline(h=qstr, lty=2)
        proportionalLabel(0.9, 0.85, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'E', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # reset Colors
     COLS  <-  c('red', '#252525')

## Panel H
# deleterious allele frequencies for r = 5
     # Make plot
     delMax  <-  r.1.Dat$qt.Xf.W[1]*1.05
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(qt.Xf.W[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.Xf.W[1:length(r.1.rel.w.YI.t)]), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qt.Xf.D[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.Xf.D[1:length(r.1.rel.w.YI.t)]), lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qt.YI.W[1:length(r.1.rel.w.YI.t)] ~ seq_along(qt.YI.W[1:length(r.1.rel.w.YI.t)]), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
        # Points indicating inversion extinctions
        points(qt.Xf.W[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(qt.Xf.D[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(qt.YI.W[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.1.Dat)
        # axes
        axis(1, las=1, cex=1.3)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'H', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(d)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        




    ###########
    # U = 0.1
    U  <-  0.1
    r.0.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 0,  x = x, h = h, generations = 5000)
    r.2.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 2,  x = x, h = h, generations = 5000)
    r.4.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 4, x = x, h = h, generations = 5000)
    r.6.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 6, x = x, h = h, generations = 5000)

    r.0.rel.w.YI.t  <-  r.0.Dat$rel.w.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.2.rel.w.YI.t  <-  r.2.Dat$rel.w.YI.t[round(r.2.Dat$YI.t,digits=4) > 0]
    r.4.rel.w.YI.t  <-  r.4.Dat$rel.w.YI.t[round(r.4.Dat$YI.t,digits=4) > 0]
    r.6.rel.w.YI.t  <-  r.6.Dat$rel.w.YI.t[round(r.6.Dat$YI.t,digits=4) > 0]

    r.0.Hap.Dat  <-  makeDeterministicSimHaploidData(U = U, s = s, r = 0,  x = x, h = h, generations = 5000)
    r.2.Hap.Dat  <-  makeDeterministicSimHaploidData(U = U, s = s, r = 2,  x = x, h = h, generations = 5000)
    r.4.Hap.Dat  <-  makeDeterministicSimHaploidData(U = U, s = s, r = 4, x = x, h = h, generations = 5000)
    r.6.Hap.Dat  <-  makeDeterministicSimHaploidData(U = U, s = s, r = 6, x = x, h = h, generations = 5000)

   # Colors
    COLS     <-  colorRampPalette(wes_palette("Zissou1"))(4)

    relFitMin  <-  0.985
    relFitMax  <-  1 + U*x*1.25

## Panel C
# Effect of initial del. mut. load on inversion relative fitness 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(relFitMin,relFitMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Inversion relative fitnesses
        lines(r.0.rel.w.YI.t ~ seq_along(r.0.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1))
        lines(r.2.rel.w.YI.t ~ seq_along(r.2.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r.4.rel.w.YI.t ~ seq_along(r.4.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.6.rel.w.YI.t ~ seq_along(r.6.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1))
        # Inversion relative fitness -- Haploid model
        if(wHap) {
            lines(rel.w.YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.0.Hap.Dat)
            lines(rel.w.YI.t[1:length(r.2.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.2.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.2.Hap.Dat)
            lines(rel.w.YI.t[1:length(r.4.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.4.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.4.Hap.Dat)
            lines(rel.w.YI.t[1:length(r.6.rel.w.YI.t)] ~ seq_along(rel.w.YI.t)[1:length(r.6.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.6.Hap.Dat)
        }
        # Benchmakr relative fitness of 1
        abline(h=1, lwd=1, lty=2, col=1)
        # Points depicting (arbitrary) extinction frequency of 0.00001
        points(r.0.rel.w.YI.t[length(r.0.rel.w.YI.t)] ~ length(r.0.rel.w.YI.t), pch=8, col='#252525')
        points(r.2.rel.w.YI.t[length(r.2.rel.w.YI.t)] ~ length(r.2.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r.4.rel.w.YI.t[length(r.4.rel.w.YI.t)] ~ length(r.4.rel.w.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
        points(r.6.rel.w.YI.t[length(r.6.rel.w.YI.t)] ~ length(r.6.rel.w.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        # Annotate 1 + Ux approximation
        points((1+U*x) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, expression(paste(italic(Ux/(hs)==8))), cex=2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(d)==0),
                            expression(italic(d)==2),
                            expression(italic(d)==4),
                            expression(italic(d)==6)),
               lty     =  1,
               lwd     =  2,
               col     =  c(transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[4], opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


qstr  <-  qstar(U=U, q.0=r.0.Dat$YI.t[1], x=x, s=s, h=h)
t          <-  seq_along(r.0.rel.w.YI.t)
YI.gensol  <-  YI.genSol(U=U, YI.init=r.0.Dat$YI.t[1], x=x, s=s, h=h, t=t)

## Panel F
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,(qstr*1.25)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.2.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.2.Dat)
        lines(YI.t[1:length(r.4.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.4.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.4.Dat)
        lines(YI.t[1:length(r.6.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.6.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.6.Dat)
            # General solution for r = 0
            lines(YI.gensol ~ t, lwd=2, lty=3, col="tomato")
        if(wHap) {
            # Inversion frequencies -- Haploid model
            lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.0.Hap.Dat)
            lines(YI.t[1:length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.2.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.2.Hap.Dat)
            lines(YI.t[1:length(r.4.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.4.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.4.Hap.Dat)
            lines(YI.t[1:length(r.6.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.6.rel.w.YI.t)], lwd=2, lty=3, col=transparentColor("#252525", opacity=1), data=r.6.Hap.Dat)
        }
        # Points
        points(YI.t[length(r.0.rel.w.YI.t)] ~ length(r.0.rel.w.YI.t), pch=8, col='#252525', data=r.0.Dat)
        points(YI.t[length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.2.rel.w.YI.t)], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.2.Dat)
        points(YI.t[length(r.4.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.4.rel.w.YI.t)], pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.4.Dat)
        points(YI.t[length(r.6.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.6.rel.w.YI.t)], pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.6.Dat)
        # benchmark qstar
        abline(h=qstr, lty=2)
        proportionalLabel(0.9, 0.85, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'F', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # reset Colors
     COLS  <-  c('red', '#252525')

## Panel I
# deleterious allele frequencies for r = 5
     # Make plot
    delMax  <-  r.2.Dat$qt.Xf.W[1]*1.05
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(qt.Xf.W[1:length(r.2.rel.w.YI.t)] ~ seq_along(qt.Xf.W[1:length(r.2.rel.w.YI.t)]), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.2.Dat)
        lines(qt.Xf.D[1:length(r.2.rel.w.YI.t)] ~ seq_along(qt.Xf.D[1:length(r.2.rel.w.YI.t)]), lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.2.Dat)
        lines(qt.YI.W[1:length(r.2.rel.w.YI.t)] ~ seq_along(qt.YI.W[1:length(r.2.rel.w.YI.t)]), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.2.Dat)
        # Points indicating inversion extinctions
        points(qt.Xf.W[length(r.2.rel.w.YI.t)] ~ length(r.2.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.2.Dat)
        points(qt.Xf.D[length(r.2.rel.w.YI.t)] ~ length(r.2.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.2.Dat)
        points(qt.YI.W[length(r.2.rel.w.YI.t)] ~ length(r.2.rel.w.YI.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.2.Dat)
        # axes
        axis(1, las=1, cex=1.2)
        axis(2, las=1, cex=1.2)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'I', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(d)==2), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        

}







############################
##  Perspective Piece Figs

perspectiveInvRelFit  <-  function() {

    # set plot layout
    layout.mat <- matrix(c(1,1,1,1,
                           1,1,1,1,
                           2,2,2,2,
                           2,2,2,2), nrow=4, ncol=4, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

    # Make data for plotting 
    s        <-  0.01
    h        <-  0.25
    nTot     <-  10^4 

    ###########
    # U = 0.02
    U        <-  0.05
    u        <-  U/nTot
    x        <-  0.2

    relFitMin  <-  0.997
    relFitMax  <-  1.005


    ###########
    # h = 0.25
    # Nei-style SLR model
    r.0.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 0, x = x, h = h, generations = 2500)
    r.1.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 1, x = x, h = h, generations = 2500)
    r.2.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 2, x = x, h = h, generations = 2500)
#    r.3.Dat  <-  makeDeterministicFigSimData_SLR_Nei(U = U, s = s, r = 3, x = x, h = h, generations = 2500)

    r.0.rel.w.YI.t  <-  r.0.Dat$rel.w.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.1.rel.w.YI.t  <-  r.1.Dat$rel.w.YI.t[round(r.1.Dat$YI.t,digits=4) > 0]
    r.2.rel.w.YI.t  <-  r.2.Dat$rel.w.YI.t[round(r.2.Dat$YI.t,digits=4) > 0]
#    r.3.rel.w.YI.t  <-  r.3.Dat$rel.w.YI.t[round(r.3.Dat$YI.t,digits=4) > 0]

   # Colors
    COLS     <-  colorRampPalette(wes_palette("Zissou1"))(3)

    relFitMin  <-  0.9925
    relFitMax  <-  1 + U*x*1.25

    # Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
    par(omi=c(0.5, 0.5, 0.5, 0.0), mar = c(3,6,1,0), bty='o', xaxt='s', yaxt='s')
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(relFitMin,relFitMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
#        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
#        plotGrid(lineCol='grey80')
#        box()
        # lines
#        lines(r.3.rel.w.YI.t ~ seq_along(r.3.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1))
        lines(r.2.rel.w.YI.t ~ seq_along(r.2.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.1.rel.w.YI.t ~ seq_along(r.1.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r.0.rel.w.YI.t ~ seq_along(r.0.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1))
        # Benchmark inversion relative fitness of 1
        abline(h=1, lwd=1, lty=2, col=1)
        # Points depicting (arbitrary) extinction frequency of 0.00001
        points(r.1.rel.w.YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r.2.rel.w.YI.t[length(r.2.rel.w.YI.t)] ~ length(r.2.rel.w.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
#        points(r.3.rel.w.YI.t[length(r.3.rel.w.YI.t)] ~ length(r.3.rel.w.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        # Annotate 1 + Ux approximation
#        points((1 + U*x) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, cex=1.2)
        # Plot labels etc.
#        proportionalLabel(0.04, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, expression(paste("Mean # mutations"), cex=2, adj=c(0.5, 0.5), xpd=NA))
        proportionalLabel(-0.15,  0.5,   expression(paste("Inversion relative fitness")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  1.012,
               legend  =  c(
                            expression(italic(d)==0),
                            expression(italic(d)==1),
                            expression(italic(d)==2)),
#                            expression(italic(d)==3)),
               lty     =  1,
               lwd     =  2,
               col     =  c(transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[3], opacity=1)),
#                            transparentColor(COLS[4], opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


    # Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,(0.025)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
#        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
#        plotGrid(lineCol='grey80')
#        box()
        # lines
        lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(YI.t[1:length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.2.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.2.Dat)
#        lines(YI.t[1:length(r.3.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.3.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.3.Dat)
        # Points
        points(YI.t[length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.1.rel.w.YI.t)], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(YI.t[length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.2.rel.w.YI.t)], pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.2.Dat)
#        points(YI.t[length(r.3.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.3.rel.w.YI.t)], pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.3.Dat)
        # benchmark qstar
#        abline(h=qstr, lty=2)
#        proportionalLabel(0.9, 0.85, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, cex=1.2)
        axis(2, las=1, cex=1.2)
        # Plot labels etc.
        proportionalLabel(-0.15,  0.5,   expression(paste("Inversion frequency")), cex=1.5, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)        

}