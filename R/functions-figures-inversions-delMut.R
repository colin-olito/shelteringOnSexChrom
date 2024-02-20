#  Functions to create figures for inversions expanding SLR 
#  in the presence of deleterious mutational variation
source('./R/functions-figures.R')
source('./R/simulations-inversions-delMut.R')



 qstar  <-  function(YI.init, U, x, s, h) {
     (YI.init*exp((U*x)/(1 - exp(-s*h)))) / (1 - YI.init*(1 - exp((U*x)/(1 - exp(-h*s)))))
 }

# qstar  <-  function(YI.init, U, x, s, h) {
#     YI.init*exp((U*x)/(1 - exp(-s*h)))
# }

# FIGURE 1. Illustration of deterministic dynamics
#           for different dominance values
deterministicDominanceIllustration  <-  function() {

    # set plot layout
    layout.mat <- matrix(c(1:9), nrow=3, ncol=3, byrow=FALSE)
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
    r.0.Dat  <-  makeDeterministicFigSimData(U = U, s = s, r = 0, x = x, h = h, generations = 2500)
    r.1.Dat  <-  makeDeterministicFigSimData(U = U, s = s, r = 1, x = x, h = h, generations = 2500)
#    r.2.Dat  <-  makeDeterministicFigSimData(s = s, r = 2, x = x, h = h, generations = 2500)
#    r.3.Dat  <-  makeDeterministicFigSimData(s = s, r = 3, x = x, h = h, generations = 2500)

    r.0.rel.w.YI.t  <-  r.0.Dat$rel.w.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.1.rel.w.YI.t  <-  r.1.Dat$rel.w.YI.t[round(r.1.Dat$YI.t,digits=4) > 0]
#    r.2.rel.w.YI.t  <-  r.2.Dat$rel.w.YI.t[round(r.2.Dat$YI.t,digits=4) > 0]
#    r.3.rel.w.YI.t  <-  r.3.Dat$rel.w.YI.t[round(r.3.Dat$YI.t,digits=4) > 0]

    autotime  <-  c(1:2500)
    q.wt.Auto  <-  p.t(u=10^-6, sdHom=0.01, h=h, t=autotime) 

   # Colors
#    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
#    COLS     <-  rev(colfunc(5))
#COLS  <-  colorRampPalette(wes_palette("Zissou1"), interpolate = c("linear"))(4)
COLS  <- colorRampPalette(wes_palette("Zissou1"))(4)
    # Axes limits
    relFitMin  <-  0.996
    relFitMax  <-  1.005


## Panel A
# Effect of initial del. mut. load on inversion relative fitness 
    # make plot
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(relFitMin, relFitMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        abline(h=1, lwd=1, lty=2, col=1)
#        abline(h=(1+U*x), lwd=1, lty=3, col=2)
#        proportionalLabel(0.15, 0.89, expression(italic(1+italic(Ux))), cex=1, adj=c(0.5, 0.5), xpd=NA)
        lines(r.0.rel.w.YI.t ~ seq_along(r.0.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1))
        lines(r.1.rel.w.YI.t ~ seq_along(r.1.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1))
#        lines(r.2.rel.w.YI.t ~ seq_along(r.2.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1))
#        lines(r.3.rel.w.YI.t ~ seq_along(r.3.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1))
        # Points depicting (arbitrary) extinction frequency of 0.00001
#        points(r.0.rel.w.YI.t[length(r.0.rel.w.YI.t)] ~ length(r.0.rel.w.YI.t), pch=8, col='#252525')
        points(r.1.rel.w.YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
#        points(r.2.rel.w.YI.t[length(r.2.rel.w.YI.t)] ~ length(r.2.rel.w.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
#        points(r.3.rel.w.YI.t[length(r.3.rel.w.YI.t)] ~ length(r.3.rel.w.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        # Annotate 1 + Ux approximation
        points((1 + U*x) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
#        abline(h=(1+U*x), lwd=1, lty=3, col=2)
        proportionalArrows(0.125, 0.86, 0.06, 0.86, cex=1, adj=c(0.5, 0.5), xpd=NA, log=FALSE, length=0.05)
        proportionalLabel(0.25, 0.86, expression(1+italic(Ux)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.5, 1.15, substitute(italic(h)==xx, list(xx=h)), cex=2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, expression(paste(italic(Ux/hs==1.6))), cex=2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,   expression(paste("Inversion Relative Fitness")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(r)==0),
                            expression(italic(r)==1)),
#                            expression(italic(r)==2)),
#                            expression(italic(r)==3)),
               lty     =  1,
               lwd     =  2,
               col     =  c(transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1)),
#                            transparentColor(COLS[3], opacity=1)),
#                            transparentColor(COLS[2], opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

qstr  <-  qstar(U=U, YI.init=r.0.Dat$YI.t[1], x=x, s=s, h=h)

## Panel D
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,0.1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
#        lines(YI.t[1:length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.2.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.2.Dat)
#        lines(YI.t[1:length(r.3.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.3.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.3.Dat)
        # Points
#        points(YI.t[length(r.0.rel.w.YI.t)] ~ length(r.0.rel.w.YI.t), pch=8, col='#252525', data=r.0.Dat)
        points(YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
#        points(YI.t[length(r.2.rel.w.YI.t)] ~ length(r.2.rel.w.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.2.Dat)
#        points(YI.t[length(r.3.rel.w.YI.t)] ~ length(r.3.rel.w.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.3.Dat)
        # benchmark qstar
        abline(h=qstr, lty=2)
        proportionalLabel(0.9, 0.1, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,   expression(paste("Inversion Frequency (", italic(Y[I]), ")")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        

    # reset Colors
     COLS  <-  c('red',
                 '#252525')
## Panel G
# deleterious allele frequencies for r = 1
     # Make plot
    delMax  <-  r.1.Dat$XaOv.wt.t[1]*1.05
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(XaOv.wt.t ~ seq_along(XaOv.wt.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(XaOv.del.t ~ seq_along(XaOv.del.t), lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qI.wt.t ~ seq_along(qI.wt.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
#        lines(q.wt.Auto ~ autotime, lwd=2, lty=1, col="tomato")
        # axes
        axis(1, las=1, cex=1.3)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'G', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(r)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,   expression(paste("Deleterious Allele Frequency")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4]*0.9,
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
    # h = 0.1
#   s        <-  0.025
#    h        <-  0.1
#    Ufactor  <-  2
#    nTot     <-  10^4 
#    U        <-  Ufactor*s
#    u        <-  U/nTot
#    x        <-  0.2


    # make data for plotting
#U  <-  0.05
h  <-  0.1
    r.0.Dat  <-  makeDeterministicFigSimData(U = U, s = s, r = 0, x = x, h = h, generations = 8000)
    r.1.Dat  <-  makeDeterministicFigSimData(U = U, s = s, r = 1, x = x, h = h, generations = 8000)
    r.2.Dat  <-  makeDeterministicFigSimData(U = U, s = s, r = 2, x = x, h = h, generations = 8000)
    r.3.Dat  <-  makeDeterministicFigSimData(U = U, s = s, r = 3, x = x, h = h, generations = 8000)

    r.0.rel.w.YI.t  <-  r.0.Dat$rel.w.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.1.rel.w.YI.t  <-  r.1.Dat$rel.w.YI.t[round(r.1.Dat$YI.t,digits=4) > 0]
    r.2.rel.w.YI.t  <-  r.2.Dat$rel.w.YI.t[round(r.2.Dat$YI.t,digits=4) > 0]
    r.3.rel.w.YI.t  <-  r.3.Dat$rel.w.YI.t[round(r.3.Dat$YI.t,digits=4) > 0]

    autotime  <-  c(1:8000)
     q.wt.Auto  <-  p.t(u=10^-6, sdHom=0.01, h=h, t=autotime) 

   # Colors
#    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
#    COLS     <-  rev(colfunc(5))
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
        abline(h=1, lwd=1, lty=2, col=1)
        lines(r.3.rel.w.YI.t ~ seq_along(r.3.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1))
        lines(r.2.rel.w.YI.t ~ seq_along(r.2.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.1.rel.w.YI.t ~ seq_along(r.1.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r.0.rel.w.YI.t ~ seq_along(r.0.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1))
        # Points depicting (arbitrary) extinction frequency of 0.00001
#        points(r.0.rel.w.YI.t[length(r.0.rel.w.YI.t)] ~ length(r.0.rel.w.YI.t), pch=8, col='#252525')
        points(r.1.rel.w.YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r.2.rel.w.YI.t[length(r.2.rel.w.YI.t)] ~ length(r.2.rel.w.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
        points(r.3.rel.w.YI.t[length(r.3.rel.w.YI.t)] ~ length(r.3.rel.w.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        # Annotate 1 + Ux approximation
        points((1 + U*x) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
#        abline(h=(1+U*x), lwd=1, lty=3, col=2)
#        proportionalArrows(0.125, 0.86, 0.06, 0.86, cex=1, adj=c(0.5, 0.5), xpd=NA, log=FALSE, length=0.05)
#        proportionalLabel(0.2, 0.86, expression(italic(1+italic(Ux))), cex=1, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#       proportionalLabel(0.5, 1.15, substitute(italic(h)==xx, list(xx=h)), cex=2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, expression(paste(italic(Ux/hs==4))), cex=2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(r)==0),
                            expression(italic(r)==1),
                            expression(italic(r)==2),
                            expression(italic(r)==3)),
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



qstr  <-  qstar(U=U, YI.init=r.0.Dat$YI.t[1], x=x, s=s, h=h)

## Panel E
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
#    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,0.1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(YI.t[1:length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.2.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.2.Dat)
        lines(YI.t[1:length(r.3.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.3.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.3.Dat)
        # Points
#        points(YI.t[length(r.0.rel.w.YI.t)] ~ length(r.0.rel.w.YI.t), pch=8, col='#252525', data=r.0.Dat)
        points(YI.t[length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.1.rel.w.YI.t)], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(YI.t[length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.2.rel.w.YI.t)], pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.2.Dat)
        points(YI.t[length(r.3.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.3.rel.w.YI.t)], pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.3.Dat)
        # benchmark qstar
        abline(h=qstr, lty=2)
        proportionalLabel(0.9, 0.28, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'E', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # reset Colors
     COLS  <-  c('red',
                 '#252525')
## Panel H
# deleterious allele frequencies for r = 5
     # Make plot
     delMax  <-  r.1.Dat$XaOv.wt.t[1]*1.05
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(XaOv.wt.t ~ seq_along(XaOv.wt.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(XaOv.del.t ~ seq_along(XaOv.del.t), lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qI.wt.t ~ seq_along(qI.wt.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
#        lines(q.wt.Auto ~ autotime, lwd=2, lty=1, col="tomato")
        # axes
        axis(1, las=1, cex=1.3)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'H', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(r)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        




    ###########
    # h = 0.01
#   s        <-  0.25
#    h        <-  0.01
#    Ufactor  <-  2
#    nTot     <-  10^4 
#    U        <-  Ufactor*s
#    u        <-  U/nTot
#    x        <-  0.2

    
    # make data for plotting
#U  <-  0.1
 h  <-  0.01
    r.0.Dat   <-  makeDeterministicFigSimData(U = U, s = s, r = 0,  x = x, h = h, generations = 40000)
    r.5.Dat   <-  makeDeterministicFigSimData(U = U, s = s, r = 5,  x = x, h = h, generations = 40000)
    r.10.Dat  <-  makeDeterministicFigSimData(U = U, s = s, r = 10, x = x, h = h, generations = 40000)
    r.15.Dat  <-  makeDeterministicFigSimData(U = U, s = s, r = 15, x = x, h = h, generations = 40000)
    r.20.Dat  <-  makeDeterministicFigSimData(U = U, s = s, r = 20, x = x, h = h, generations = 40000)

    r.0.rel.w.YI.t   <-  r.0.Dat$rel.w.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.5.rel.w.YI.t   <-  r.5.Dat$rel.w.YI.t[round(r.5.Dat$YI.t,digits=4) > 0]
    r.10.rel.w.YI.t  <-  r.10.Dat$rel.w.YI.t[round(r.10.Dat$YI.t,digits=4) > 0]
    r.15.rel.w.YI.t  <-  r.15.Dat$rel.w.YI.t[round(r.15.Dat$YI.t,digits=4) > 0]
    r.20.rel.w.YI.t  <-  r.20.Dat$rel.w.YI.t[round(r.20.Dat$YI.t,digits=4) > 0]

    autotime  <-  c(1:40000)
    q.wt.Auto  <-  p.t(u=10^-6, sdHom=0.01, h=h, t=autotime) 

   # Colors
#    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
#    COLS     <-  rev(colfunc(5))
    COLS     <-  colorRampPalette(wes_palette("Zissou1"))(4)

    qHat  <-  r.0.Dat$qY.wt.t[1]
    rBar  <-  round(nTot*x*qHat)


## Panel C
# Effect of initial del. mut. load on inversion relative fitness 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(relFitMin,relFitMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        abline(h=1, lwd=1, lty=2, col=1)
#        abline(h=(1+U*x), lwd=1, lty=3, col=2)
#        proportionalLabel(0.15, 0.89, expression(italic(1+italic(Ux))), cex=1, adj=c(0.5, 0.5), xpd=NA)
        lines(r.0.rel.w.YI.t ~ seq_along(r.0.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1))
        lines(r.5.rel.w.YI.t ~ seq_along(r.5.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r.10.rel.w.YI.t ~ seq_along(r.10.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.15.rel.w.YI.t ~ seq_along(r.15.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1))
#        lines(r.20.rel.w.YI.t ~ seq_along(r.20.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[5], opacity=1))
        # Points depicting (arbitrary) extinction frequency of 0.00001
        points(r.0.rel.w.YI.t[length(r.0.rel.w.YI.t)] ~ length(r.0.rel.w.YI.t), pch=8, col='#252525')
        points(r.5.rel.w.YI.t[length(r.5.rel.w.YI.t)] ~ length(r.5.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r.10.rel.w.YI.t[length(r.10.rel.w.YI.t)] ~ length(r.10.rel.w.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
        points(r.15.rel.w.YI.t[length(r.15.rel.w.YI.t)] ~ length(r.15.rel.w.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
#        points(r.20.rel.w.YI.t[length(r.20.rel.w.YI.t)] ~ length(r.20.rel.w.YI.t), pch=21, col=transparentColor(COLS[5], opacity=1), bg=transparentColor(COLS[5], opacity=0.7))
        # Annotate 1 + Ux approximation
        points((1 + U*x) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
#        abline(h=(1+U*x), lwd=1, lty=3, col=2)
#        proportionalArrows(0.125, 0.86, 0.06, 0.86, cex=1, adj=c(0.5, 0.5), xpd=NA, log=FALSE, length=0.05)
#        proportionalLabel(0.2, 0.86, expression(italic(1+italic(Ux))), cex=1, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.5, 1.15, substitute(italic(h)==xx, list(xx=h)), cex=2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, expression(paste(italic(Ux/hs==40))), cex=2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(r)==0),
                            expression(italic(r)==5),
                            expression(italic(r)==10),
                            expression(italic(r)==15)),
#                            expression(italic(r)==20)),
               lty     =  1,
               lwd     =  2,
               col     =  c(transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[4], opacity=1)),
#                            transparentColor(COLS[5], opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

qstr  <-  qstar(U=U, YI.init=r.0.Dat$YI.t[1], x=x, s=s, h=h)

## Panel F
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.5.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.5.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.5.Dat)
        lines(YI.t[1:length(r.10.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.10.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.10.Dat)
        lines(YI.t[1:length(r.15.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.15.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.15.Dat)
#        lines(YI.t[1:length(r.20.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.20.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[5], opacity=1), data=r.20.Dat)
        # Points
        points(YI.t[length(r.0.rel.w.YI.t)] ~ length(r.0.rel.w.YI.t), pch=8, col='#252525', data=r.0.Dat)
        points(YI.t[length(r.5.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.5.rel.w.YI.t)], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.5.Dat)
        points(YI.t[length(r.10.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.10.rel.w.YI.t)], pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.10.Dat)
        points(YI.t[length(r.15.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.15.rel.w.YI.t)], pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.15.Dat)
#        points(YI.t[length(r.20.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.20.rel.w.YI.t)], pch=21, col=transparentColor(COLS[5], opacity=1), bg=transparentColor(COLS[5], opacity=0.7), data=r.20.Dat)
       # benchmark qstar
        abline(h=qstr, lty=2)
        proportionalLabel(0.9, 0.96, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'F', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # reset Colors
     COLS  <-  c('red',
                 '#252525')

## Panel I
# deleterious allele frequencies for r = 5
     # Make plot
    delMax  <-  r.10.Dat$XaOv.wt.t[1]*1.05
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(XaOv.wt.t ~ seq_along(XaOv.wt.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.10.Dat)
        lines(XaOv.del.t ~ seq_along(XaOv.del.t), lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.10.Dat)
        lines(qI.wt.t ~ seq_along(qI.wt.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.10.Dat)
#        lines(q.wt.Auto ~ autotime, lwd=1.5, lty=1, col="tomato")
        # axes
        axis(1, las=1, cex=1.2)
        axis(2, las=1, cex=1.2)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'I', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(r)==10), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        

}




# Supp FIGURE . Illustration of deterministic dynamics
#           for different dominance values
deterministicMutRateIllustration  <-  function() {

    # set plot layout
    layout.mat <- matrix(c(1:9), nrow=3, ncol=3, byrow=FALSE)
    layout     <- layout(layout.mat,respect=TRUE)

    # Make data for plotting 
    s        <-  0.01
    h        <-  0.25
    nTot     <-  10^4 
    U        <-  0.02
    u        <-  U/nTot
    x        <-  0.8

    ###########
    # h = 0.25

    # make data for plotting
    r.0.Dat  <-  makeDeterministicFigSimData(U = U, s = s, r = 0, x = x, h = h, generations = 2500)
    r.1.Dat  <-  makeDeterministicFigSimData(U = U, s = s, r = 1, x = x, h = h, generations = 2500)
#    r.2.Dat  <-  makeDeterministicFigSimData(s = s, r = 2, x = x, h = h, generations = 2500)
#    r.3.Dat  <-  makeDeterministicFigSimData(s = s, r = 3, x = x, h = h, generations = 2500)

    r.0.rel.w.YI.t  <-  r.0.Dat$rel.w.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.1.rel.w.YI.t  <-  r.1.Dat$rel.w.YI.t[round(r.1.Dat$YI.t,digits=4) > 0]
#    r.2.rel.w.YI.t  <-  r.2.Dat$rel.w.YI.t[round(r.2.Dat$YI.t,digits=4) > 0]
#    r.3.rel.w.YI.t  <-  r.3.Dat$rel.w.YI.t[round(r.3.Dat$YI.t,digits=4) > 0]

    autotime  <-  c(1:2500)
    q.wt.Auto  <-  p.t(u=10^-6, sdHom=0.01, h=h, t=autotime) 

   # Colors
#    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
#    COLS     <-  rev(colfunc(5))
#COLS  <-  colorRampPalette(wes_palette("Zissou1"), interpolate = c("linear"))(4)
COLS  <- colorRampPalette(wes_palette("Zissou1"))(4)
    # Axes limits
    relFitMin  <-  0.996
    relFitMax  <-  1.005


## Panel A
# Effect of initial del. mut. load on inversion relative fitness 
    # make plot
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(relFitMin, relFitMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        abline(h=1, lwd=1, lty=2, col=1)
#        abline(h=(1+U*x), lwd=1, lty=3, col=2)
#        proportionalLabel(0.15, 0.89, expression(italic(1+italic(Ux))), cex=1, adj=c(0.5, 0.5), xpd=NA)
        lines(r.0.rel.w.YI.t ~ seq_along(r.0.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1))
        lines(r.1.rel.w.YI.t ~ seq_along(r.1.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1))
#        lines(r.2.rel.w.YI.t ~ seq_along(r.2.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1))
#        lines(r.3.rel.w.YI.t ~ seq_along(r.3.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1))
        # Points depicting (arbitrary) extinction frequency of 0.00001
#        points(r.0.rel.w.YI.t[length(r.0.rel.w.YI.t)] ~ length(r.0.rel.w.YI.t), pch=8, col='#252525')
        points(r.1.rel.w.YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
#        points(r.2.rel.w.YI.t[length(r.2.rel.w.YI.t)] ~ length(r.2.rel.w.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
#        points(r.3.rel.w.YI.t[length(r.3.rel.w.YI.t)] ~ length(r.3.rel.w.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        # Annotate 1 + Ux approximation
        points((1 + U*x) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
#        abline(h=(1+U*x), lwd=1, lty=3, col=2)
        proportionalArrows(0.125, 0.86, 0.06, 0.86, cex=1, adj=c(0.5, 0.5), xpd=NA, log=FALSE, length=0.05)
        proportionalLabel(0.25, 0.86, expression(1+italic(Ux)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.5, 1.15, substitute(italic(h)==xx, list(xx=h)), cex=2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, expression(paste(italic(Ux/hs==1.6))), cex=2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,   expression(paste("Inversion Relative Fitness")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(r)==0),
                            expression(italic(r)==1)),
#                            expression(italic(r)==2)),
#                            expression(italic(r)==3)),
               lty     =  1,
               lwd     =  2,
               col     =  c(transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1)),
#                            transparentColor(COLS[3], opacity=1)),
#                            transparentColor(COLS[2], opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

qstr  <-  qstar(U=U, YI.init=r.0.Dat$YI.t[1], x=x, s=s, h=h)

## Panel D
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,0.1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
#        lines(YI.t[1:length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.2.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.2.Dat)
#        lines(YI.t[1:length(r.3.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.3.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.3.Dat)
        # Points
#        points(YI.t[length(r.0.rel.w.YI.t)] ~ length(r.0.rel.w.YI.t), pch=8, col='#252525', data=r.0.Dat)
        points(YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
#        points(YI.t[length(r.2.rel.w.YI.t)] ~ length(r.2.rel.w.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.2.Dat)
#        points(YI.t[length(r.3.rel.w.YI.t)] ~ length(r.3.rel.w.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.3.Dat)
        # benchmark qstar
        abline(h=qstr, lty=2)
        proportionalLabel(0.9, 0.1, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,   expression(paste("Inversion Frequency (", italic(Y[I]), ")")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        

    # reset Colors
     COLS  <-  c('red',
                 '#252525')
## Panel G
# deleterious allele frequencies for r = 1
     # Make plot
    delMax  <-  r.1.Dat$XaOv.wt.t[1]*1.05
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(XaOv.wt.t ~ seq_along(XaOv.wt.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(XaOv.del.t ~ seq_along(XaOv.del.t), lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qI.wt.t ~ seq_along(qI.wt.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
#        lines(q.wt.Auto ~ autotime, lwd=2, lty=1, col="tomato")
        # axes
        axis(1, las=1, cex=1.3)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'G', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(r)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,   expression(paste("Deleterious Allele Frequency")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4]*0.9,
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
    # make data for plotting
U  <-  0.05
    r.0.Dat  <-  makeDeterministicFigSimData(U = U, s = s, r = 0, x = x, h = h, generations = 5000)
    r.1.Dat  <-  makeDeterministicFigSimData(U = U, s = s, r = 1, x = x, h = h, generations = 5000)
    r.2.Dat  <-  makeDeterministicFigSimData(U = U, s = s, r = 2, x = x, h = h, generations = 5000)
    r.3.Dat  <-  makeDeterministicFigSimData(U = U, s = s, r = 3, x = x, h = h, generations = 5000)

    r.0.rel.w.YI.t  <-  r.0.Dat$rel.w.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.1.rel.w.YI.t  <-  r.1.Dat$rel.w.YI.t[round(r.1.Dat$YI.t,digits=4) > 0]
    r.2.rel.w.YI.t  <-  r.2.Dat$rel.w.YI.t[round(r.2.Dat$YI.t,digits=4) > 0]
    r.3.rel.w.YI.t  <-  r.3.Dat$rel.w.YI.t[round(r.3.Dat$YI.t,digits=4) > 0]

    autotime  <-  c(1:8000)
     q.wt.Auto  <-  p.t(u=10^-6, sdHom=0.01, h=h, t=autotime) 

   # Colors
#    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
#    COLS     <-  rev(colfunc(5))
    COLS     <-  colorRampPalette(wes_palette("Zissou1"))(4)

    relFitMin  <-  min(r.3.Dat$rel.w.YI.t)
    relFitMax  <-  1 + U*x*1.05

## Panel B
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(relFitMin,relFitMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        abline(h=1, lwd=1, lty=2, col=1)
        lines(r.3.rel.w.YI.t ~ seq_along(r.3.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1))
        lines(r.2.rel.w.YI.t ~ seq_along(r.2.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.1.rel.w.YI.t ~ seq_along(r.1.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r.0.rel.w.YI.t ~ seq_along(r.0.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1))
        # Points depicting (arbitrary) extinction frequency of 0.00001
#        points(r.0.rel.w.YI.t[length(r.0.rel.w.YI.t)] ~ length(r.0.rel.w.YI.t), pch=8, col='#252525')
        points(r.1.rel.w.YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r.2.rel.w.YI.t[length(r.2.rel.w.YI.t)] ~ length(r.2.rel.w.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
        points(r.3.rel.w.YI.t[length(r.3.rel.w.YI.t)] ~ length(r.3.rel.w.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        # Annotate 1 + Ux approximation
        points((1 + U*x) ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
#        abline(h=(1+U*x), lwd=1, lty=3, col=2)
#        proportionalArrows(0.125, 0.86, 0.06, 0.86, cex=1, adj=c(0.5, 0.5), xpd=NA, log=FALSE, length=0.05)
#        proportionalLabel(0.2, 0.86, expression(italic(1+italic(Ux))), cex=1, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#       proportionalLabel(0.5, 1.15, substitute(italic(h)==xx, list(xx=h)), cex=2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, expression(paste(italic(Ux/hs==4))), cex=2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(r)==0),
                            expression(italic(r)==1),
                            expression(italic(r)==2),
                            expression(italic(r)==3)),
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



qstr  <-  qstar(U=U, YI.init=r.0.Dat$YI.t[1], x=x, s=s, h=h)

## Panel E
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
#    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,0.1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(YI.t[1:length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.2.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.2.Dat)
        lines(YI.t[1:length(r.3.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.3.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.3.Dat)
        # Points
#        points(YI.t[length(r.0.rel.w.YI.t)] ~ length(r.0.rel.w.YI.t), pch=8, col='#252525', data=r.0.Dat)
        points(YI.t[length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.1.rel.w.YI.t)], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(YI.t[length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.2.rel.w.YI.t)], pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.2.Dat)
        points(YI.t[length(r.3.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.3.rel.w.YI.t)], pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.3.Dat)
        # benchmark qstar
        abline(h=qstr, lty=2)
        proportionalLabel(0.9, 0.28, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'E', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # reset Colors
     COLS  <-  c('red',
                 '#252525')
## Panel H
# deleterious allele frequencies for r = 5
     # Make plot
     delMax  <-  r.1.Dat$XaOv.wt.t[1]*1.05
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(XaOv.wt.t ~ seq_along(XaOv.wt.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(XaOv.del.t ~ seq_along(XaOv.del.t), lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qI.wt.t ~ seq_along(qI.wt.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
#        lines(q.wt.Auto ~ autotime, lwd=2, lty=1, col="tomato")
        # axes
        axis(1, las=1, cex=1.3)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'H', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(r)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        




    ###########
    # U = 0.1
    
    # make data for plotting
U  <-  0.1
#U  <- 0.004
# h  <-  0.1
    r.0.Dat   <-  makeDeterministicFigSimData(U = U, s = s, r = 0,  x = x, h = h, generations = 5000)
    r.5.Dat   <-  makeDeterministicFigSimData(U = U, s = s, r = 5,  x = x, h = h, generations = 5000)
    r.10.Dat  <-  makeDeterministicFigSimData(U = U, s = s, r = 10, x = x, h = h, generations = 5000)
    r.15.Dat  <-  makeDeterministicFigSimData(U = U, s = s, r = 15, x = x, h = h, generations = 5000)
    r.20.Dat  <-  makeDeterministicFigSimData(U = U, s = s, r = 20, x = x, h = h, generations = 5000)

    r.0.rel.w.YI.t   <-  r.0.Dat$rel.w.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.5.rel.w.YI.t   <-  r.5.Dat$rel.w.YI.t[round(r.5.Dat$YI.t,digits=4) > 0]
    r.10.rel.w.YI.t  <-  r.10.Dat$rel.w.YI.t[round(r.10.Dat$YI.t,digits=4) > 0]
    r.15.rel.w.YI.t  <-  r.15.Dat$rel.w.YI.t[round(r.15.Dat$YI.t,digits=4) > 0]
    r.20.rel.w.YI.t  <-  r.20.Dat$rel.w.YI.t[round(r.20.Dat$YI.t,digits=4) > 0]

    autotime  <-  c(1:40000)
    q.wt.Auto  <-  p.t(u=10^-6, sdHom=0.01, h=h, t=autotime) 

   # Colors
#    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
#    COLS     <-  rev(colfunc(5))
    COLS     <-  colorRampPalette(wes_palette("Zissou1"))(4)

    qHat  <-  r.0.Dat$qY.wt.t[1]
    rBar  <-  round(nTot*x*qHat)

    relFitMin  <-  0.97
    relFitMax  <-  1 + U*x*1.05

## Panel C
# Effect of initial del. mut. load on inversion relative fitness 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(relFitMin,relFitMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        abline(h=1, lwd=1, lty=2, col=1)
#        abline(h=(1+U*x), lwd=1, lty=3, col=2)
#        proportionalLabel(0.15, 0.89, expression(italic(1+italic(Ux))), cex=1, adj=c(0.5, 0.5), xpd=NA)
        lines(r.0.rel.w.YI.t ~ seq_along(r.0.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1))
        lines(r.5.rel.w.YI.t ~ seq_along(r.5.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r.10.rel.w.YI.t ~ seq_along(r.10.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.15.rel.w.YI.t ~ seq_along(r.15.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1))
#        lines(r.20.rel.w.YI.t ~ seq_along(r.20.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[5], opacity=1))
        # Points depicting (arbitrary) extinction frequency of 0.00001
#        points(r.0.rel.w.YI.t[length(r.0.rel.w.YI.t)] ~ length(r.0.rel.w.YI.t), pch=8, col='#252525')
        points(r.5.rel.w.YI.t[length(r.5.rel.w.YI.t)] ~ length(r.5.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r.10.rel.w.YI.t[length(r.10.rel.w.YI.t)] ~ length(r.10.rel.w.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
        points(r.15.rel.w.YI.t[length(r.15.rel.w.YI.t)] ~ length(r.15.rel.w.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
#        points(r.20.rel.w.YI.t[length(r.20.rel.w.YI.t)] ~ length(r.20.rel.w.YI.t), pch=21, col=transparentColor(COLS[5], opacity=1), bg=transparentColor(COLS[5], opacity=0.7))
        # Annotate 1 + Ux approximation
        points(r.0.rel.w.YI.t[1] ~ 1, pch=23, bg=transparentColor('#252525', opacity=0.7), col='#252525')
#        abline(h=(1+U*x), lwd=1, lty=3, col=2)
#        proportionalArrows(0.125, 0.86, 0.06, 0.86, cex=1, adj=c(0.5, 0.5), xpd=NA, log=FALSE, length=0.05)
#        proportionalLabel(0.2, 0.86, expression(italic(1+italic(Ux))), cex=1, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
#        proportionalLabel(0.5, 1.15, substitute(italic(h)==xx, list(xx=h)), cex=2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.15, expression(paste(italic(Ux/hs==8))), cex=2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(r)==0),
                            expression(italic(r)==5),
                            expression(italic(r)==10),
                            expression(italic(r)==15)),
#                            expression(italic(r)==20)),
               lty     =  1,
               lwd     =  2,
               col     =  c(transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[4], opacity=1)),
#                            transparentColor(COLS[5], opacity=1)),
               cex     =  1.5,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


qstr  <-  qstar(U=U, YI.init=r.0.Dat$YI.t[1], x=x, s=s, h=h)

## Panel F
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.5.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.5.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.5.Dat)
        lines(YI.t[1:length(r.10.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.10.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.10.Dat)
        lines(YI.t[1:length(r.15.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.15.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.15.Dat)
#        lines(YI.t[1:length(r.20.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.20.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[5], opacity=1), data=r.20.Dat)
        # Points
#        points(YI.t[length(r.0.rel.w.YI.t)] ~ length(r.0.rel.w.YI.t), pch=8, col='#252525', data=r.0.Dat)
        points(YI.t[length(r.5.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.5.rel.w.YI.t)], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.5.Dat)
        points(YI.t[length(r.10.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.10.rel.w.YI.t)], pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.10.Dat)
        points(YI.t[length(r.15.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.15.rel.w.YI.t)], pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.15.Dat)
#        points(YI.t[length(r.20.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.20.rel.w.YI.t)], pch=21, col=transparentColor(COLS[5], opacity=1), bg=transparentColor(COLS[5], opacity=0.7), data=r.20.Dat)
        # benchmark qstar
        abline(h=qstr, lty=2)
        proportionalLabel(0.9, 0.58, expression(paste(italic(q),"*")), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'F', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # reset Colors
     COLS  <-  c('red',
                 '#252525')

## Panel I
# deleterious allele frequencies for r = 5
     # Make plot
    delMax  <-  r.5.Dat$XaOv.wt.t[1]*1.05
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(XaOv.wt.t ~ seq_along(XaOv.wt.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.5.Dat)
        lines(XaOv.del.t ~ seq_along(XaOv.del.t), lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.5.Dat)
        lines(qI.wt.t ~ seq_along(qI.wt.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.5.Dat)
#        lines(q.wt.Auto ~ autotime, lwd=1.5, lty=1, col="tomato")
        # axes
        axis(1, las=1, cex=1.2)
        axis(2, las=1, cex=1.2)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'I', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(r)==5), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        

}

###########




#  FIXED SH
deterministicDom_Fixed_sh  <-  function() {

    # set plot layout
    layout.mat <- matrix(c(1:9), nrow=3, ncol=3, byrow=FALSE)
    layout     <- layout(layout.mat,respect=TRUE)

    # Make data for plotting 
    s        <-  0.01
    h        <-  0.25
    sh       <-  s*h
    nTot     <-  10^4 
    U        <-  0.02
    u        <-  U/nTot
    x        <-  0.2


    ###########
    # h = 0.25

    # make data for plotting
    r.0.Dat  <-  makeDeterministicFigSimData_Fixed_sh(s = s, r = 0, x = x, h = h, generations = 2500)
    r.1.Dat  <-  makeDeterministicFigSimData_Fixed_sh(s = s, r = 1, x = x, h = h, generations = 2500)
    r.2.Dat  <-  makeDeterministicFigSimData_Fixed_sh(s = s, r = 2, x = x, h = h, generations = 2500)
    r.3.Dat  <-  makeDeterministicFigSimData_Fixed_sh(s = s, r = 3, x = x, h = h, generations = 2500)

    r.0.rel.w.YI.t  <-  r.0.Dat$rel.w.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.1.rel.w.YI.t  <-  r.1.Dat$rel.w.YI.t[round(r.1.Dat$YI.t,digits=4) > 0]
    r.2.rel.w.YI.t  <-  r.2.Dat$rel.w.YI.t[round(r.2.Dat$YI.t,digits=4) > 0]
    r.3.rel.w.YI.t  <-  r.3.Dat$rel.w.YI.t[round(r.3.Dat$YI.t,digits=4) > 0]

    autotime  <-  c(1:2500)
    q.wt.Auto  <-  p.t(u=10^-6, sdHom=0.01, h=h, t=autotime) 

   # Colors
#    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
#    COLS     <-  rev(colfunc(5))
#COLS  <-  colorRampPalette(wes_palette("Zissou1"), interpolate = c("linear"))(4)
COLS  <- colorRampPalette(wes_palette("Zissou1"))(4)
    # Axes limits
    relFitMin  <-  0.992
    relFitMax  <-  1.005


## Panel A
# Effect of initial del. mut. load on inversion relative fitness 
    # make plot
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(relFitMin, relFitMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        abline(h=1, lwd=1, lty=2, col=1)
        lines(r.0.rel.w.YI.t ~ seq_along(r.0.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1))
        lines(r.1.rel.w.YI.t ~ seq_along(r.1.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r.2.rel.w.YI.t ~ seq_along(r.2.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.3.rel.w.YI.t ~ seq_along(r.3.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1))
        # Points depicting (arbitrary) extinction frequency of 0.00001
#        points(r.0.rel.w.YI.t[length(r.0.rel.w.YI.t)] ~ length(r.0.rel.w.YI.t), pch=8, col='#252525')
        points(r.1.rel.w.YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r.2.rel.w.YI.t[length(r.2.rel.w.YI.t)] ~ length(r.2.rel.w.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
        points(r.3.rel.w.YI.t[length(r.3.rel.w.YI.t)] ~ length(r.3.rel.w.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.25, 1.15, substitute(italic(h)==xx, list(xx=h)), cex=2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.75, 1.15, substitute(italic(sh)==xx, list(xx=sh)), cex=2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,   expression(paste("Inversion Relative Fitness")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(r)==0),
                            expression(italic(r)==1),
                            expression(italic(r)==2),
                            expression(italic(r)==3)),
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

## Panel D
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(YI.t[1:length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.2.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.2.Dat)
        lines(YI.t[1:length(r.3.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.3.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.3.Dat)
        # Points
#        points(YI.t[length(r.0.rel.w.YI.t)] ~ length(r.0.rel.w.YI.t), pch=8, col='#252525', data=r.0.Dat)
        points(YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(YI.t[length(r.2.rel.w.YI.t)] ~ length(r.2.rel.w.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.2.Dat)
        points(YI.t[length(r.3.rel.w.YI.t)] ~ length(r.3.rel.w.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.3.Dat)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,   expression(paste("Inversion Frequency (", italic(Y[I]), ")")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        

    # reset Colors
     COLS  <-  c('red',
                 '#252525')
## Panel G
# deleterious allele frequencies for r = 1
     # Make plot
    delMax  <-  r.1.Dat$XaOv.wt.t[1]*1.05
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(XaOv.wt.t ~ seq_along(XaOv.wt.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(XaOv.del.t ~ seq_along(XaOv.del.t), lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qI.wt.t ~ seq_along(qI.wt.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
#        lines(q.wt.Auto ~ autotime, lwd=2, lty=1, col="tomato")
        # axes
        axis(1, las=1, cex=1.3)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'G', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(r)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,   expression(paste("Deleterious Allele Frequency")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4]*0.9,
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
    # h = 0.1

    # make data for plotting
    s  <-  0.025
    h  <-  0.1
    sh <-  s*h
    r.0.Dat  <-  makeDeterministicFigSimData_Fixed_sh(s = s, r = 0, x = x, h = h, generations = 2500)
    r.1.Dat  <-  makeDeterministicFigSimData_Fixed_sh(s = s, r = 1, x = x, h = h, generations = 2500)
    r.2.Dat  <-  makeDeterministicFigSimData_Fixed_sh(s = s, r = 2, x = x, h = h, generations = 2500)
    r.3.Dat  <-  makeDeterministicFigSimData_Fixed_sh(s = s, r = 3, x = x, h = h, generations = 2500)

    r.0.rel.w.YI.t  <-  r.0.Dat$rel.w.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.1.rel.w.YI.t  <-  r.1.Dat$rel.w.YI.t[round(r.1.Dat$YI.t,digits=4) > 0]
    r.2.rel.w.YI.t  <-  r.2.Dat$rel.w.YI.t[round(r.2.Dat$YI.t,digits=4) > 0]
    r.3.rel.w.YI.t  <-  r.3.Dat$rel.w.YI.t[round(r.3.Dat$YI.t,digits=4) > 0]

    autotime  <-  c(1:8000)
     q.wt.Auto  <-  p.t(u=10^-6, sdHom=0.01, h=h, t=autotime) 

   # Colors
#    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
#    COLS     <-  rev(colfunc(5))
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
        abline(h=1, lwd=1, lty=2, col=1)
        lines(r.3.rel.w.YI.t ~ seq_along(r.3.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1))
        lines(r.2.rel.w.YI.t ~ seq_along(r.2.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.1.rel.w.YI.t ~ seq_along(r.1.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r.0.rel.w.YI.t ~ seq_along(r.0.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1))
        # Points depicting (arbitrary) extinction frequency of 0.00001
#        points(r.0.rel.w.YI.t[length(r.0.rel.w.YI.t)] ~ length(r.0.rel.w.YI.t), pch=8, col='#252525')
        points(r.1.rel.w.YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r.2.rel.w.YI.t[length(r.2.rel.w.YI.t)] ~ length(r.2.rel.w.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
        points(r.3.rel.w.YI.t[length(r.3.rel.w.YI.t)] ~ length(r.3.rel.w.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.25, 1.15, substitute(italic(h)==xx, list(xx=h)), cex=2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.75, 1.15, substitute(italic(sh)==xx, list(xx=sh)), cex=2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(r)==0),
                            expression(italic(r)==1),
                            expression(italic(r)==2),
                            expression(italic(r)==3)),
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

## Panel E
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
#    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(YI.t[1:length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.2.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.2.Dat)
        lines(YI.t[1:length(r.3.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.3.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.3.Dat)
        # Points
#        points(YI.t[length(r.0.rel.w.YI.t)] ~ length(r.0.rel.w.YI.t), pch=8, col='#252525', data=r.0.Dat)
        points(YI.t[length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.1.rel.w.YI.t)], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(YI.t[length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.2.rel.w.YI.t)], pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.2.Dat)
        points(YI.t[length(r.3.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.3.rel.w.YI.t)], pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.3.Dat)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'E', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # reset Colors
     COLS  <-  c('red',
                 '#252525')
## Panel H
# deleterious allele frequencies for r = 5
     # Make plot
     delMax  <-  r.1.Dat$XaOv.wt.t[1]*1.05
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(XaOv.wt.t ~ seq_along(XaOv.wt.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(XaOv.del.t ~ seq_along(XaOv.del.t), lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qI.wt.t ~ seq_along(qI.wt.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
#        lines(q.wt.Auto ~ autotime, lwd=2, lty=1, col="tomato")
        # axes
        axis(1, las=1, cex=1.3)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'H', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(r)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        




    ###########
    # h = 0.01
#   s        <-  0.25
#    h        <-  0.01
#    Ufactor  <-  2
#    nTot     <-  10^4 
#    U        <-  Ufactor*s
#    u        <-  U/nTot
#    x        <-  0.2
    
    # make data for plotting
    s  <-  0.25
    h  <-  0.01
    sh <-  s*h
    r.0.Dat  <-  makeDeterministicFigSimData_Fixed_sh(s = s, r = 0, x = x, h = h, generations = 2500)
    r.1.Dat  <-  makeDeterministicFigSimData_Fixed_sh(s = s, r = 1, x = x, h = h, generations = 2500)
    r.2.Dat  <-  makeDeterministicFigSimData_Fixed_sh(s = s, r = 2, x = x, h = h, generations = 2500)
    r.3.Dat  <-  makeDeterministicFigSimData_Fixed_sh(s = s, r = 3, x = x, h = h, generations = 2500)

    r.0.rel.w.YI.t  <-  r.0.Dat$rel.w.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.1.rel.w.YI.t  <-  r.1.Dat$rel.w.YI.t[round(r.1.Dat$YI.t,digits=4) > 0]
    r.2.rel.w.YI.t  <-  r.2.Dat$rel.w.YI.t[round(r.2.Dat$YI.t,digits=4) > 0]
    r.3.rel.w.YI.t  <-  r.3.Dat$rel.w.YI.t[round(r.3.Dat$YI.t,digits=4) > 0]


    autotime  <-  c(1:40000)
    q.wt.Auto  <-  p.t(u=10^-6, sdHom=0.01, h=h, t=autotime) 

   # Colors
#    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
#    COLS     <-  rev(colfunc(5))
    COLS     <-  colorRampPalette(wes_palette("Zissou1"))(4)

    qHat  <-  r.0.Dat$qY.wt.t[1]
    rBar  <-  round(nTot*x*qHat)


## Panel C
# Effect of initial del. mut. load on inversion relative fitness 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(relFitMin,relFitMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        abline(h=1, lwd=1, lty=2, col=1)
        lines(r.3.rel.w.YI.t ~ seq_along(r.3.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1))
        lines(r.2.rel.w.YI.t ~ seq_along(r.2.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.1.rel.w.YI.t ~ seq_along(r.1.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r.0.rel.w.YI.t ~ seq_along(r.0.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1))
        # Points depicting (arbitrary) extinction frequency of 0.00001
#        points(r.0.rel.w.YI.t[length(r.0.rel.w.YI.t)] ~ length(r.0.rel.w.YI.t), pch=8, col='#252525')
        points(r.1.rel.w.YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r.2.rel.w.YI.t[length(r.2.rel.w.YI.t)] ~ length(r.2.rel.w.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
        points(r.3.rel.w.YI.t[length(r.3.rel.w.YI.t)] ~ length(r.3.rel.w.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.25, 1.15, substitute(italic(h)==xx, list(xx=h)), cex=2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.75, 1.15, substitute(italic(sh)==xx, list(xx=sh)), cex=2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(r)==0),
                            expression(italic(r)==1),
                            expression(italic(r)==2),
                            expression(italic(r)==3)),
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

## Panel F
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(YI.t[1:length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.2.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.2.Dat)
        lines(YI.t[1:length(r.3.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.3.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.3.Dat)
        # Points
#        points(YI.t[length(r.0.rel.w.YI.t)] ~ length(r.0.rel.w.YI.t), pch=8, col='#252525', data=r.0.Dat)
        points(YI.t[length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.1.rel.w.YI.t)], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(YI.t[length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.2.rel.w.YI.t)], pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.2.Dat)
        points(YI.t[length(r.3.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.3.rel.w.YI.t)], pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.3.Dat)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'F', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # reset Colors
     COLS  <-  c('red',
                 '#252525')
## Panel I
# deleterious allele frequencies for r = 5
     # Make plot
    delMax  <-  r.1.Dat$XaOv.wt.t[1]*1.05
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(XaOv.wt.t ~ seq_along(XaOv.wt.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(XaOv.del.t ~ seq_along(XaOv.del.t), lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qI.wt.t ~ seq_along(qI.wt.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
#        lines(q.wt.Auto ~ autotime, lwd=1.5, lty=1, col="tomato")
        # axes
        axis(1, las=1, cex=1.2)
        axis(2, las=1, cex=1.2)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'I', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(r)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        

}








#  
deterministic_U_sh  <-  function() {

    # set plot layout
    layout.mat <- matrix(c(1:9), nrow=3, ncol=3, byrow=FALSE)
    layout     <- layout(layout.mat,respect=TRUE)

    # Make data for plotting 
    s        <-  0.01
    h        <-  0.25
    sh       <-  s*h
    nTot     <-  10^4 
    U        <-  0.02
    u        <-  U/nTot
    x        <-  0.2


    ###########
    # U  <-  0.02

    # make data for plotting
    r.0.Dat  <-  makeDeterministicFigSimData_Fixed_sh(U = U, s = s, r = 0, x = x, h = h, generations = 2500)
    r.1.Dat  <-  makeDeterministicFigSimData_Fixed_sh(U = U, s = s, r = 1, x = x, h = h, generations = 2500)
    r.2.Dat  <-  makeDeterministicFigSimData_Fixed_sh(U = U, s = s, r = 2, x = x, h = h, generations = 2500)
    r.3.Dat  <-  makeDeterministicFigSimData_Fixed_sh(U = U, s = s, r = 3, x = x, h = h, generations = 2500)

    r.0.rel.w.YI.t  <-  r.0.Dat$rel.w.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.1.rel.w.YI.t  <-  r.1.Dat$rel.w.YI.t[round(r.1.Dat$YI.t,digits=4) > 0]
    r.2.rel.w.YI.t  <-  r.2.Dat$rel.w.YI.t[round(r.2.Dat$YI.t,digits=4) > 0]
    r.3.rel.w.YI.t  <-  r.3.Dat$rel.w.YI.t[round(r.3.Dat$YI.t,digits=4) > 0]

    autotime  <-  c(1:2500)
    q.wt.Auto  <-  p.t(u=10^-6, sdHom=0.01, h=h, t=autotime) 

   # Colors
#    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
#    COLS     <-  rev(colfunc(5))
#COLS  <-  colorRampPalette(wes_palette("Zissou1"), interpolate = c("linear"))(4)
COLS  <- colorRampPalette(wes_palette("Zissou1"))(4)
    # Axes limits
    relFitMin  <-  0.992
    relFitMax  <-  1.005


## Panel A
# Effect of initial del. mut. load on inversion relative fitness 
    # make plot
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(relFitMin, relFitMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        abline(h=1, lwd=1, lty=2, col=1)
        lines(r.0.rel.w.YI.t ~ seq_along(r.0.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1))
        lines(r.1.rel.w.YI.t ~ seq_along(r.1.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r.2.rel.w.YI.t ~ seq_along(r.2.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.3.rel.w.YI.t ~ seq_along(r.3.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1))
        # Points depicting (arbitrary) extinction frequency of 0.00001
#        points(r.0.rel.w.YI.t[length(r.0.rel.w.YI.t)] ~ length(r.0.rel.w.YI.t), pch=8, col='#252525')
        points(r.1.rel.w.YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r.2.rel.w.YI.t[length(r.2.rel.w.YI.t)] ~ length(r.2.rel.w.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
        points(r.3.rel.w.YI.t[length(r.3.rel.w.YI.t)] ~ length(r.3.rel.w.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.25, 1.15, substitute(italic(h)==xx, list(xx=h)), cex=2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.75, 1.15, substitute(italic(sh)==xx, list(xx=sh)), cex=2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,   expression(paste("Inversion Relative Fitness")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(r)==0),
                            expression(italic(r)==1),
                            expression(italic(r)==2),
                            expression(italic(r)==3)),
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

## Panel D
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(YI.t[1:length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.2.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.2.Dat)
        lines(YI.t[1:length(r.3.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.3.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.3.Dat)
        # Points
#        points(YI.t[length(r.0.rel.w.YI.t)] ~ length(r.0.rel.w.YI.t), pch=8, col='#252525', data=r.0.Dat)
        points(YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(YI.t[length(r.2.rel.w.YI.t)] ~ length(r.2.rel.w.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.2.Dat)
        points(YI.t[length(r.3.rel.w.YI.t)] ~ length(r.3.rel.w.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.3.Dat)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,   expression(paste("Inversion Frequency (", italic(Y[I]), ")")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        

    # reset Colors
     COLS  <-  c('red',
                 '#252525')
## Panel G
# deleterious allele frequencies for r = 1
     # Make plot
    delMax  <-  r.1.Dat$XaOv.wt.t[1]*1.05
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(XaOv.wt.t ~ seq_along(XaOv.wt.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(XaOv.del.t ~ seq_along(XaOv.del.t), lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qI.wt.t ~ seq_along(qI.wt.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
#        lines(q.wt.Auto ~ autotime, lwd=2, lty=1, col="tomato")
        # axes
        axis(1, las=1, cex=1.3)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'G', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(r)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.3,  0.5,   expression(paste("Deleterious Allele Frequency")), cex=1.75, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4]*0.9,
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
    # U = 0.2

    

    # make data for plotting
    U  <-  0.2
    r.0.Dat  <-  makeDeterministicFigSimData_Fixed_sh(U = U, s = s, r = 0, x = x, h = h, generations = 2500)
    r.1.Dat  <-  makeDeterministicFigSimData_Fixed_sh(U = U, s = s, r = 1, x = x, h = h, generations = 2500)
    r.2.Dat  <-  makeDeterministicFigSimData_Fixed_sh(U = U, s = s, r = 2, x = x, h = h, generations = 2500)
    r.3.Dat  <-  makeDeterministicFigSimData_Fixed_sh(U = U, s = s, r = 3, x = x, h = h, generations = 2500)

    r.0.rel.w.YI.t  <-  r.0.Dat$rel.w.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.1.rel.w.YI.t  <-  r.1.Dat$rel.w.YI.t[round(r.1.Dat$YI.t,digits=4) > 0]
    r.2.rel.w.YI.t  <-  r.2.Dat$rel.w.YI.t[round(r.2.Dat$YI.t,digits=4) > 0]
    r.3.rel.w.YI.t  <-  r.3.Dat$rel.w.YI.t[round(r.3.Dat$YI.t,digits=4) > 0]

    autotime  <-  c(1:8000)
     q.wt.Auto  <-  p.t(u=10^-6, sdHom=0.01, h=h, t=autotime) 

   # Colors
    COLS     <-  colorRampPalette(wes_palette("Zissou1"))(4)

    # Analytic expectation of initial fitness for r = 0
    ylim  <-  exp(-U*x)/exp(-2*U*x)
    ylim  <-  1.05

## Panel B
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(relFitMin,ylim*1.01), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        abline(h=1, lwd=1, lty=2, col=1)
        lines(r.3.rel.w.YI.t ~ seq_along(r.3.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1))
        lines(r.2.rel.w.YI.t ~ seq_along(r.2.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.1.rel.w.YI.t ~ seq_along(r.1.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r.0.rel.w.YI.t ~ seq_along(r.0.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1))
        # Points depicting (arbitrary) extinction frequency of 0.00001
#        points(r.0.rel.w.YI.t[length(r.0.rel.w.YI.t)] ~ length(r.0.rel.w.YI.t), pch=8, col='#252525')
        points(r.1.rel.w.YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r.2.rel.w.YI.t[length(r.2.rel.w.YI.t)] ~ length(r.2.rel.w.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
        points(r.3.rel.w.YI.t[length(r.3.rel.w.YI.t)] ~ length(r.3.rel.w.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.25, 1.15, substitute(italic(h)==xx, list(xx=h)), cex=2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.75, 1.15, substitute(italic(sh)==xx, list(xx=sh)), cex=2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(r)==0),
                            expression(italic(r)==1),
                            expression(italic(r)==2),
                            expression(italic(r)==3)),
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

## Panel E
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
#    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(YI.t[1:length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.2.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.2.Dat)
        lines(YI.t[1:length(r.3.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.3.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.3.Dat)
        # Points
#        points(YI.t[length(r.0.rel.w.YI.t)] ~ length(r.0.rel.w.YI.t), pch=8, col='#252525', data=r.0.Dat)
        points(YI.t[length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.1.rel.w.YI.t)], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(YI.t[length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.2.rel.w.YI.t)], pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.2.Dat)
        points(YI.t[length(r.3.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.3.rel.w.YI.t)], pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.3.Dat)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'E', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # reset Colors
     COLS  <-  c('red',
                 '#252525')
## Panel H
# deleterious allele frequencies for r = 5
     # Make plot
     delMax  <-  r.1.Dat$XaOv.wt.t[1]*1.05
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(XaOv.wt.t ~ seq_along(XaOv.wt.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(XaOv.del.t ~ seq_along(XaOv.del.t), lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qI.wt.t ~ seq_along(qI.wt.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
#        lines(q.wt.Auto ~ autotime, lwd=2, lty=1, col="tomato")
        # axes
        axis(1, las=1, cex=1.3)
        axis(2, las=1, cex=1.3)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'H', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(r)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        




    ###########
    # h = 0.01
#   s        <-  0.25
#    h        <-  0.01
#    Ufactor  <-  2
#    nTot     <-  10^4 
#    U        <-  Ufactor*s
#    u        <-  U/nTot
#    x        <-  0.2
    
    # make data for plotting
    s  <-  0.25
    h  <-  0.01
    sh <-  s*h
    r.0.Dat  <-  makeDeterministicFigSimData_Fixed_sh(s = s, r = 0, x = x, h = h, generations = 2500)
    r.1.Dat  <-  makeDeterministicFigSimData_Fixed_sh(s = s, r = 1, x = x, h = h, generations = 2500)
    r.2.Dat  <-  makeDeterministicFigSimData_Fixed_sh(s = s, r = 2, x = x, h = h, generations = 2500)
    r.3.Dat  <-  makeDeterministicFigSimData_Fixed_sh(s = s, r = 3, x = x, h = h, generations = 2500)

    r.0.rel.w.YI.t  <-  r.0.Dat$rel.w.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.1.rel.w.YI.t  <-  r.1.Dat$rel.w.YI.t[round(r.1.Dat$YI.t,digits=4) > 0]
    r.2.rel.w.YI.t  <-  r.2.Dat$rel.w.YI.t[round(r.2.Dat$YI.t,digits=4) > 0]
    r.3.rel.w.YI.t  <-  r.3.Dat$rel.w.YI.t[round(r.3.Dat$YI.t,digits=4) > 0]


    autotime  <-  c(1:40000)
    q.wt.Auto  <-  p.t(u=10^-6, sdHom=0.01, h=h, t=autotime) 

   # Colors
#    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
#    COLS     <-  rev(colfunc(5))
    COLS     <-  colorRampPalette(wes_palette("Zissou1"))(4)

    qHat  <-  r.0.Dat$qY.wt.t[1]
    rBar  <-  round(nTot*x*qHat)


## Panel C
# Effect of initial del. mut. load on inversion relative fitness 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(relFitMin,relFitMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        abline(h=1, lwd=1, lty=2, col=1)
        lines(r.3.rel.w.YI.t ~ seq_along(r.3.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1))
        lines(r.2.rel.w.YI.t ~ seq_along(r.2.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.1.rel.w.YI.t ~ seq_along(r.1.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r.0.rel.w.YI.t ~ seq_along(r.0.rel.w.YI.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1))
        # Points depicting (arbitrary) extinction frequency of 0.00001
#        points(r.0.rel.w.YI.t[length(r.0.rel.w.YI.t)] ~ length(r.0.rel.w.YI.t), pch=8, col='#252525')
        points(r.1.rel.w.YI.t[length(r.1.rel.w.YI.t)] ~ length(r.1.rel.w.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r.2.rel.w.YI.t[length(r.2.rel.w.YI.t)] ~ length(r.2.rel.w.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
        points(r.3.rel.w.YI.t[length(r.3.rel.w.YI.t)] ~ length(r.3.rel.w.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.25, 1.15, substitute(italic(h)==xx, list(xx=h)), cex=2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.75, 1.15, substitute(italic(sh)==xx, list(xx=sh)), cex=2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(r)==0),
                            expression(italic(r)==1),
                            expression(italic(r)==2),
                            expression(italic(r)==3)),
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

## Panel F
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        lines(YI.t[1:length(r.0.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.0.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.1.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(YI.t[1:length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.2.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.2.Dat)
        lines(YI.t[1:length(r.3.rel.w.YI.t)] ~ seq_along(YI.t)[1:length(r.3.rel.w.YI.t)], lwd=2, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.3.Dat)
        # Points
#        points(YI.t[length(r.0.rel.w.YI.t)] ~ length(r.0.rel.w.YI.t), pch=8, col='#252525', data=r.0.Dat)
        points(YI.t[length(r.1.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.1.rel.w.YI.t)], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.1.Dat)
        points(YI.t[length(r.2.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.2.rel.w.YI.t)], pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.2.Dat)
        points(YI.t[length(r.3.rel.w.YI.t)] ~ seq_along(YI.t)[length(r.3.rel.w.YI.t)], pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.3.Dat)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'F', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # reset Colors
     COLS  <-  c('red',
                 '#252525')
## Panel I
# deleterious allele frequencies for r = 5
     # Make plot
    delMax  <-  r.1.Dat$XaOv.wt.t[1]*1.05
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,delMax), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(XaOv.wt.t ~ seq_along(XaOv.wt.t), lwd=2, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(XaOv.del.t ~ seq_along(XaOv.del.t), lwd=2, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qI.wt.t ~ seq_along(qI.wt.t), lwd=2, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
#        lines(q.wt.Auto ~ autotime, lwd=1.5, lty=1, col="tomato")
        # axes
        axis(1, las=1, cex=1.2)
        axis(2, las=1, cex=1.2)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'I', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(italic(r)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.75, adj=c(0.5, 0.5), xpd=NA)        

}









#' Deterministic dynamics of inversions expanding the SDR on a Y chromosome:
#' Illustration of: (A) inversion relative fitness
#'                  (B) Inversion frequency dynamics, and 
#'                  (C) del. allele frequency dynamics.
deterministicFig  <-  function() {

	# make data for plotting
	r0Dat  <-  makeDeterministicFigSimData(r = 0, x=0.5, h=0.01)
	r1Dat  <-  makeDeterministicFigSimData(r = 1, x=0.5, h=0.01)
	r2Dat  <-  makeDeterministicFigSimData(r = 2, x=0.5, h=0.01)
	r3Dat  <-  makeDeterministicFigSimData(r = 3, x=0.5, h=0.01)
	r4Dat  <-  makeDeterministicFigSimData(r = 4, x=0.5, h=0.01)

	# Colors
	colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
	COLS     <-  colfunc(5)

    # set plot layout
    layout.mat <- matrix(c(1:3), nrow=3, ncol=1, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

## Panel A
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r0Dat)), ylim = c(0.99,1.021), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        r0.wbar.YI.t  <-  r0Dat$wbar.YI.t[round(r0Dat$YI.t,digits=4) < 1]
        r1.wbar.YI.t  <-  r1Dat$wbar.YI.t[round(r1Dat$YI.t,digits=4) > 0]
        r2.wbar.YI.t  <-  r2Dat$wbar.YI.t[round(r2Dat$YI.t,digits=4) > 0]
        r3.wbar.YI.t  <-  r3Dat$wbar.YI.t[round(r3Dat$YI.t,digits=4) > 0]
        r4.wbar.YI.t  <-  r4Dat$wbar.YI.t[round(r4Dat$YI.t,digits=4) > 0]
        lines(r4.wbar.YI.t ~ seq_along(r4.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1))
        lines(r3.wbar.YI.t ~ seq_along(r3.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1))
        lines(r2.wbar.YI.t ~ seq_along(r2.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r1.wbar.YI.t ~ seq_along(r1.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r0.wbar.YI.t ~ seq_along(r0.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1))
		abline(h=1, lwd=1, lty=2, col=2)
		# Points depicting (arbitrary) extinction frequency of 0.00001
		points(r4.wbar.YI.t[length(r4.wbar.YI.t)] ~ length(r4.wbar.YI.t), pch=21, col=transparentColor(COLS[5], opacity=1), bg=transparentColor(COLS[5], opacity=0.7))
		points(r3.wbar.YI.t[length(r3.wbar.YI.t)] ~ length(r3.wbar.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        points(r2.wbar.YI.t[length(r2.wbar.YI.t)] ~ length(r2.wbar.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
#        points(r1.wbar.YI.t[length(r1.wbar.YI.t)] ~ length(r1.wbar.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r0.wbar.YI.t[length(r0.wbar.YI.t)] ~ length(r0.wbar.YI.t), pch=8, col='tomato')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.225,  0.5,   expression(paste("Inversion Relative Fitness")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(paste(italic(r), " = ", 0)),
                            expression(paste(italic(r), " = ", 1)),
                            expression(paste(italic(r), " = ", 2)),
                            expression(paste(italic(r), " = ", 3)),
                            expression(paste(italic(r), " = ", 4))),
               lty     =  1,
               lwd     =  1.5,
               col     =  c(transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[4], opacity=1),
                            transparentColor(COLS[5], opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

## Panel B
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
#    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r0Dat)), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        r0.YI.t  <-  r0Dat$YI.t[round(r0Dat$YI.t,digits=4) < 1]
        r1.YI.t  <-  r1Dat$YI.t[round(r1Dat$YI.t,digits=4) > 0]
        r2.YI.t  <-  r2Dat$YI.t[round(r2Dat$YI.t,digits=4) > 0]
        r3.YI.t  <-  r3Dat$YI.t[round(r3Dat$YI.t,digits=4) > 0]
        r4.YI.t  <-  r4Dat$YI.t[round(r4Dat$YI.t,digits=4) > 0]
        lines(r4.YI.t ~ seq_along(r4.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1))
        lines(r3.YI.t ~ seq_along(r3.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1))
        lines(r2.YI.t ~ seq_along(r2.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r1.YI.t ~ seq_along(r1.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r0.YI.t ~ seq_along(r0.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1))
        # Points depicting (arbitrary) extinction frequency of 0.00001
        points(r4.YI.t[length(r4.YI.t)] ~ length(r4.YI.t), pch=21, col=transparentColor(COLS[5], opacity=1), bg=transparentColor(COLS[5], opacity=0.7))
        points(r3.YI.t[length(r3.YI.t)] ~ length(r3.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        points(r2.YI.t[length(r2.YI.t)] ~ length(r2.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
#        points(r1.YI.t[length(r1.YI.t)] ~ length(r1.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r0.YI.t[length(r0.YI.t)] ~ length(r0.YI.t), pch=8, col='tomato')

        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.225,  0.5,   expression(paste("Inversion Frequency (", italic(Y[I]), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(paste(italic(r), " = ", 0)),
                            expression(paste(italic(r), " = ", 1)),
                            expression(paste(italic(r), " = ", 2)),
                            expression(paste(italic(r), " = ", 3)),
                            expression(paste(italic(r), " = ", 4))),
               lty     =  1,
               lwd     =  1.5,
               col     =  c(transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[4], opacity=1),
                            transparentColor(COLS[5], opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

    # reset Colors
     COLS  <-  c('red',
      			 '#252525')
## Panel C
# deleterious allele frequencies for r = 3
     # Make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r2Dat)), ylim = c(0,0.002), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(XaOv.wt.t ~ seq_along(XaOv.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r1Dat)
        lines(XaOv.del.t ~ seq_along(XaOv.del.t), lwd=1.5, lty=2, col=transparentColor(COLS[2], opacity=1), data=r1Dat)
        lines(qI.wt.t ~ seq_along(qI.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r1Dat)
#        lines(XaSp.wt.t ~ seq_along(XaSp.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1), data=r1Dat)
#        lines(qY.del.t ~ seq_along(qY.del.t), lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1), data=r1Dat)
#        lines(XaSp.del.t ~ seq_along(XaSp.del.t), lwd=1.5, lty=1, col=transparentColor(COLS[7], opacity=1), data=r1Dat)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(paste(italic(r), " = ", 1)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.225,  0.5,   expression(paste("Deleterious Allele Frequency")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.2,  expression(paste("Generations")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4]*0.25,
			   legend = c(expression(italic(q)[italic(X[f])]^italic(wt)),
                          expression(italic(q)[italic(Y[I])]^italic(wt)), 
						  expression(italic(q)[italic(X[f])]^italic(del))),
               lty     =  c(1,1,2),
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

}

#' Deterministic dynamics of inversions expanding the SDR on a Y chromosome:
#' Illustration of: (A) inversion relative fitness
#'                  (B) Inversion frequency dynamics, and 
#'                  (C) del. allele frequency dynamics.
deterministicFigLessRecessive  <-  function() {

    # make data for plotting
    r0Dat  <-  makeDeterministicFigSimData(r = 0, x=0.5, h=0.25)
    r1Dat  <-  makeDeterministicFigSimData(r = 1, x=0.5, h=0.25)
    r2Dat  <-  makeDeterministicFigSimData(r = 2, x=0.5, h=0.25)
    r3Dat  <-  makeDeterministicFigSimData(r = 3, x=0.5, h=0.25)
    r4Dat  <-  makeDeterministicFigSimData(r = 4, x=0.5, h=0.25)

    # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  colfunc(5)

    # set plot layout
    layout.mat <- matrix(c(1:3), nrow=3, ncol=1, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

## Panel A
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,2500), ylim = c(0.99,1.021), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        r0.wbar.YI.t  <-  r0Dat$wbar.YI.t[round(r0Dat$YI.t,digits=4) < 1]
        r1.wbar.YI.t  <-  r1Dat$wbar.YI.t[round(r1Dat$YI.t,digits=4) > 0]
        r2.wbar.YI.t  <-  r2Dat$wbar.YI.t[round(r2Dat$YI.t,digits=4) > 0]
        r3.wbar.YI.t  <-  r3Dat$wbar.YI.t[round(r3Dat$YI.t,digits=4) > 0]
        r4.wbar.YI.t  <-  r4Dat$wbar.YI.t[round(r4Dat$YI.t,digits=4) > 0]
        lines(r4.wbar.YI.t ~ seq_along(r4.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1))
        lines(r3.wbar.YI.t ~ seq_along(r3.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1))
        lines(r2.wbar.YI.t ~ seq_along(r2.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r1.wbar.YI.t ~ seq_along(r1.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r0.wbar.YI.t ~ seq_along(r0.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1))
        abline(h=1, lwd=1, lty=2, col=2)
        # Points depicting (arbitrary) extinction frequency of 0.00001
        points(r4.wbar.YI.t[length(r4.wbar.YI.t)] ~ length(r4.wbar.YI.t), pch=21, col=transparentColor(COLS[5], opacity=1), bg=transparentColor(COLS[5], opacity=0.7))
        points(r3.wbar.YI.t[length(r3.wbar.YI.t)] ~ length(r3.wbar.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        points(r2.wbar.YI.t[length(r2.wbar.YI.t)] ~ length(r2.wbar.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
        points(r1.wbar.YI.t[length(r1.wbar.YI.t)] ~ length(r1.wbar.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r0.wbar.YI.t[length(r0.wbar.YI.t)] ~ length(r0.wbar.YI.t), pch=8, col='tomato')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.225,  0.5,   expression(paste("Inversion Relative Fitness")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(paste(italic(r), " = ", 0)),
                            expression(paste(italic(r), " = ", 1)),
                            expression(paste(italic(r), " = ", 2)),
                            expression(paste(italic(r), " = ", 3)),
                            expression(paste(italic(r), " = ", 4))),
               lty     =  1,
               lwd     =  1.5,
               col     =  c(transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[4], opacity=1),
                            transparentColor(COLS[5], opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

## Panel B
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
#    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,2500), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        r0.YI.t  <-  r0Dat$YI.t[round(r0Dat$YI.t,digits=4) < 1]
        r1.YI.t  <-  r1Dat$YI.t[round(r1Dat$YI.t,digits=4) > 0]
        r2.YI.t  <-  r2Dat$YI.t[round(r2Dat$YI.t,digits=4) > 0]
        r3.YI.t  <-  r3Dat$YI.t[round(r3Dat$YI.t,digits=4) > 0]
        r4.YI.t  <-  r4Dat$YI.t[round(r4Dat$YI.t,digits=4) > 0]
        lines(r4.YI.t ~ seq_along(r4.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1))
        lines(r3.YI.t ~ seq_along(r3.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1))
        lines(r2.YI.t ~ seq_along(r2.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r1.YI.t ~ seq_along(r1.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r0.YI.t ~ seq_along(r0.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1))
        # Points depicting (arbitrary) extinction frequency of 0.00001
        points(r4.YI.t[length(r4.YI.t)] ~ length(r4.YI.t), pch=21, col=transparentColor(COLS[5], opacity=1), bg=transparentColor(COLS[5], opacity=0.7))
        points(r3.YI.t[length(r3.YI.t)] ~ length(r3.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        points(r2.YI.t[length(r2.YI.t)] ~ length(r2.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
        points(r1.YI.t[length(r1.YI.t)] ~ length(r1.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r0.YI.t[length(r0.YI.t)] ~ length(r0.YI.t), pch=8, col='tomato')

        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.225,  0.5,   expression(paste("Inversion Frequency (", italic(Y[I]), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(paste(italic(r), " = ", 0)),
                            expression(paste(italic(r), " = ", 1)),
                            expression(paste(italic(r), " = ", 2)),
                            expression(paste(italic(r), " = ", 3)),
                            expression(paste(italic(r), " = ", 4))),
               lty     =  1,
               lwd     =  1.5,
               col     =  c(transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[4], opacity=1),
                            transparentColor(COLS[5], opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

    # reset Colors
     COLS  <-  c('red',
                 '#252525')
## Panel C
# deleterious allele frequencies for r = 3
     # Make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,2500), ylim = c(0,0.0009), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(XaOv.wt.t ~ seq_along(XaOv.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r1Dat)
        lines(XaOv.del.t ~ seq_along(XaOv.del.t), lwd=1.5, lty=2, col=transparentColor(COLS[2], opacity=1), data=r1Dat)
        lines(qI.wt.t ~ seq_along(qI.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r1Dat)
#        lines(XaSp.wt.t ~ seq_along(XaSp.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1), data=r1Dat)
#        lines(qY.del.t ~ seq_along(qY.del.t), lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1), data=r1Dat)
#        lines(XaSp.del.t ~ seq_along(XaSp.del.t), lwd=1.5, lty=1, col=transparentColor(COLS[7], opacity=1), data=r1Dat)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(paste(italic(r), " = ", 1)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.225,  0.5,   expression(paste("Deleterious Allele Frequency")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.2,  expression(paste("Generations")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4]*0.25,
               legend = c(expression(italic(q)[italic(X[f])]^italic(wt)),
                          expression(italic(q)[italic(Y[I])]^italic(wt)), 
                          expression(italic(q)[italic(X[f])]^italic(del))),
               lty     =  c(1,1,2),
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

}


#' Deterministic dynamics of inversions expanding the SDR on a Y chromosome:
#' Illustration of: (A) inversion relative fitness
#'                  (B) Inversion frequency dynamics, and 
#'                  (C) del. allele frequency dynamics.
deterministicFigStrongRecessive  <-  function() {

    # make data for plotting
    r0Dat  <-  makeDeterministicFigSimData(r = 0, h=0.05, x=0.5, generations=20000)
    r1Dat  <-  makeDeterministicFigSimData(r = 1, h=0.05, x=0.5, generations=20000)
    r2Dat  <-  makeDeterministicFigSimData(r = 2, h=0.05, x=0.5, generations=20000)
    r3Dat  <-  makeDeterministicFigSimData(r = 3, h=0.05, x=0.5, generations=20000)
    r4Dat  <-  makeDeterministicFigSimData(r = 4, h=0.05, x=0.5, generations=20000)

    # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  colfunc(5)

    # set plot layout
    layout.mat <- matrix(c(1:3), nrow=3, ncol=1, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

## Panel A
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r0Dat)), ylim = c(0.995,1.021), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        r0.wbar.YI.t  <-  r0Dat$wbar.YI.t[round(r0Dat$YI.t,digits=4) < 1]
        r1.wbar.YI.t  <-  r1Dat$wbar.YI.t[round(r1Dat$YI.t,digits=4) > 0]
        r2.wbar.YI.t  <-  r2Dat$wbar.YI.t[round(r2Dat$YI.t,digits=4) > 0]
        r3.wbar.YI.t  <-  r3Dat$wbar.YI.t[round(r3Dat$YI.t,digits=4) > 0]
        r4.wbar.YI.t  <-  r4Dat$wbar.YI.t[round(r4Dat$YI.t,digits=4) > 0]
        lines(r4.wbar.YI.t ~ seq_along(r4.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1))
        lines(r3.wbar.YI.t ~ seq_along(r3.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1))
        lines(r2.wbar.YI.t ~ seq_along(r2.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r1.wbar.YI.t ~ seq_along(r1.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r0.wbar.YI.t ~ seq_along(r0.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1))
        abline(h=1, lwd=1, lty=2, col=2)
        # Points depicting (arbitrary) extinction frequency of 0.00001
        points(r0.wbar.YI.t[length(r0.wbar.YI.t)] ~ length(r0.wbar.YI.t), pch=8, col='tomato')
        points(r4.wbar.YI.t[length(r4.wbar.YI.t)] ~ length(r4.wbar.YI.t), pch=21, col=transparentColor(COLS[5], opacity=1), bg=transparentColor(COLS[5], opacity=0.7))
        points(r3.wbar.YI.t[length(r3.wbar.YI.t)] ~ length(r3.wbar.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
#        points(r2.wbar.YI.t[length(r2.wbar.YI.t)] ~ length(r2.wbar.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
#        points(r1.wbar.YI.t[length(r1.wbar.YI.t)] ~ length(r1.wbar.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
#        points(r0.wbar.YI.t[length(r0.wbar.YI.t)] ~ length(r0.wbar.YI.t), pch=8, col='tomato')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.225,  0.5,   expression(paste("Inversion Relative Fitness")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(paste(italic(r), " = ", 0)),
                            expression(paste(italic(r), " = ", 1)),
                            expression(paste(italic(r), " = ", 2)),
                            expression(paste(italic(r), " = ", 3)),
                            expression(paste(italic(r), " = ", 4))),
               lty     =  1,
               lwd     =  1.5,
               col     =  c(transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[4], opacity=1),
                            transparentColor(COLS[5], opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

## Panel B
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
#    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r0Dat)), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        r0.YI.t  <-  r0Dat$YI.t[round(r0Dat$YI.t,digits=4) < 1]
        r1.YI.t  <-  r1Dat$YI.t[round(r1Dat$YI.t,digits=4) > 0]
        r2.YI.t  <-  r2Dat$YI.t[round(r2Dat$YI.t,digits=4) > 0]
        r3.YI.t  <-  r3Dat$YI.t[round(r3Dat$YI.t,digits=4) > 0]
        r4.YI.t  <-  r4Dat$YI.t[round(r4Dat$YI.t,digits=4) > 0]
        lines(r4.YI.t ~ seq_along(r4.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1))
        lines(r3.YI.t ~ seq_along(r3.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1))
        lines(r2.YI.t ~ seq_along(r2.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r1.YI.t ~ seq_along(r1.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r0.YI.t ~ seq_along(r0.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1))
        # Points depicting (arbitrary) extinction frequency of 0.00001
        points(r4.YI.t[length(r4.YI.t)] ~ length(r4.YI.t), pch=21, col=transparentColor(COLS[5], opacity=1), bg=transparentColor(COLS[5], opacity=0.7))
        points(r3.YI.t[length(r3.YI.t)] ~ length(r3.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        points(r0.YI.t[length(r0.YI.t)] ~ length(r0.YI.t), pch=8, col='tomato')
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.225,  0.5,   expression(paste("Inversion Frequency (", italic(Y[I]), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4]*0.7,
               legend  =  c(
                            expression(paste(italic(r), " = ", 0)),
                            expression(paste(italic(r), " = ", 1)),
                            expression(paste(italic(r), " = ", 2)),
                            expression(paste(italic(r), " = ", 3)),
                            expression(paste(italic(r), " = ", 4))),
               lty     =  1,
               lwd     =  1.5,
               col     =  c(transparentColor(COLS[1], opacity=1),
                            transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[4], opacity=1),
                            transparentColor(COLS[5], opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

    # reset Colors
     COLS  <-  c('red',
                 '#252525')
## Panel C
# deleterious allele frequencies for r = 3
     # Make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r4Dat)), ylim = c(0,0.004), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(XaOv.wt.t ~ seq_along(XaOv.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r4Dat)
        lines(XaOv.del.t ~ seq_along(XaOv.del.t), lwd=1.5, lty=2, col=transparentColor(COLS[2], opacity=1), data=r4Dat)
        lines(qI.wt.t ~ seq_along(qI.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r4Dat)
#        lines(XaSp.wt.t ~ seq_along(XaSp.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1), data=r4Dat)
#        lines(qY.del.t ~ seq_along(qY.del.t), lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1), data=r4Dat)
#        lines(XaSp.del.t ~ seq_along(XaSp.del.t), lwd=1.5, lty=1, col=transparentColor(COLS[7], opacity=1), data=r4Dat)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(paste(italic(r), " = ", 4)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.225,  0.5,   expression(paste("Deleterious Allele Frequency")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.2,  expression(paste("Generations")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4]*0.9,
               legend = c(expression(italic(q)[italic(X[f])]^italic(wt)),
                          expression(italic(q)[italic(Y[I])]^italic(wt)), 
                          expression(italic(q)[italic(X[f])]^italic(del))),
               lty     =  c(1,1,2),
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

}


#' Fixation probability of different size inversions 
#' expanding the SDR on a Y chromosome: Pr(fix | x) 
#' Illustration of Pr(fix | x) for different 
#' population size and U/s ratios.
PrFixFig  <-  function() {

	# Import data for plotting
	dat  <-  read.csv(file = './data/PrFixFig_h0.25_s0.01_N1k_deterministic_q.csv', header=TRUE)
	dat$PrFixNe  <-  dat$PrFix*(dat$N/2)
	uniqueU  <-  unique(dat$Ufac)
	datN1k   <-  dat[dat$N == 1000,]
		datN1kUf1  <-  datN1k[datN1k$Ufac == uniqueU[1],]
		datN1kUf2  <-  datN1k[datN1k$Ufac == uniqueU[2],]
		datN1kUf3  <-  datN1k[datN1k$Ufac == uniqueU[3],]

    dat  <-  read.csv(file = './data/PrFixFig_h0.25_s0.01_N10k_deterministic_q.csv', header=TRUE)
    dat$PrFixNe  <-  dat$PrFix*(dat$N/2)
    uniqueU  <-  unique(dat$Ufac)
    datN10k   <-  dat[dat$N == 10000,]
        datN10kUf1  <-  datN10k[datN10k$Ufac == uniqueU[1],]
        datN10kUf2  <-  datN10k[datN10k$Ufac == uniqueU[2],]
        datN10kUf3  <-  datN10k[datN10k$Ufac == uniqueU[3],]


    Adat  <-  read.csv(file = './data/PrFixAuto3Fig_h0.25_s0.01.csv', header=TRUE)
    Adat$PrFixNe  <-  Adat$PrFix*(2*Adat$N)
    Adat$PrFixAdj  <-  Adat$PrFix
    Adat$PrFixAdj[Adat$PrFixAdj < (1/(10^6))]  <-  (1/(10^6))
    uniqueU    <-  unique(Adat$Ufac)
    AdatN1k    <-  Adat[Adat$N == 1000,]
        AdatN1kUf1  <-  AdatN1k[AdatN1k$Ufac == uniqueU[1],]
        AdatN1kUf2  <-  AdatN1k[AdatN1k$Ufac == uniqueU[2],]
        AdatN1kUf3  <-  AdatN1k[AdatN1k$Ufac == uniqueU[3],]
    AdatN10k   <-  Adat[Adat$N == 10000,]
        AdatN10kUf1  <-  AdatN10k[AdatN10k$Ufac == uniqueU[1],]
        AdatN10kUf2  <-  AdatN10k[AdatN10k$Ufac == uniqueU[2],]
        AdatN10kUf3  <-  AdatN10k[AdatN10k$Ufac == uniqueU[3],]

#    IBMDat  <-  read.csv(file = './data/PrFixIBM_h0.1_s0.01_N1000_Ufac2__xTest.csv', header=TRUE)

	# Colors
	colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
	COLS     <-  colfunc(3)

    # set plot layout
    layout.mat <- matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

## Panel A
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(4,4,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c((1/(10^6)), (0.75/(2*10^2))), log='y', ylab='', xlab='', cex.lab=1.2)
#     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80', log='y')
        box()
        # Points
        abline(h=c((1/(2*10^4)), (1/(2*10^3))), lwd=1, lty=2, col=1)
#        abline(h=1, lwd=1, lty=1, col=1)
        points(PrFixAdj ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[3], opacity=0.6), data=AdatN1kUf1)
        points(PrFixAdj ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[3], opacity=0.6), data=AdatN1kUf2)
        points(PrFixAdj ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[3], opacity=0.6), data=AdatN1kUf3)
        points(PrFixAdj ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.6), data=AdatN10kUf1)
        points(PrFixAdj ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.6), data=AdatN10kUf2)
        points(PrFixAdj ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.6), data=AdatN10kUf3)
        # axes
        axis(1, las=1)
#        axis(2, las=1)
        axis(2,las=1, at=c((1/(10^6)), (1/(2*10^5)), (1/(2*10^4)), (1/(2*10^3))),
            labels=c(0, expression(1/(2%*%10^5)), expression(1/(2%*%10^4)), expression(1/(2%*%10^3))))
        proportionalLabel(-0.0175, 0.16, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(-0.0175, 0.14, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(0.04, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.5, 1.1, 'Autosomal', cex=1.25, adj=c(0.5, 0.5), xpd=NA, log='y')
        # Plot labels etc.
        proportionalLabel(-0.45, 0.5, expression(Fixation~Probability), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90, log='y')
        proportionalLabel( 0.5,  -0.225,  expression(paste("Inversion size (", italic(x),")")), log='y', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
               x       =  0.6,
               y       =  0.006,
               legend  =  c(
                            expression(paste(italic(N), " = ", 10^3)),
                            expression(paste(italic(N), " = ", 10^4))),
               pch     =  22,
               col     =  transparentColor(COLS[1], opacity=0),
               pt.bg   =  c(transparentColor(COLS[3], opacity=0.6),
                            transparentColor(COLS[1], opacity=0.6)),
               pt.cex  = 1.5,
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
        legend(
               x       =  1.02,
               y       =  0.006,
               legend  =  c(
                            expression(paste(italic(U/hs), " = ", 2)),
                            expression(paste(italic(U/hs), " = ", 5)),
                            expression(paste(italic(U/hs), " = ", 10))),
               pch     =  c(21,22,24),
               col     =  transparentColor(COLS[1], opacity=1),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


## Panel B
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
#    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c((10^-6), 2/(10^2)), log='y', ylab='', xlab='', cex.lab=1.2)
#     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c(0, 20), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80', log='y')
#        plotGrid(lineCol='grey80')
        box()
        # Points
        abline(h=(2/1000), lwd=1, lty=2, col=1)
#        abline(h=1, lwd=1, lty=1, col=1)
        points(PrFix ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[3], opacity=0.6), data=datN1kUf1)
        points(PrFix ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[3], opacity=0.6), data=datN1kUf2)
        points(PrFix ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[3], opacity=0.6), data=datN1kUf3)
        abline(h=(2/10000), lwd=1, lty=2, col=1)
        points(PrFix ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.6), data=datN10kUf1)
        points(PrFix ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.6), data=datN10kUf2)
        points(PrFix ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.6), data=datN10kUf3)

#        points(PrFix ~ x, pch=21, col=transparentColor('tomato', opacity=1), bg=transparentColor('tomato', opacity=0.6), data=IBMDat)
        # axes
        axis(1, las=1)
        axis(2, las=1, at=c(2/10^6, 2/10^5, 2/10^4, 2/10^3, 2/10^2, 2/10^1), 
        		labels=c(expression(2/10^6), expression(2/10^5), expression(2/10^4), expression(2/10^3), expression(2/10^2), expression(2/10^1)))
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.5, 1.1, 'Expand SLR', cex=1.25, adj=c(0.5, 0.5), xpd=NA, log='y')
#        proportionalLabel(-0.25,  0.5,   expression(paste("Fixation Probability")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90, log='y')        
        proportionalLabel( 0.5,  -0.225,  expression(paste("Inversion size (", italic(x),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')

}




#' Fixation probability of different size inversions 
#' expanding the SDR on a Y chromosome: Pr(fix | x) 
#' Illustration of Pr(fix | x) for different 
#' population size and U/s ratios.
PrFixFig4Panel  <-  function() {

    # Import data for plotting

    ## Autosomal Inversion Data
    Adat  <-  read.csv(file = './data/PrFixAuto3Fig_h0.25_s0.01.csv', header=TRUE)
    Adat$PrFixNe  <-  Adat$PrFix*(2*Adat$N)
    Adat$PrFixAdj  <-  Adat$PrFix
#    Adat$PrFixAdj[Adat$PrFixAdj < (1/(2*10^7))]  <-  (1/(10^7))
    Adat$PrFixAdj[Adat$PrFixAdj == 0]  <-  (1/(2*10^7))
    uniqueU    <-  unique(Adat$Ufac)
    AdatN10k   <-  Adat[Adat$N == 10000,]
        AdatN10kUf1  <-  AdatN10k[AdatN10k$Ufac == uniqueU[1],]
        AdatN10kUf2  <-  AdatN10k[AdatN10k$Ufac == uniqueU[2],]
        AdatN10kUf3  <-  AdatN10k[AdatN10k$Ufac == uniqueU[3],]

    Adat  <-  read.csv(file = './data/PrFixAutoFig_h0.25_s0.01_N100k.csv', header=TRUE)
    Adat$PrFixNe  <-  Adat$PrFix*(Adat$N/2)
    uniqueU  <-  unique(Adat$Ufac)
    Adat$PrFixAdj  <-  Adat$PrFix
    Adat$PrFixAdj[Adat$PrFixAdj < (1/(2*10^7))]  <-  (1/(2*10^7))
    AdatN100k   <-  Adat[Adat$N == 100000,]
        AdatN100kUf1  <-  AdatN100k[AdatN100k$Ufac == uniqueU[1],]
        AdatN100kUf2  <-  AdatN100k[AdatN100k$Ufac == uniqueU[2],]
        AdatN100kUf3  <-  AdatN100k[AdatN100k$Ufac == uniqueU[3],]

    ## SLR Expanding Inversion data
    ## N = 100k data
    datN100kUs1  <-  read.csv(file = './data/PrFixFig_h0.25_s0.01_N100k_Us2_deterministic_q.csv', header=TRUE)
    datN100kUs2  <-  read.csv(file = './data/PrFixFig_h0.25_s0.01_N100k_Us5_deterministic_q.csv', header=TRUE)
    datN100kUs3  <-  read.csv(file = './data/PrFixFig_h0.25_s0.01_N100k_Us10_deterministic_q.csv', header=TRUE)
    datN100kUs1$PrFixNe  <-  datN100kUs1$PrFix*(datN100kUs1$N/2)
    datN100kUs2$PrFixNe  <-  datN100kUs2$PrFix*(datN100kUs2$N/2)
    datN100kUs3$PrFixNe  <-  datN100kUs3$PrFix*(datN100kUs3$N/2)
    datN100kUs1$PrFix[datN100kUs1$PrFix < (2/(10^8))]  <-  (2/(10^8))
    datN100kUs2$PrFix[datN100kUs2$PrFix < (2/(10^8))]  <-  (2/(10^8))
    datN100kUs3$PrFix[datN100kUs3$PrFix < (2/(10^8))]  <-  (2/(10^8))
    rdatN100kUs1  <-  read.csv(file = './data/PrFixFig_rFixedInv_h0.25_s0.01_N100k_Us2_deterministic_q.csv', header=TRUE)
    rdatN100kUs2  <-  read.csv(file = './data/PrFixFig_rFixedInv_h0.25_s0.01_N100k_Us5_deterministic_q.csv', header=TRUE)
    rdatN100kUs3  <-  read.csv(file = './data/PrFixFig_rFixedInv_h0.25_s0.01_N100k_Us10_deterministic_q.csv', header=TRUE)
    ## N = 1 mil data
    datN1milUs1  <-  read.csv(file = './data/PrFixFig_h0.25_s0.01_N1mil_Us2_deterministic_q.csv', header=TRUE)
    datN1milUs2  <-  read.csv(file = './data/PrFixFig_h0.25_s0.01_N1mil_Us5_deterministic_q.csv', header=TRUE)
    datN1milUs3  <-  read.csv(file = './data/PrFixFig_h0.25_s0.01_N1mil_Us10_deterministic_q.csv', header=TRUE)
    datN1milUs1$PrFixNe  <-  datN1milUs1$PrFix*(datN1milUs1$N/2)
    datN1milUs2$PrFixNe  <-  datN1milUs2$PrFix*(datN1milUs2$N/2)
    datN1milUs3$PrFixNe  <-  datN1milUs3$PrFix*(datN1milUs3$N/2)
    datN1milUs1$PrFix[datN1milUs1$PrFix < (2/(10^8))]  <-  (2/(10^8))
    datN1milUs2$PrFix[datN1milUs2$PrFix < (2/(10^8))]  <-  (2/(10^8))
    datN1milUs3$PrFix[datN1milUs3$PrFix < (2/(10^8))]  <-  (2/(10^8))
    rdatN1milUs1  <-  read.csv(file = './data/PrFixFig_rFixedInv_h0.25_s0.01_N1mil_Us2_deterministic_q.csv', header=TRUE)
    rdatN1milUs2  <-  read.csv(file = './data/PrFixFig_rFixedInv_h0.25_s0.01_N1mil_Us5_deterministic_q.csv', header=TRUE)
    rdatN1milUs3  <-  read.csv(file = './data/PrFixFig_rFixedInv_h0.25_s0.01_N1mil_Us10_deterministic_q.csv', header=TRUE)

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
#     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
         usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80', log='y')
        box()
        # Points
        abline(h=(1/(2*10^4)), lwd=1, lty=2, col=1)
#        abline(h=1, lwd=1, lty=1, col=1)
        points(PrFixAdj ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=AdatN10kUf1)
        points(PrFixAdj ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=AdatN10kUf2)
        points(PrFixAdj ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=AdatN10kUf3)
        # axes
        axis(1, las=1, labels=NA)
#        axis(2,las=1, at=c((1/(10^6)), (1/(2*10^5)), (1/(2*10^4)), (1/(2*10^3)), (1/(2*10^2))),
#            labels=NA)
         axis(2, las=1, at=c((1/(2*10^7)), (1/(2*10^6)), (1/(2*10^5)), (1/(2*10^4)), (1/(2*10^3))), 
                labels=c(0, expression(1/(2%*%10^6)), expression(1/(2%*%10^5)), expression(1/(2%*%10^4)), expression(1/(2%*%10^3))))
        proportionalLabel(1.15, 1.15, 'Autosomal', cex=1.5, adj=c(0.5, 0.5), xpd=NA, log='y')
       proportionalLabel(-0.0175, 0.16, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(-0.0175, 0.14, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(0.04, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.5, 1.06, expression(italic(N)==10^4) , cex=1.25, adj=c(0.5, 0.5), xpd=NA, log='y')
        # Plot labels etc.
        proportionalLabel(-0.4, 0.5, expression(Fixation~Probability), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90, log='y')
#        proportionalLabel( 0.5,  -0.225,  expression(paste("Inversion size (", italic(x),")")), log='y', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        # Legend

## Panel B
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c((1/(2*10^7)), (4/(2*10^4))), log='y', ylab='', xlab='', cex.lab=1.2)
#     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80', log='y')
        box()
        # Points
        abline(h=(1/(2*10^5)), lwd=1, lty=2, col=1)
        abline(h=1, lwd=1, lty=1, col=1)
        points(PrFixAdj ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=AdatN100kUf1)
        points(PrFixAdj ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=AdatN100kUf2)
        points(PrFixAdj ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=AdatN100kUf3)
        # axes
        axis(1, las=1, labels=NA)
#        axis(2, las=1)
         axis(2, las=1, at=c((1/(2*10^7)), (1/(2*10^6)), (1/(2*10^5)), (1/(2*10^4)), (1/(2*10^3))), 
                labels=NA)
        proportionalLabel(-0.0175, 0.16, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(-0.0175, 0.14, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(0.04, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.5, 1.06, expression(italic(N)==10^5) , cex=1.25, adj=c(0.5, 0.5), xpd=NA, log='y')
        # Plot labels etc.
        # Legend
        legend(
               x       =  1.02,
               y       =  0.0003,
               legend  =  c(
                            expression(paste(italic(U/hs), " = ", 8)),
                            expression(paste(italic(U/hs), " = ", 20)),
                            expression(paste(italic(U/hs), " = ", 40))),
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
#    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c((2/10^8), 2/(10^4)), log='y', ylab='', xlab='', cex.lab=1.2)
#     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c(0, 20), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80', log='y')
#        plotGrid(lineCol='grey80')
        box()
        # Benchmark 2/N
        abline(h=(2/100000), lwd=1, lty=2, col=1)
        # PrFix Points
        points(PrFix ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN100kUs1)
        points(PrFix ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN100kUs2)
        points(PrFix ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN100kUs3)
        # r0 PrFix Points
        points(Freq ~ rInvSize, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rdatN100kUs1)
        points(Freq ~ rInvSize, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rdatN100kUs2)
        points(Freq ~ rInvSize, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rdatN100kUs3)
        # f0*qStar
qStars  <-  qstar(YI.init=(2/datN100kUs1$N[1]), U=datN100kUs1$U[1], x=datN100kUs1$x, s=datN100kUs1$s[1], h=datN100kUs1$h[1])        
p0  <-  exp(-(datN100kUs1$U[1] * datN100kUs1$x)/(datN100kUs1$s[1]*datN100kUs1$h[1]))
points((p0*qStars) ~ datN100kUs1$x, pch=21)
qStars  <-  qstar(YI.init=(2/datN100kUs2$N[1]), U=datN100kUs2$U[1], x=datN100kUs2$x, s=datN100kUs2$s[1], h=datN100kUs2$h[1])        
p0  <-  exp(-(datN100kUs2$U[1] * datN100kUs2$x)/(datN100kUs2$s[1]*datN100kUs2$h[1]))
points((p0*qStars) ~ datN100kUs2$x, pch=22)
qStars  <-  qstar(YI.init=(2/datN100kUs3$N[1]), U=datN100kUs3$U[1], x=datN100kUs3$x, s=datN100kUs3$s[1], h=datN100kUs3$h[1])        
p0  <-  exp(-(datN100kUs3$U[1] * datN100kUs3$x)/(datN100kUs3$s[1]*datN100kUs3$h[1]))
points((p0*qStars) ~ datN100kUs3$x, pch=24)
        # axes
        axis(1, las=1)
        axis(2, las=1, at=c(2/10^8, 2/10^7, 2/10^6, 2/10^5, 2/10^4), 
                labels=c(0, expression(2/10^7), expression(2/10^6), expression(2/10^5), expression(2/10^4)))
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
#head(datN1milUs2)
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c((2/10^8), 2/(10^4)), log='y', ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80', log='y')
        box()
        # Benchmark 2/N
        abline(h=(2/1000000), lwd=1, lty=2, col=1)
        # PrFix Points
        points(PrFix ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN1milUs1)
        points(PrFix ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN1milUs2)
        points(PrFix ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN1milUs3)
        # r0 PrFix Points
        points(Freq ~ rInvSize, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rdatN1milUs1)
        points(Freq ~ rInvSize, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rdatN1milUs2)
        points(Freq ~ rInvSize, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor('tomato', opacity=0.6), data=rdatN1milUs3)
        # f0*qStar
qStars  <-  qstar(YI.init=(2/datN1milUs1$N[1]), U=datN1milUs1$U[1], x=datN1milUs1$x, s=datN1milUs1$s[1], h=datN1milUs1$h[1])        
p0  <-  exp(-(datN1milUs1$U[1] * datN1milUs1$x)/(datN1milUs1$s[1]*datN1milUs1$h[1]))
points((p0*qStars) ~ datN1milUs1$x, pch=21)
qStars  <-  qstar(YI.init=(2/datN1milUs2$N[1]), U=datN1milUs2$U[1], x=datN1milUs2$x, s=datN1milUs2$s[1], h=datN1milUs2$h[1])        
p0  <-  exp(-(datN1milUs2$U[1] * datN1milUs2$x)/(datN1milUs2$s[1]*datN1milUs2$h[1]))
points((p0*qStars) ~ datN1milUs2$x, pch=22)
qStars  <-  qstar(YI.init=(2/datN1milUs3$N[1]), U=datN1milUs3$U[1], x=datN1milUs3$x, s=datN1milUs3$s[1], h=datN1milUs3$h[1])        
p0  <-  exp(-(datN1milUs3$U[1] * datN1milUs3$x)/(datN1milUs3$s[1]*datN1milUs3$h[1]))
points((p0*qStars) ~ datN1milUs3$x, pch=24)
        # axes
        axis(1, las=1)
#        axis(2, las=1, at=c(2/10^6, 2/10^5, 2/10^4, 2/10^3, 2/10^2, 2/10^1), labels=NA)
        axis(2, las=1, at=c(2/10^7, 2/10^6, 2/10^5, 2/10^4, 2/10^3, 2/10^2, 2/10^1), labels=NA)
        # Plot labels etc.
        proportionalLabel(0.5, 1.06, expression(italic(N)==10^6) , cex=1.25, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.04, 1.05, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(-0.0175, 0.16, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(-0.0175, 0.14, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel( 0.5,  -0.225,  expression(paste("Inversion size (", italic(x),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        # Legend
        legend(
               x       =  1.02,
               y       =  0.0003,
               legend  =  c(
                            expression(paste("All ", italic(r))),
                            expression(paste(italic(r), " = ", 0)),
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
               y       =  0.00031,
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
               y       =  0.00031,
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





#' Fixation probability of different size inversions 
#' expanding the SDR on a Y chromosome: Pr(fix | x) 
#' Illustration of Pr(fix | x) for different 
#' population size and U/s ratios.
PrFixFigN500k  <-  function() {

    dat  <-  read.csv(file = './data/PrFixFig_h0.1_s0.01_N500k_deterministic_q.csv', header=TRUE)
    dat$PrFixNe  <-  dat$PrFix*(dat$N/2)
    uniqueU  <-  unique(dat$Ufac)
    dat$PrFix[dat$PrFix < (2/(10^6))]  <-  (2/(10^7))
    datN500k   <-  dat[dat$N == 500000,]
        datN500kUf1  <-  datN500k[datN500k$Ufac == uniqueU[1],]
        datN500kUf2  <-  datN500k[datN500k$Ufac == uniqueU[2],]
        datN500kUf3  <-  datN500k[datN500k$Ufac == uniqueU[3],]

    # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  colfunc(3)

    # set plot layout
#    layout.mat <- matrix(c(1:4), nrow=2, ncol=2, byrow=TRUE)
#    layout     <- layout(layout.mat,respect=TRUE)

## Panel A
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(5,4,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c((2/10^7), 2/(10^5)), log='y', ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80', log='y')
        box()
        # Points
        abline(h=(2/500000), lwd=1, lty=2, col=1)
        points(PrFix ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN500kUf1)
        points(PrFix ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN500kUf2)
        points(PrFix ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN500kUf3)
#        points(PrFix ~ x, pch=21, col=transparentColor('tomato', opacity=1), bg=transparentColor('tomato', opacity=0.6), data=IBMDat)
        # axes
        axis(1, las=1)
#        axis(2, las=1, at=c(2/10^6, 2/10^5, 2/10^4, 2/10^3, 2/10^2, 2/10^1), labels=NA)
        axis(2, las=1, at=c(2/10^7, 2/10^6, 2/10^5, 2/10^4, 2/10^3, 2/10^2, 2/10^1), labels=NA)
        # Plot labels etc.
        proportionalLabel(0.5, 1.06, expression(italic(N)==5%*%10^5) , cex=1.25, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.04, 1.05, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(-0.0175, 0.16, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(-0.0175, 0.14, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel( 0.5,  -0.225,  expression(paste("Inversion size (", italic(x),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')

        # Points
        abline(h=(2/500000), lwd=1, lty=2, col=1)
        points(PrFix ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN500kUf1)
        points(PrFix ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN500kUf2)
        points(PrFix ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN500kUf3)
#        points(PrFix ~ x, pch=21, col=transparentColor('tomato', opacity=1), bg=transparentColor('tomato', opacity=0.6), data=IBMDat)
        # axes
        axis(1, las=1)
#        axis(2, las=1, at=c(2/10^6, 2/10^5, 2/10^4, 2/10^3, 2/10^2, 2/10^1), labels=NA)
        axis(2, las=1, at=c(2/10^7, 2/10^6, 2/10^5, 2/10^4, 2/10^3, 2/10^2, 2/10^1), labels=NA)
        # Plot labels etc.
        proportionalLabel(0.5, 1.06, expression(italic(N)==50^5) , cex=1.25, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.04, 1.05, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(-0.0175, 0.16, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel(-0.0175, 0.14, '_', cex=1.2, adj=c(0.5, 0.5), xpd=NA,log='y', srt=30)
        proportionalLabel( 0.5,  -0.225,  expression(paste("Inversion size (", italic(x),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')

}











#' Fixation probability of different size inversions 
#' expanding the SDR on a Y chromosome: Pr(fix | x) 
#' Illustration of Pr(fix | x) for different 
#' population size and U/s ratios.
PrFixNeFig  <-  function() {

    # Import data for plotting
    dat  <-  read.csv(file = './data/PrFixFig_h0.1_s0.01_N1k_deterministic_q.csv', header=TRUE)
    dat$PrFixNe  <-  dat$PrFix*(dat$N/2)
    uniqueU  <-  unique(dat$Ufac)
    datN1k   <-  dat[dat$N == 1000,]
        datN1kUf1  <-  datN1k[datN1k$Ufac == uniqueU[1],]
        datN1kUf2  <-  datN1k[datN1k$Ufac == uniqueU[2],]
        datN1kUf3  <-  datN1k[datN1k$Ufac == uniqueU[3],]

    dat  <-  read.csv(file = './data/PrFixFig_h0.1_s0.01_N10k_deterministic_q.csv', header=TRUE)
    dat$PrFixNe  <-  dat$PrFix*(dat$N/2)
    uniqueU  <-  unique(dat$Ufac)
    datN10k   <-  dat[dat$N == 10000,]
        datN10kUf1  <-  datN10k[datN10k$Ufac == uniqueU[1],]
        datN10kUf2  <-  datN10k[datN10k$Ufac == uniqueU[2],]
        datN10kUf3  <-  datN10k[datN10k$Ufac == uniqueU[3],]

    IBMDat  <-  read.csv(file = './data/PrFixIBM_h0.1_s0.01_N1000_Ufac2__xTest.csv', header=TRUE)
    IBMDat$PrFixNe  <-  IBMDat$PrFix*(IBMDat$N/2)


    Adat  <-  read.csv(file = './data/PrFixAuto3Fig_h0.1_s0.01.csv', header=TRUE)
    Adat$PrFixNe  <-  Adat$PrFix*(2*Adat$N)
    Adat$PrFixAdj  <-  Adat$PrFix
    Adat$PrFixAdj[Adat$PrFixAdj < (1/(10^6))]  <-  (1/(10^6))
    uniqueU    <-  unique(Adat$Ufac)
    AdatN1k    <-  Adat[Adat$N == 1000,]
        AdatN1kUf1  <-  AdatN1k[AdatN1k$Ufac == uniqueU[1],]
        AdatN1kUf2  <-  AdatN1k[AdatN1k$Ufac == uniqueU[2],]
        AdatN1kUf3  <-  AdatN1k[AdatN1k$Ufac == uniqueU[3],]
    AdatN10k   <-  Adat[Adat$N == 10000,]
        AdatN10kUf1  <-  AdatN10k[AdatN10k$Ufac == uniqueU[1],]
        AdatN10kUf2  <-  AdatN10k[AdatN10k$Ufac == uniqueU[2],]
        AdatN10kUf3  <-  AdatN10k[AdatN10k$Ufac == uniqueU[3],]


    # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  colfunc(3)

    # set plot layout
    layout.mat <- matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

## Panel A
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(4,4,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        abline(h=1, lwd=1, lty=1, col=1)
        points(PrFixNe ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[3], opacity=0.6), data=AdatN1kUf1)
        points(PrFixNe ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[3], opacity=0.6), data=AdatN1kUf2)
        points(PrFixNe ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[3], opacity=0.6), data=AdatN1kUf3)
        points(PrFixNe ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.6), data=AdatN10kUf1)
        points(PrFixNe ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.6), data=AdatN10kUf2)
        points(PrFixNe ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.6), data=AdatN10kUf3)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        proportionalLabel(0.04, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.1, 'Autosomal', cex=1.25, adj=c(0.5, 0.5), xpd=NA)
        # Plot labels etc.
        proportionalLabel(-0.45, 0.5, expression(Fixation~Probability/italic(2*N[e])), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel( 0.5,  -0.225,  expression(paste("Inversion size (", italic(x),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(paste(italic(N), " = ", 10^3)),
                            expression(paste(italic(N), " = ", 10^4))),
               pch     =  22,
               col     =  transparentColor(COLS[1], opacity=0),
               pt.bg   =  c(transparentColor(COLS[3], opacity=0.6),
                            transparentColor(COLS[1], opacity=0.6)),
               pt.cex  = 1.5,
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
        legend(
               x       =  usr[2],
               y       =  usr[4]*0.8,
               legend  =  c(
                            expression(paste(italic(U/hs), " = ", 20)),
                            expression(paste(italic(U/hs), " = ", 50)),
                            expression(paste(italic(U/hs), " = ", 100))),
               pch     =  c(21,22,24),
               col     =  transparentColor(COLS[1], opacity=1),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


## Panel B
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
#    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c(0, 10), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        abline(h=1, lwd=1, lty=1, col=1)
        points(PrFixNe ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[3], opacity=0.6), data=datN1kUf1)
        points(PrFixNe ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[3], opacity=0.6), data=datN1kUf2)
        points(PrFixNe ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[3], opacity=0.6), data=datN1kUf3)
        points(PrFixNe ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.6), data=datN10kUf1)
        points(PrFixNe ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.6), data=datN10kUf2)
        points(PrFixNe ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.6), data=datN10kUf3)
        points(PrFixNe ~ x, pch=21, col=transparentColor('tomato', opacity=1), bg=transparentColor('tomato', opacity=0.6), data=IBMDat)

        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.5, 1.1, 'Expand SLR', cex=1.25, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.25,  0.5,   expression(2%*%Fixation~Probability/italic(N[e])), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.225,  expression(paste("Inversion size (", italic(x),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)

}


#' Fixation probability of different size inversions 
#' expanding the SDR on a Y chromosome: Pr(fix | x) 
#' Illustration of Pr(fix | x) for different 
#' population size and U/s ratios.
PrFixAutoFig  <-  function() {

	# Import data for plotting
	dat  <-  read.csv(file = './data/PrFixAuto3Fig_h0.1_s0.01.csv', header=TRUE)
	dat$PrFixNe  <-  dat$PrFix*(2*dat$N)
	uniqueU  <-  unique(dat$Ufac)
	datN100   <-  dat[dat$N == 100,]
		datN100Uf1  <-  datN100[datN100$Ufac == uniqueU[1],]
		datN100Uf2  <-  datN100[datN100$Ufac == uniqueU[2],]
		datN100Uf3  <-  datN100[datN100$Ufac == uniqueU[3],]
	datN1k    <-  dat[dat$N == 1000,]
		datN1kUf1  <-  datN1k[datN1k$Ufac == uniqueU[1],]
		datN1kUf2  <-  datN1k[datN1k$Ufac == uniqueU[2],]
		datN1kUf3  <-  datN1k[datN1k$Ufac == uniqueU[3],]
	datN10k   <-  dat[dat$N == 10000,]
		datN10kUf1  <-  datN10k[datN10k$Ufac == uniqueU[1],]
		datN10kUf2  <-  datN10k[datN10k$Ufac == uniqueU[2],]
		datN10kUf3  <-  datN10k[datN10k$Ufac == uniqueU[3],]

	# Colors
	colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
	COLS     <-  colfunc(3)

    # set plot layout
    layout.mat <- matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

## Panel A
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c(0,5), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
		abline(h=1, lwd=1.5, lty=1, col=1)
        points(PrFixNe ~ x, pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.6), data=datN100Uf1)
        points(PrFixNe ~ x, pch=22, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.6), data=datN100Uf2)
        points(PrFixNe ~ x, pch=23, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.6), data=datN100Uf3)
        points(PrFixNe ~ x, pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN1kUf1)
        points(PrFixNe ~ x, pch=22, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN1kUf2)
        points(PrFixNe ~ x, pch=23, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN1kUf3)
        points(PrFixNe ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.6), data=datN10kUf1)
        points(PrFixNe ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.6), data=datN10kUf2)
        points(PrFixNe ~ x, pch=23, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.6), data=datN10kUf3)
        # axes
        axis(1, las=1)
        axis(2, las=1)
#        axis(2, las=1, at=c((1/2)*10^-4, (1/2)*10^-3, (1/2)*10^-2, (1/2)*10^-1), 
#        		labels=c(expression(1/2%*%10^-4), expression(1/2%*%10^-3), expression(1/2%*%10^-2), expression(2%*%10^-1)))
        # Plot labels etc.
        proportionalLabel(-0.25,  0.5,   expression(Fixation~Probability%*%2*italic(N)), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.25,  expression(paste("Inversion size (", italic(x),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
               x       =  usr[2]*0.965,
               y       =  usr[4],
               legend  =  c(
                            expression(paste(italic(N), " = ", 10^4)),
                            expression(paste(italic(N), " = ", 10^3)),
                            expression(paste(italic(N), " = ", 10^2))),
               pch     =  22,
               col     =  transparentColor(COLS[1], opacity=0),
               pt.bg   =  c(transparentColor(COLS[1], opacity=0.6),
                            transparentColor(COLS[2], opacity=0.6),
                            transparentColor(COLS[3], opacity=0.6)),
               pt.cex  = 1.5,
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
        legend(
               x       =  usr[2],
               y       =  usr[4]*0.775
               ,
               legend  =  c(
                            expression(paste(italic(U/hs), " = ", 20)),
                            expression(paste(italic(U/hs), " = ", 50)),
                            expression(paste(italic(U/hs), " = ", 100))),
               pch     =  c(21,22,23),
               col     =  transparentColor(COLS[1], opacity=1),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        abline(h=1, lwd=1.5, lty=1, col=1)
        points(PrFixNe ~ x, pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.6), data=datN100Uf1)
        points(PrFixNe ~ x, pch=22, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.6), data=datN100Uf2)
        points(PrFixNe ~ x, pch=23, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.6), data=datN100Uf3)
        points(PrFixNe ~ x, pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN1kUf1)
        points(PrFixNe ~ x, pch=22, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN1kUf2)
        points(PrFixNe ~ x, pch=23, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=datN1kUf3)
        points(PrFixNe ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.6), data=datN10kUf1)
        points(PrFixNe ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.6), data=datN10kUf2)
        points(PrFixNe ~ x, pch=23, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.6), data=datN10kUf3)
        # axes
        axis(1, las=1)
        axis(2, las=1)
#        axis(2, las=1, at=c((1/2)*10^-4, (1/2)*10^-3, (1/2)*10^-2, (1/2)*10^-1), 
#               labels=c(expression(1/2%*%10^-4), expression(1/2%*%10^-3), expression(1/2%*%10^-2), expression(2%*%10^-1)))
        # Plot labels etc.
        proportionalLabel( 0.5,  -0.25,  expression(paste("Inversion size (", italic(x),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
#        legend(
#               x       =  usr[2]*0.99,
#               y       =  usr[4],
#               legend  =  c(
#                            expression(paste(italic(N[e]), " = ", 10^2)),
#                            expression(paste(italic(N[e]), " = ", 10^3)),
#                            expression(paste(italic(N[e]), " = ", 10^4))),
#               pch     =  22,
#               col     =  transparentColor(COLS[1], opacity=0),
#               pt.bg   =  c(transparentColor(COLS[1], opacity=0.6),
#                            transparentColor(COLS[2], opacity=0.6),
#                            transparentColor(COLS[3], opacity=0.6)),
#               pt.cex  = 1.5,
#               cex     =  1,
#               xjust   =  1,
#               yjust   =  1,
#               bty     =  'n',
#               border  =  NA
#               )
#        legend(
#               x       =  usr[2],
#               y       =  usr[4]*0.85,
#               legend  =  c(
#                            expression(paste(italic(U/s), " = ", 10)),
#                            expression(paste(italic(U/s), " = ", 5)),
#                            expression(paste(italic(U/s), " = ", 2))),
#               pch     =  c(23,22,21),
#               col     =  transparentColor(COLS[1], opacity=1),
#               cex     =  1,
#               xjust   =  1,
#               yjust   =  1,
#               bty     =  'n',
#               border  =  NA
#               )

}




#########################################
## SUPPLEMENTARY FIGURES

#' Deterministic dynamics of inversions expanding the SDR on a Y chromosome:
#' Illustration of: (A) inversion relative fitness
#'                  (B) Inversion frequency dynamics, and 
#'                  (C) del. allele frequency dynamics.
deterministicSuppFig  <-  function() {


    # set plot layout
    layout.mat <- matrix(c(1:20), nrow=4, ncol=5, byrow=FALSE)
    layout     <- layout(layout.mat,respect=TRUE)

    # Make data for plotting (x = 0.1; 1st column of panels)
    h        <-  0.1
    s        <-  0.01
    Ufactor  <-  2
    nTot     <-  10^4 
    U        <-  Ufactor*s
    u        <-  U/nTot
    qHat     <-  (U/(nTot*h*s))
    x        <-  0.05
    Prs  <-  round(
                c(dpois(x=0, lambda=U*x/(s*h)),
                  dpois(x=1, lambda=U*x/(s*h))
                  ), digits=2)

    # make data for plotting
    r.0.Dat  <-  makeDeterministicFigSimData(r = 0, x = x, generations = 8000)
    r.1.Dat  <-  makeDeterministicFigSimData(r = 1, x = x, generations = 8000)

    r.0.wbar.YI.t   <-  r.0.Dat$wbar.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.1.wbar.YI.t   <-  r.1.Dat$wbar.YI.t[round(r.1.Dat$YI.t,digits=4) > 0]

    # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  rev(colfunc(5))

## Panel A
# Distribution of captured deleterious alleles (r) given inversion size (x)
# with benchmarks for 5th, 10th, 25th, and 40th quantiles, as well as the mean.
    pDens  <-  dpois(x=c(0:10), lambda=U*x/(s*h))
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,10), ylim = c(0, 0.4), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(pDens~c(0:10), lwd=1.5, lty=1, col=transparentColor("#252525", opacity=1))
        points(pDens~c(0:10), lwd=1.5, lty=1, pch=21, col=transparentColor("#252525", opacity=1), bg=transparentColor("#252525", opacity=0.25))
        abline(v=(U*x/(s*h)), col=1, lwd=1)
        abline(v=c(0,1), col=2, lwd=1, lty=2)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.5, 1.1, substitute(italic(x)==xx, list(xx=x)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.04, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(paste(Pr(italic(r)~"|"~italic(x)))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel( 0.5,  -0.2,  expression(italic(r)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(r)=="{0,1}"),
                            expression(bar(italic(r))==frac(italic(U*x),italic(h*s)))),
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
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(min(na.omit(r.1.wbar.YI.t)), max(na.omit(r.0.wbar.YI.t))), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(r.0.wbar.YI.t ~ seq_along(r.0.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1))
        lines(r.1.wbar.YI.t ~ seq_along(r.1.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1))
        abline(h=1, lwd=1, lty=2, col=2)
        # Points depicting (arbitrary) extinction frequency of 0.00001
#        points(r.0.wbar.YI.t[length(r.0.wbar.YI.t)] ~ length(r.0.wbar.YI.t), pch=21, col=transparentColor(COLS[5], opacity=1), bg=transparentColor(COLS[5], opacity=0.7))
        points(r.1.wbar.YI.t[length(r.1.wbar.YI.t)] ~ length(r.1.wbar.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
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
                            expression(italic(r)==0),
                            expression(italic(r)==1)),
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
    # make plot
#    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        lines(YI.t[1:length(r.0.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.0.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.1.wbar.YI.t)] ~ seq_along(YI.t[1:length(r.1.wbar.YI.t)]), lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.1.Dat)
        # Points
#        points(YI.t[length(r.0.wbar.YI.t)] ~ length(r.0.wbar.YI.t), pch=21, col=transparentColor(COLS[5], opacity=1), bg=transparentColor(COLS[5], opacity=0.7), data=r.0.Dat)
        points(YI.t[length(r.1.wbar.YI.t)] ~ seq_along(YI.t)[length(r.1.wbar.YI.t)], pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.1.Dat)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(paste("Inversion Frequency (", italic(Y[I]), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        


    # reset Colors
     COLS  <-  c('red',
                 '#252525')
## Panel D
# deleterious allele frequencies for r = 5
     # Make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,0.002), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(XaOv.wt.t ~ seq_along(XaOv.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.0.Dat)
        lines(XaOv.del.t ~ seq_along(XaOv.del.t), lwd=1.5, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.0.Dat)
        lines(qI.wt.t ~ seq_along(qI.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(italic(r)==0), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(paste("Deleterious Allele Frequency")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4]*0.9,
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
    Prs  <-  round(
                c(dpois(x=0, lambda=U*x/(s*h)),
                  dpois(x=1, lambda=U*x/(s*h)),
                  dpois(x=2, lambda=U*x/(s*h))
                  ), digits=2)

    # make data for plotting
    r.0.Dat  <-  makeDeterministicFigSimData(r = 0, x = x, generations = 8000)
    r.1.Dat  <-  makeDeterministicFigSimData(r = 1, x = x, generations = 8000)
    r.2.Dat  <-  makeDeterministicFigSimData(r = 2, x = x, generations = 8000)

    r.0.wbar.YI.t  <-  r.0.Dat$wbar.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.1.wbar.YI.t  <-  r.1.Dat$wbar.YI.t[round(r.1.Dat$YI.t,digits=4) > 0]
    r.2.wbar.YI.t  <-  r.2.Dat$wbar.YI.t[round(r.2.Dat$YI.t,digits=4) > 0]

   # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  rev(colfunc(5))

## Panel E
# Distribution of captured deleterious alleles (r) given inversion size (x)
# with benchmarks for 5th, 10th, 25th, and 40th quantiles, as well as the mean.
   pDens  <-  dpois(x=c(0:10), lambda=U*x/(s*h))
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,10), ylim = c(0, 0.4), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(pDens~c(0:10), lwd=1.5, lty=1, col=transparentColor("#252525", opacity=1))
        points(pDens~c(0:10), lwd=1.5, lty=1, pch=21, col=transparentColor("#252525", opacity=1), bg=transparentColor("#252525", opacity=0.25))
        abline(v=(U*x/(s*h)), col=1, lwd=1)
        abline(v=c(0:2), col=2, lwd=1, lty=2)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.5, 1.1, substitute(italic(x)==xx, list(xx=x)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.04, 1.05, 'E', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.2,  expression(italic(r)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(r)=="{0,1,2}"),
                            expression(bar(italic(r))==frac(italic(U*x),italic(h*s)))),
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
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(min(na.omit(r.2.wbar.YI.t)), max(na.omit(r.0.wbar.YI.t))), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(r.0.wbar.YI.t ~ seq_along(r.0.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1))
        lines(r.1.wbar.YI.t ~ seq_along(r.1.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1))
        lines(r.2.wbar.YI.t ~ seq_along(r.2.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1))
        abline(h=1, lwd=1, lty=2, col=2)
        # Points depicting (arbitrary) extinction frequency of 0.00001
#        points(r.0.wbar.YI.t[length(r.0.wbar.YI.t)] ~ length(r.0.wbar.YI.t), pch=8, col='tomato')
        points(r.1.wbar.YI.t[length(r.1.wbar.YI.t)] ~ length(r.1.wbar.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        points(r.2.wbar.YI.t[length(r.2.wbar.YI.t)] ~ length(r.2.wbar.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
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
                            expression(italic(r)==0),
                            expression(italic(r)==1),
                            expression(italic(r)==2)),
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
    # make plot
#    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(YI.t[1:length(r.0.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.0.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.1.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.1.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.1.Dat)
        lines(YI.t[1:length(r.2.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.2.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.2.Dat)
        # Points
#        points(YI.t[length(r.0.wbar.YI.t)] ~ length(r.0.wbar.YI.t), pch=8, col='tomato', data=r.0.Dat)
        points(YI.t[length(r.1.wbar.YI.t)] ~ seq_along(YI.t)[length(r.1.wbar.YI.t)], pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.1.Dat)
        points(YI.t[length(r.2.wbar.YI.t)] ~ seq_along(YI.t)[length(r.2.wbar.YI.t)], pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.2.Dat)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'G', cex=1.2, adj=c(0.5, 0.5), xpd=NA)


    # reset Colors
     COLS  <-  c('red',
                 '#252525')
## Panel H
# deleterious allele frequencies for r = 5
     # Make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.1.Dat)), ylim = c(0,0.002), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(XaOv.wt.t ~ seq_along(XaOv.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(XaOv.del.t ~ seq_along(XaOv.del.t), lwd=1.5, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qI.wt.t ~ seq_along(qI.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'H', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1, expression(italic(r)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        




## 3rd Column of Panels: x = 0.2
    # Set inversion size, r values
    x        <-  0.2
    Prs  <-  round(
                c(dpois(x=0, lambda=U*x/(s*h)),
                  dpois(x=1, lambda=U*x/(s*h)),
                  dpois(x=2, lambda=U*x/(s*h)),
                  dpois(x=3, lambda=U*x/(s*h))
                  ), digits=2)

    # make data for plotting
    r.0.Dat  <-  makeDeterministicFigSimData(r = 0, x = x, generations = 8000)
    r.1.Dat  <-  makeDeterministicFigSimData(r = 1, x = x, generations = 8000)
    r.2.Dat  <-  makeDeterministicFigSimData(r = 2, x = x, generations = 8000)
    r.3.Dat  <-  makeDeterministicFigSimData(r = 3, x = x, generations = 8000)

    r.0.wbar.YI.t  <-  r.0.Dat$wbar.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.1.wbar.YI.t  <-  r.1.Dat$wbar.YI.t[round(r.1.Dat$YI.t,digits=4) > 0]
    r.2.wbar.YI.t  <-  r.2.Dat$wbar.YI.t[round(r.2.Dat$YI.t,digits=4) > 0]
    r.3.wbar.YI.t  <-  r.3.Dat$wbar.YI.t[round(r.3.Dat$YI.t,digits=4) > 0]

   # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  rev(colfunc(5))

## Panel I
# Distribution of captured deleterious alleles (r) given inversion size (x)
# with benchmarks for 5th, 10th, 25th, and 40th quantiles, as well as the mean.
   pDens  <-  dpois(x=c(0:10), lambda=U*x/(s*h))
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,10), ylim = c(0, 0.4), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(pDens~c(0:10), lwd=1.5, lty=1, col=transparentColor("#252525", opacity=1))
        points(pDens~c(0:10), lwd=1.5, lty=1, pch=21, col=transparentColor("#252525", opacity=1), bg=transparentColor("#252525", opacity=0.25))
        abline(v=(U*x/(s*h)), col=1, lwd=1)
        abline(v=c(0:3), col=2, lwd=1, lty=2)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.5, 1.1, substitute(italic(x)==xx, list(xx=x)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.04, 1.05, 'I', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.2,  expression(italic(r)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(r)=="{0,1,2,3}"),
                            expression(bar(italic(r))==frac(italic(U*x),italic(h*s)))),
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
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(min(na.omit(r.3.wbar.YI.t)), max(na.omit(r.0.wbar.YI.t))), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(r.0.wbar.YI.t ~ seq_along(r.0.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1))
        lines(r.1.wbar.YI.t ~ seq_along(r.1.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1))
        lines(r.2.wbar.YI.t ~ seq_along(r.2.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.3.wbar.YI.t ~ seq_along(r.3.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1))
        abline(h=1, lwd=1, lty=2, col=2)
        # Points depicting (arbitrary) extinction frequency of 0.00001
#        points(r.0.wbar.YI.t[length(r.0.wbar.YI.t)] ~ length(r.0.wbar.YI.t), pch=8, col='tomato')
        points(r.1.wbar.YI.t[length(r.1.wbar.YI.t)] ~ length(r.1.wbar.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        points(r.2.wbar.YI.t[length(r.2.wbar.YI.t)] ~ length(r.2.wbar.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
        points(r.3.wbar.YI.t[length(r.3.wbar.YI.t)] ~ length(r.3.wbar.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
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
                            expression(italic(r)==0),
                            expression(italic(r)==1),
                            expression(italic(r)==2),
                            expression(italic(r)==3)),
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
    # make plot
#    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(YI.t[1:length(r.0.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.0.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.1.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.1.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.1.Dat)
        lines(YI.t[1:length(r.2.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.2.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.2.Dat)
        lines(YI.t[1:length(r.3.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.3.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.3.Dat)
        # Points
#        points(YI.t[length(r.0.wbar.YI.t)] ~ length(r.0.wbar.YI.t), pch=8, col='tomato', data=r.0.Dat)
        points(YI.t[length(r.1.wbar.YI.t)] ~ seq_along(YI.t)[length(r.1.wbar.YI.t)], pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.1.Dat)
        points(YI.t[length(r.2.wbar.YI.t)] ~ seq_along(YI.t)[length(r.2.wbar.YI.t)], pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.2.Dat)
        points(YI.t[length(r.3.wbar.YI.t)] ~ seq_along(YI.t)[length(r.3.wbar.YI.t)], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.3.Dat)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'K', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # reset Colors
     COLS  <-  c('red',
                 '#252525')
## Panel L
# deleterious allele frequencies for r = 5
     # Make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,0.002), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(XaOv.wt.t ~ seq_along(XaOv.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(XaOv.del.t ~ seq_along(XaOv.del.t), lwd=1.5, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qI.wt.t ~ seq_along(qI.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'L', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(italic(r)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        



## 4th Column of Panels: x = 0.1
    # Set inversion size, r values
    x        <-  0.5
    Prs  <-  round(
                c(dpois(x=0, lambda=U*x/(s*h)),
                  dpois(x=2, lambda=U*x/(s*h)),
                  dpois(x=4, lambda=U*x/(s*h)),
                  dpois(x=6, lambda=U*x/(s*h)),
                  dpois(x=8, lambda=U*x/(s*h))
                  ), digits=2)

    # make data for plotting
    r.0.Dat  <-  makeDeterministicFigSimData(r = 0, x = x, generations = 8000)
    r.2.Dat  <-  makeDeterministicFigSimData(r = 2, x = x, generations = 8000)
    r.4.Dat  <-  makeDeterministicFigSimData(r = 4, x = x, generations = 8000)
    r.6.Dat  <-  makeDeterministicFigSimData(r = 6, x = x, generations = 8000)
    r.8.Dat  <-  makeDeterministicFigSimData(r = 8, x = x, generations = 8000)

    r.0.wbar.YI.t   <-  r.0.Dat$wbar.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.2.wbar.YI.t   <-  r.2.Dat$wbar.YI.t[round(r.2.Dat$YI.t,digits=4) > 0]
    r.4.wbar.YI.t   <-  r.4.Dat$wbar.YI.t[round(r.4.Dat$YI.t,digits=4) > 0]
    r.6.wbar.YI.t   <-  r.6.Dat$wbar.YI.t[round(r.6.Dat$YI.t,digits=4) > 0]   
    r.8.wbar.YI.t   <-  r.8.Dat$wbar.YI.t[round(r.8.Dat$YI.t,digits=4) > 0]   

   # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  rev(colfunc(5))

## Panel M
# Distribution of captured deleterious alleles (r) given inversion size (x)
# with benchmarks for 5th, 10th, 25th, and 40th quantiles, as well as the mean.
   pDens  <-  dpois(x=c(0:25), lambda=U*x/(s*h))
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,35), ylim = c(0, 0.4), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(pDens~c(0:25), lwd=1.5, lty=1, col=transparentColor("#252525", opacity=1))
        points(pDens~c(0:25), lwd=1.5, lty=1, pch=21, col=transparentColor("#252525", opacity=1), bg=transparentColor("#252525", opacity=0.25))
        abline(v=(U*x/(s*h)), col=1, lwd=1)
        abline(v=c(0,2,4,6,8), col=2, lwd=1, lty=2)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.5, 1.1, substitute(italic(x)==xx, list(xx=x)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.04, 1.05, 'M', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.2,  expression(italic(r)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(r)=="{0,2,4,6,8}"),
                            expression(bar(italic(r))==frac(italic(U*x),italic(h*s)))),
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
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(min(na.omit(r.8.wbar.YI.t)), max(na.omit(r.0.wbar.YI.t))), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(r.0.wbar.YI.t ~ seq_along(r.0.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1))
        lines(r.2.wbar.YI.t ~ seq_along(r.2.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1))
        lines(r.4.wbar.YI.t ~ seq_along(r.4.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.6.wbar.YI.t ~ seq_along(r.6.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r.8.wbar.YI.t ~ seq_along(r.8.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1))
        abline(h=1, lwd=1, lty=2, col=2)
        # Points depicting (arbitrary) extinction frequency of 0.00001
#        points(r.0.wbar.YI.t[length(r.0.wbar.YI.t)] ~ length(r.0.wbar.YI.t), pch=8, col='tomato')
        points(r.2.wbar.YI.t[length(r.2.wbar.YI.t)] ~ length(r.2.wbar.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        points(r.4.wbar.YI.t[length(r.4.wbar.YI.t)] ~ length(r.4.wbar.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
        points(r.6.wbar.YI.t[length(r.6.wbar.YI.t)] ~ length(r.6.wbar.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r.8.wbar.YI.t[length(r.8.wbar.YI.t)] ~ length(r.8.wbar.YI.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7))
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
                            expression(italic(r)==0),
                            expression(italic(r)==2),
                            expression(italic(r)==4),
                            expression(italic(r)==6),
                            expression(italic(r)==8)),
               lty     =  1,
               lwd     =  1.5,
               col     =  c(transparentColor(COLS[5], opacity=1),
                            transparentColor(COLS[4], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[1], opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

## Panel O
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
#    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(YI.t[1:length(r.0.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.0.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.2.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.2.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.2.Dat)
        lines(YI.t[1:length(r.4.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.4.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.4.Dat)
        lines(YI.t[1:length(r.6.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.6.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.6.Dat)
        lines(YI.t[1:length(r.8.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.8.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.8.Dat)
        # Points
#        points(YI.t[length(r.0.wbar.YI.t)] ~ length(r.0.wbar.YI.t), pch=8, col='tomato', data=r.0.Dat)
        points(YI.t[length(r.2.wbar.YI.t)] ~ seq_along(YI.t)[length(r.2.wbar.YI.t)], pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.2.Dat)
        points(YI.t[length(r.4.wbar.YI.t)] ~ seq_along(YI.t)[length(r.4.wbar.YI.t)], pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.4.Dat)
        points(YI.t[length(r.6.wbar.YI.t)] ~ seq_along(YI.t)[length(r.6.wbar.YI.t)], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.6.Dat)
        points(YI.t[length(r.8.wbar.YI.t)] ~ seq_along(YI.t)[length(r.8.wbar.YI.t)], pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.8.Dat)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'O', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        # Legend

    # reset Colors
     COLS  <-  c('red',
                 '#252525')
## Panel P
# deleterious allele frequencies for r = 5
     # Make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,0.002), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(XaOv.wt.t ~ seq_along(XaOv.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.2.Dat)
        lines(XaOv.del.t ~ seq_along(XaOv.del.t), lwd=1.5, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.2.Dat)
        lines(qI.wt.t ~ seq_along(qI.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.2.Dat)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'P', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(italic(r)==2), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        



## 5th Column of Panels: x = 0.8
    # Set inversion size, r values
    x    <-  0.8
    Prs  <-  round(
                c(dpois(x=0, lambda=U*x/(s*h)),
                  dpois(x=1, lambda=U*x/(s*h)),
                  dpois(x=4, lambda=U*x/(s*h)),
                  dpois(x=8, lambda=U*x/(s*h)),
                  dpois(x=12, lambda=U*x/(s*h))
                  ), digits=2)

    # make data for plotting
    r.0.Dat  <-  makeDeterministicFigSimData(r = 0, x = x, generations = 8000)
    r.1.Dat  <-  makeDeterministicFigSimData(r = 1, x = x, generations = 8000)
    r.4.Dat  <-  makeDeterministicFigSimData(r = 4, x = x, generations = 8000)
    r.8.Dat  <-  makeDeterministicFigSimData(r = 8, x = x, generations = 8000)
    r.12.Dat <-  makeDeterministicFigSimData(r = 12, x = x, generations = 8000)

    r.0.wbar.YI.t  <-  r.0.Dat$wbar.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.1.wbar.YI.t  <-  r.1.Dat$wbar.YI.t[round(r.1.Dat$YI.t,digits=4) > 0]
    r.4.wbar.YI.t  <-  r.4.Dat$wbar.YI.t[round(r.4.Dat$YI.t,digits=4) > 0]
    r.8.wbar.YI.t  <-  r.8.Dat$wbar.YI.t[round(r.8.Dat$YI.t,digits=4) > 0]   
    r.12.wbar.YI.t <-  r.12.Dat$wbar.YI.t[round(r.12.Dat$YI.t,digits=4) > 0]   

   # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  rev(colfunc(5))

## Panel Q
# Distribution of captured deleterious alleles (r) given inversion size (x)
# with benchmarks for 5th, 10th, 25th, and 40th quantiles, as well as the mean.
   pDens  <-  dpois(x=c(0:40), lambda=U*x/(s*h))
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,55), ylim = c(0, 0.4), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(pDens~c(0:40), lwd=1.5, lty=1, col=transparentColor("#252525", opacity=1))
        points(pDens~c(0:40), lwd=1.5, lty=1, pch=21, col=transparentColor("#252525", opacity=1), bg=transparentColor("#252525", opacity=0.25))
        abline(v=(U*x/(s*h)), col=1, lwd=1)
        abline(v=c(0,1,4,8,12), col=2, lwd=1, lty=2)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.5, 1.1, substitute(italic(x)==xx, list(xx=x)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.04, 1.05, 'Q', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.2,  expression(italic(r)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(r)=="{0,1,4,8,12}"),
                            expression(bar(italic(r))==frac(italic(U*x),italic(h*s)))),
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
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(min(na.omit(r.12.wbar.YI.t)), max(na.omit(r.0.wbar.YI.t))), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(r.0.wbar.YI.t ~ seq_along(r.0.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1))
        lines(r.1.wbar.YI.t ~ seq_along(r.1.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1))
        lines(r.4.wbar.YI.t ~ seq_along(r.4.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.8.wbar.YI.t ~ seq_along(r.8.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r.12.wbar.YI.t ~ seq_along(r.12.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1))
        abline(h=1, lwd=1, lty=2, col=2)
        # Points depicting (arbitrary) extinction frequency of 0.00001
#        points(r.0.wbar.YI.t[length(r.0.wbar.YI.t)] ~ length(r.0.wbar.YI.t), pch=8, col='tomato')
#        points(r.1.wbar.YI.t[length(r.1.wbar.YI.t)] ~ length(r.1.wbar.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        points(r.4.wbar.YI.t[length(r.4.wbar.YI.t)] ~ length(r.4.wbar.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
        points(r.8.wbar.YI.t[length(r.8.wbar.YI.t)] ~ length(r.8.wbar.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r.12.wbar.YI.t[length(r.12.wbar.YI.t)] ~ length(r.12.wbar.YI.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7))
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
                            expression(italic(r)==0),
                            expression(italic(r)==1),
                            expression(italic(r)==4),
                            expression(italic(r)==8),
                            expression(italic(r)==12)),
               lty     =  1,
               lwd     =  1.5,
               col     =  c(transparentColor(COLS[5], opacity=1),
                            transparentColor(COLS[4], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[1], opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

## Panel S
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
#    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        # Points
        lines(YI.t[1:length(r.0.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.0.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.1.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.1.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.1.Dat)
        lines(YI.t[1:length(r.4.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.4.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.4.Dat)
        lines(YI.t[1:length(r.8.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.8.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.8.Dat)
        lines(YI.t[1:length(r.12.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.12.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.12.Dat)
        # Points
#        points(YI.t[length(r.0.wbar.YI.t)] ~ length(r.0.wbar.YI.t), pch=8, col='tomato', data=r.0.Dat)
#        points(YI.t[length(r.1.wbar.YI.t)] ~ length(r.1.wbar.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        points(YI.t[length(r.4.wbar.YI.t)] ~ length(r.4.wbar.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.4.Dat)
        points(YI.t[length(r.8.wbar.YI.t)] ~ length(r.8.wbar.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.8.Dat)
        points(YI.t[length(r.12.wbar.YI.t)] ~ length(r.12.wbar.YI.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.12.Dat)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'S', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        # Legend

    # reset Colors
     COLS  <-  c('red',
                 '#252525')
## Panel T
# deleterious allele frequencies for r = 5
     # Make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.4.Dat)), ylim = c(0,0.002), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(XaOv.wt.t ~ seq_along(XaOv.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.4.Dat)
        lines(XaOv.del.t ~ seq_along(XaOv.del.t), lwd=1.5, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.4.Dat)
        lines(qI.wt.t ~ seq_along(qI.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.4.Dat)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'T', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(italic(r)==4), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        

}




deterministicSuppFig_h0.25  <-  function() {


    # set plot layout
    layout.mat <- matrix(c(1:20), nrow=4, ncol=5, byrow=FALSE)
    layout     <- layout(layout.mat,respect=TRUE)

    # Make data for plotting (x = 0.1; 1st column of panels)
    h        <-  0.25
    s        <-  0.01
    U        <-  0.02
    nTot     <-  10^4 
    u        <-  U/nTot
    qHat     <-  (U/(nTot*h*s))
    x        <-  0.05

    # make data for plotting
    r.0.Dat  <-  makeDeterministicFigSimData(r = 0, U = U, x = x, h = h, generations = 3500)
    r.1.Dat  <-  makeDeterministicFigSimData(r = 1, U = U, x = x, h = h, generations = 3500)

    r.0.wbar.YI.t   <-  r.0.Dat$wbar.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.1.wbar.YI.t   <-  r.1.Dat$wbar.YI.t[round(r.1.Dat$YI.t,digits=4) > 0]

    # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  rev(colfunc(5))


## Panel A
# Distribution of captured deleterious alleles (r) given inversion size (x)
# with benchmarks for 5th, 10th, 25th, and 40th quantiles, as well as the mean.
    qHat  <-  r.0.Dat$qY.wt.t[1]
    pDens  <-  dpois(x=c(0:10), lambda=(U*x/(h*s)))
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,10), ylim = c(0, 0.7), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(pDens~c(0:10), lwd=1.5, lty=1, col=transparentColor("#252525", opacity=1))
        points(pDens~c(0:10), lwd=1.5, lty=1, pch=21, col=transparentColor("#252525", opacity=1), bg=transparentColor("#252525", opacity=0.25))
        abline(v=(U*x/(h*s)), col=1, lwd=1)
        abline(v=c(0,1), col=2, lwd=1, lty=2)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.5, 1.1, substitute(italic(x)==xx, list(xx=x)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.04, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(paste(Pr(italic(r)~"|"~italic(x)))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel( 0.5,  -0.2,  expression(italic(r)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(r)=="{0,1}"),
                            expression(bar(italic(r))==frac(italic(U*x),italic(h*s)))),
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
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(min(na.omit(r.1.wbar.YI.t)), max(na.omit(r.0.wbar.YI.t))), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(r.0.wbar.YI.t ~ seq_along(r.0.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1))
        lines(r.1.wbar.YI.t ~ seq_along(r.1.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1))
        abline(h=1, lwd=1, lty=2, col=2)
        # Points depicting (arbitrary) extinction frequency of 0.00001
        #points(r.0.wbar.YI.t[length(r.0.wbar.YI.t)] ~ length(r.0.wbar.YI.t), pch=8, col='tomato')
        points(r.1.wbar.YI.t[length(r.1.wbar.YI.t)] ~ length(r.1.wbar.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
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
                            expression(italic(r)==0),
                            expression(italic(r)==1)),
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
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        lines(YI.t[1:length(r.0.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.0.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.1.wbar.YI.t)] ~ seq_along(YI.t[1:length(r.1.wbar.YI.t)]), lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.1.Dat)
        # Points
#        points(YI.t[length(r.0.wbar.YI.t)] ~ length(r.0.wbar.YI.t), pch=8, col='tomato', data=r.0.Dat)
        points(YI.t[length(r.1.wbar.YI.t)] ~ seq_along(YI.t)[length(r.1.wbar.YI.t)], pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.1.Dat)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(paste("Inversion Frequency (", italic(Y[I]), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        


    # reset Colors
     COLS  <-  c('red',
                 '#252525')
## Panel D
# deleterious allele frequencies for r = 5
     # Make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,0.00095), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(XaOv.wt.t ~ seq_along(XaOv.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.0.Dat)
        lines(XaOv.del.t ~ seq_along(XaOv.del.t), lwd=1.5, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.0.Dat)
        lines(qI.wt.t ~ seq_along(qI.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.0.Dat)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(italic(r)==0), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(paste("Deleterious Allele Frequency")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4]*0.98,
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

    # make data for plotting
    r.0.Dat  <-  makeDeterministicFigSimData(r = 0, U = U, x = x, h = h, generations = 3500)
    r.1.Dat  <-  makeDeterministicFigSimData(r = 1, U = U, x = x, h = h, generations = 3500)
    r.2.Dat  <-  makeDeterministicFigSimData(r = 2, U = U, x = x, h = h, generations = 3500)

    r.0.wbar.YI.t  <-  r.0.Dat$wbar.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.1.wbar.YI.t  <-  r.1.Dat$wbar.YI.t[round(r.1.Dat$YI.t,digits=4) > 0]
    r.2.wbar.YI.t  <-  r.2.Dat$wbar.YI.t[round(r.2.Dat$YI.t,digits=4) > 0]

   # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  rev(colfunc(5))

## Panel E
# Distribution of captured deleterious alleles (r) given inversion size (x)
# with benchmarks for 5th, 10th, 25th, and 40th quantiles, as well as the mean.
    qHat  <-  r.0.Dat$qY.wt.t[1]
    rBar  <-  round(nTot*x*qHat, digits=2)
    pDens  <-  dpois(x=c(0:10), lambda=U*x/(h*s))
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,10), ylim = c(0, 0.7), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(pDens~c(0:10), lwd=1.5, lty=1, col=transparentColor("#252525", opacity=1))
        points(pDens~c(0:10), lwd=1.5, lty=1, pch=21, col=transparentColor("#252525", opacity=1), bg=transparentColor("#252525", opacity=0.25))
        abline(v=(U*x/(h*s)), col=1, lwd=1)
        abline(v=c(0:2), col=2, lwd=1, lty=2)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.5, 1.1, substitute(italic(x)==xx, list(xx=x)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.04, 1.05, 'E', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.2,  expression(italic(r)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(r)=="{0,1,2}"),
                            expression(bar(italic(r))==frac(italic(U*x),italic(h*s)))),
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
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(min(na.omit(r.2.wbar.YI.t)), max(na.omit(r.0.wbar.YI.t))), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(r.0.wbar.YI.t ~ seq_along(r.0.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1))
        lines(r.1.wbar.YI.t ~ seq_along(r.1.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1))
        lines(r.2.wbar.YI.t ~ seq_along(r.2.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1))
        abline(h=1, lwd=1, lty=2, col=2)
        # Points depicting (arbitrary) extinction frequency of 0.00001
#        points(r.0.wbar.YI.t[length(r.0.wbar.YI.t)] ~ length(r.0.wbar.YI.t), pch=8, col='tomato')
        points(r.1.wbar.YI.t[length(r.1.wbar.YI.t)] ~ length(r.1.wbar.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        points(r.2.wbar.YI.t[length(r.2.wbar.YI.t)] ~ length(r.2.wbar.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
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
                            expression(italic(r)==0),
                            expression(italic(r)==1),
                            expression(italic(r)==2)),
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
    # make plot
#    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        lines(YI.t[1:length(r.0.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.0.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.1.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.1.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.1.Dat)
        lines(YI.t[1:length(r.2.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.2.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.2.Dat)
        # Points
#        points(YI.t[length(r.0.wbar.YI.t)] ~ length(r.0.wbar.YI.t), pch=8, col='tomato', data=r.0.Dat)
        points(YI.t[length(r.1.wbar.YI.t)] ~ length(r.1.wbar.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.1.Dat)
        points(YI.t[length(r.2.wbar.YI.t)] ~ length(r.2.wbar.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.2.Dat)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'G', cex=1.2, adj=c(0.5, 0.5), xpd=NA)


    # reset Colors
     COLS  <-  c('red',
                 '#252525')
## Panel H
# deleterious allele frequencies for r = 5
     # Make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.1.Dat)), ylim = c(0,0.00095), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(XaOv.wt.t ~ seq_along(XaOv.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(XaOv.del.t ~ seq_along(XaOv.del.t), lwd=1.5, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qI.wt.t ~ seq_along(qI.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'H', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1, expression(italic(r)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        




## 3rd Column of Panels: x = 0.2
    # Set inversion size, r values
    x        <-  0.2

    # make data for plotting
    r.0.Dat  <-  makeDeterministicFigSimData(r = 0, U = U, x = x, h = h, generations = 3500)
    r.1.Dat  <-  makeDeterministicFigSimData(r = 1, U = U, x = x, h = h, generations = 3500)
    r.2.Dat  <-  makeDeterministicFigSimData(r = 2, U = U, x = x, h = h, generations = 3500)
    r.3.Dat  <-  makeDeterministicFigSimData(r = 3, U = U, x = x, h = h, generations = 3500)

    r.0.wbar.YI.t  <-  r.0.Dat$wbar.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.1.wbar.YI.t  <-  r.1.Dat$wbar.YI.t[round(r.1.Dat$YI.t,digits=4) > 0]
    r.2.wbar.YI.t  <-  r.2.Dat$wbar.YI.t[round(r.2.Dat$YI.t,digits=4) > 0]
    r.3.wbar.YI.t  <-  r.3.Dat$wbar.YI.t[round(r.3.Dat$YI.t,digits=4) > 0]

   # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  rev(colfunc(5))

## Panel I
# Distribution of captured deleterious alleles (r) given inversion size (x)
# with benchmarks for 5th, 10th, 25th, and 40th quantiles, as well as the mean.
    qHat  <-  r.0.Dat$qY.wt.t[1]
    rBar  <-  round(nTot*x*qHat, digits=2)
    pDens  <-  dpois(x=c(0:10), lambda=U*x/(h*s))
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,10), ylim = c(0, 0.7), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(pDens~c(0:10), lwd=1.5, lty=1, col=transparentColor("#252525", opacity=1))
        points(pDens~c(0:10), lwd=1.5, lty=1, pch=21, col=transparentColor("#252525", opacity=1), bg=transparentColor("#252525", opacity=0.25))
        abline(v=rBar, col=1, lwd=1)
        abline(v=c(0:3), col=2, lwd=1, lty=2)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.5, 1.1, substitute(italic(x)==xx, list(xx=x)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.04, 1.05, 'I', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.2,  expression(italic(r)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(r)=="{0,1,2,3}"),
                            expression(bar(italic(r))==frac(italic(U*x),italic(h*s)))),
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
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(min(na.omit(r.3.wbar.YI.t)), max(na.omit(r.0.wbar.YI.t))), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(r.0.wbar.YI.t ~ seq_along(r.0.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1))
        lines(r.1.wbar.YI.t ~ seq_along(r.1.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1))
        lines(r.2.wbar.YI.t ~ seq_along(r.2.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.3.wbar.YI.t ~ seq_along(r.3.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1))
        abline(h=1, lwd=1, lty=2, col=2)
        # Points depicting (arbitrary) extinction frequency of 0.00001
#        points(r.0.wbar.YI.t[length(r.0.wbar.YI.t)] ~ length(r.0.wbar.YI.t), pch=8, col='tomato')
        points(r.1.wbar.YI.t[length(r.1.wbar.YI.t)] ~ length(r.1.wbar.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        points(r.2.wbar.YI.t[length(r.2.wbar.YI.t)] ~ length(r.2.wbar.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
        points(r.3.wbar.YI.t[length(r.3.wbar.YI.t)] ~ length(r.3.wbar.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
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
                            expression(italic(r)==0),
                            expression(italic(r)==1),
                            expression(italic(r)==2),
                            expression(italic(r)==3)),
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
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        lines(YI.t[1:length(r.0.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.0.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.1.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.1.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.1.Dat)
        lines(YI.t[1:length(r.2.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.2.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.2.Dat)
        lines(YI.t[1:length(r.3.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.3.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.3.Dat)
        # Points
#        points(YI.t[length(r.0.wbar.YI.t)] ~ length(r.0.wbar.YI.t), pch=8, col='tomato', data=r.0.Dat)
        points(YI.t[length(r.1.wbar.YI.t)] ~ length(r.1.wbar.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.1.Dat)
        points(YI.t[length(r.2.wbar.YI.t)] ~ length(r.2.wbar.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.2.Dat)
        points(YI.t[length(r.3.wbar.YI.t)] ~ length(r.3.wbar.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.3.Dat)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'K', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # reset Colors
     COLS  <-  c('red',
                 '#252525')
## Panel L
# deleterious allele frequencies for r = 1
     # Make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,0.00095), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(XaOv.wt.t ~ seq_along(XaOv.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(XaOv.del.t ~ seq_along(XaOv.del.t), lwd=1.5, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qI.wt.t ~ seq_along(qI.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'L', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(italic(r)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        



## 4th Column of Panels: x = 0.1
    # Set inversion size, r values
    x        <-  0.5

    # make data for plotting
    r.0.Dat  <-  makeDeterministicFigSimData(r = 0, U = U, x = x, h = h, generations = 3500)
    r.1.Dat  <-  makeDeterministicFigSimData(r = 1, U = U, x = x, h = h, generations = 3500)
    r.2.Dat  <-  makeDeterministicFigSimData(r = 2, U = U, x = x, h = h, generations = 3500)
    r.3.Dat  <-  makeDeterministicFigSimData(r = 3, U = U, x = x, h = h, generations = 3500)
    r.4.Dat  <-  makeDeterministicFigSimData(r = 4, U = U, x = x, h = h, generations = 3500)

    r.0.wbar.YI.t   <-  r.0.Dat$wbar.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.1.wbar.YI.t   <-  r.1.Dat$wbar.YI.t[round(r.1.Dat$YI.t,digits=4) > 0]
    r.2.wbar.YI.t   <-  r.2.Dat$wbar.YI.t[round(r.2.Dat$YI.t,digits=4) > 0]
    r.3.wbar.YI.t   <-  r.3.Dat$wbar.YI.t[round(r.3.Dat$YI.t,digits=4) > 0]   
    r.4.wbar.YI.t   <-  r.4.Dat$wbar.YI.t[round(r.4.Dat$YI.t,digits=4) > 0]   

   # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  rev(colfunc(5))

## Panel M
# Distribution of captured deleterious alleles (r) given inversion size (x)
# with benchmarks for 5th, 10th, 25th, and 40th quantiles, as well as the mean.
    qHat  <-  r.0.Dat$qY.wt.t[1]
    rBar  <-  round(nTot*x*qHat, digits=2)
    pDens  <-  dpois(x=c(0:10), lambda=U*x/(h*s))
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,10), ylim = c(0, 0.4), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(pDens~c(0:10), lwd=1.5, lty=1, col=transparentColor("#252525", opacity=1))
        points(pDens~c(0:10), lwd=1.5, lty=1, pch=21, col=transparentColor("#252525", opacity=1), bg=transparentColor("#252525", opacity=0.25))
        abline(v=(U*x/(h*s)), col=1, lwd=1)
        abline(v=c(0,1,2,3,4), col=2, lwd=1, lty=2)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.5, 1.1, substitute(italic(x)==xx, list(xx=x)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.04, 1.05, 'M', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.2,  expression(italic(r)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(r)=="{0,1,2,3,4}"),
                            expression(bar(italic(r))==frac(italic(U*x),italic(h*s)))),
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
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(min(na.omit(r.4.wbar.YI.t)), max(na.omit(r.0.wbar.YI.t))), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(r.0.wbar.YI.t ~ seq_along(r.0.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1))
        lines(r.1.wbar.YI.t ~ seq_along(r.1.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1))
        lines(r.2.wbar.YI.t ~ seq_along(r.2.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.3.wbar.YI.t ~ seq_along(r.3.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r.4.wbar.YI.t ~ seq_along(r.4.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1))
        abline(h=1, lwd=1, lty=2, col=2)
        # Points depicting (arbitrary) extinction frequency of 0.00001
#        points(r.0.wbar.YI.t[length(r.0.wbar.YI.t)] ~ length(r.0.wbar.YI.t), pch=8, col='tomato')
        points(r.1.wbar.YI.t[length(r.1.wbar.YI.t)] ~ length(r.1.wbar.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        points(r.2.wbar.YI.t[length(r.2.wbar.YI.t)] ~ length(r.2.wbar.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
        points(r.3.wbar.YI.t[length(r.3.wbar.YI.t)] ~ length(r.3.wbar.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r.4.wbar.YI.t[length(r.4.wbar.YI.t)] ~ length(r.4.wbar.YI.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7))
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
                            expression(italic(r)==0),
                            expression(italic(r)==1),
                            expression(italic(r)==2),
                            expression(italic(r)==3),
                            expression(italic(r)==4)),
               lty     =  1,
               lwd     =  1.5,
               col     =  c(transparentColor(COLS[5], opacity=1),
                            transparentColor(COLS[4], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[1], opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

## Panel O
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(YI.t[1:length(r.0.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.0.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.1.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.1.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.1.Dat)
        lines(YI.t[1:length(r.2.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.2.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.2.Dat)
        lines(YI.t[1:length(r.3.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.3.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.3.Dat)
        lines(YI.t[1:length(r.4.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.4.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.4.Dat)
        # Points depicting (arbitrary) extinction frequency of 0.00001
#        points(YI.t[length(r.0.wbar.YI.t)] ~ length(r.0.wbar.YI.t), pch=8, col='tomato', data=r.0.Dat)
        points(YI.t[length(r.1.wbar.YI.t)] ~ length(r.1.wbar.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.1.Dat)
        points(YI.t[length(r.2.wbar.YI.t)] ~ length(r.2.wbar.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.2.Dat)
        points(YI.t[length(r.3.wbar.YI.t)] ~ length(r.3.wbar.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.3.Dat)
        points(YI.t[length(r.4.wbar.YI.t)] ~ length(r.4.wbar.YI.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.4.Dat)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'O', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        # Legend

    # reset Colors
     COLS  <-  c('red',
                 '#252525')
## Panel P
# deleterious allele frequencies for r = 1
     # Make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,0.00095), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(XaOv.wt.t ~ seq_along(XaOv.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(XaOv.del.t ~ seq_along(XaOv.del.t), lwd=1.5, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qI.wt.t ~ seq_along(qI.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'P', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(italic(r)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        



## 5th Column of Panels: x = 0.8
    # Set inversion size, r values
    x    <-  0.8

    # make data for plotting
    r.0.Dat  <-  makeDeterministicFigSimData(r = 0, U = U, x = x, h = h, generations = 3500)
    r.1.Dat  <-  makeDeterministicFigSimData(r = 1, U = U, x = x, h = h, generations = 3500)
    r.2.Dat  <-  makeDeterministicFigSimData(r = 2, U = U, x = x, h = h, generations = 3500)
    r.4.Dat  <-  makeDeterministicFigSimData(r = 4, U = U, x = x, h = h, generations = 3500)
    r.6.Dat  <-  makeDeterministicFigSimData(r = 6, U = U, x = x, h = h, generations = 3500)

    r.0.wbar.YI.t  <-  r.0.Dat$wbar.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.1.wbar.YI.t  <-  r.1.Dat$wbar.YI.t[round(r.1.Dat$YI.t,digits=4) > 0]
    r.2.wbar.YI.t  <-  r.2.Dat$wbar.YI.t[round(r.2.Dat$YI.t,digits=4) > 0]
    r.4.wbar.YI.t  <-  r.4.Dat$wbar.YI.t[round(r.4.Dat$YI.t,digits=4) > 0]   
    r.6.wbar.YI.t  <-  r.6.Dat$wbar.YI.t[round(r.6.Dat$YI.t,digits=4) > 0]   

   # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  rev(colfunc(5))

## Panel Q
# Distribution of captured deleterious alleles (r) given inversion size (x)
# with benchmarks for 5th, 10th, 25th, and 40th quantiles, as well as the mean.
    qHat  <-  r.0.Dat$qY.wt.t[1]
    rBar  <-  round(nTot*x*qHat, digits=2)
    pDens  <-  dpois(x=c(0:20), lambda=U*x/(h*s))
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,20), ylim = c(0, 0.7), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(pDens~c(0:20), lwd=1.5, lty=1, col=transparentColor("#252525", opacity=1))
        points(pDens~c(0:20), lwd=1.5, lty=1, pch=21, col=transparentColor("#252525", opacity=1), bg=transparentColor("#252525", opacity=0.25))
        abline(v=(U*x/(h*s)), col=1, lwd=1)
        abline(v=c(0,1,2,4,6), col=2, lwd=1, lty=2)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.5, 1.1, substitute(italic(x)==xx, list(xx=x)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.04, 1.05, 'Q', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.2,  expression(italic(r)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(r)=="{0,1,2,4,6}"),
                            expression(bar(italic(r))==frac(italic(U*x),italic(h*s)))),
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
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(min(na.omit(r.6.wbar.YI.t)), max(na.omit(r.0.wbar.YI.t))), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(r.0.wbar.YI.t ~ seq_along(r.0.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1))
        lines(r.1.wbar.YI.t ~ seq_along(r.1.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1))
        lines(r.2.wbar.YI.t ~ seq_along(r.2.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.4.wbar.YI.t ~ seq_along(r.4.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r.6.wbar.YI.t ~ seq_along(r.6.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1))
        abline(h=1, lwd=1, lty=2, col=2)
        # Points depicting (arbitrary) extinction frequency of 0.00001
#        points(r.0.wbar.YI.t[length(r.0.wbar.YI.t)] ~ length(r.0.wbar.YI.t), pch=8, col='tomato')
        points(r.1.wbar.YI.t[length(r.1.wbar.YI.t)] ~ length(r.1.wbar.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        points(r.2.wbar.YI.t[length(r.2.wbar.YI.t)] ~ length(r.2.wbar.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
        points(r.4.wbar.YI.t[length(r.4.wbar.YI.t)] ~ length(r.4.wbar.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r.6.wbar.YI.t[length(r.6.wbar.YI.t)] ~ length(r.6.wbar.YI.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7))
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
                            expression(italic(r)==0),
                            expression(italic(r)==1),
                            expression(italic(r)==2),
                            expression(italic(r)==4),
                            expression(italic(r)==6)),
               lty     =  1,
               lwd     =  1.5,
               col     =  c(transparentColor(COLS[5], opacity=1),
                            transparentColor(COLS[4], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[1], opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

## Panel S
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(YI.t[1:length(r.0.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.0.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.1.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.1.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.1.Dat)
        lines(YI.t[1:length(r.2.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.2.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.2.Dat)
        lines(YI.t[1:length(r.4.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.4.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.4.Dat)
        lines(YI.t[1:length(r.6.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.6.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.6.Dat)
        # Points depicting (arbitrary) extinction frequency of 0.00001
#        points(YI.t[length(r.0.wbar.YI.t)] ~ length(r.0.wbar.YI.t), pch=8, col='tomato', data=r.0.Dat)
        points(YI.t[length(r.1.wbar.YI.t)] ~ length(r.1.wbar.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.1.Dat)
        points(YI.t[length(r.2.wbar.YI.t)] ~ length(r.2.wbar.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.2.Dat)
        points(YI.t[length(r.4.wbar.YI.t)] ~ length(r.4.wbar.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.4.Dat)
        points(YI.t[length(r.6.wbar.YI.t)] ~ length(r.6.wbar.YI.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.6.Dat)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'S', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        # Legend

    # reset Colors
     COLS  <-  c('red',
                 '#252525')
## Panel T
# deleterious allele frequencies for r = 1
     # Make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.1.Dat)), ylim = c(0,0.00095), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(XaOv.wt.t ~ seq_along(XaOv.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(XaOv.del.t ~ seq_along(XaOv.del.t), lwd=1.5, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qI.wt.t ~ seq_along(qI.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'T', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(italic(r)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        

}


########################################
#' Deterministic dynamics of inversions expanding the SDR on a Y chromosome:
#' Illustration of: (A) inversion relative fitness
#'                  (B) Inversion frequency dynamics, and 
#'                  (C) del. allele frequency dynamics.
deterministicSuppFig_h0.01  <-  function() {


    # set plot layout
    layout.mat <- matrix(c(1:20), nrow=4, ncol=5, byrow=FALSE)
    layout     <- layout(layout.mat,respect=TRUE)

    # Make data for plotting (x = 0.1; 1st column of panels)
    h        <-  0.01
    s        <-  0.01
    U        <-  0.02
    Ufactor  <-  2
    nTot     <-  10^4 
    u        <-  U/nTot
    x        <-  0.05

    # make data for plotting
    r.0.Dat  <-  makeDeterministicFigSimData(r = 0, U = U, x = x, h = h, generations = 20000)
    r.1.Dat  <-  makeDeterministicFigSimData(r = 1, U = U, x = x, h = h, generations = 20000)
    r.2.Dat  <-  makeDeterministicFigSimData(r = 2, U = U, x = x, h = h, generations = 20000)
    r.3.Dat  <-  makeDeterministicFigSimData(r = 3, U = U, x = x, h = h, generations = 20000)
    r.4.Dat  <-  makeDeterministicFigSimData(r = 4, U = U, x = x, h = h, generations = 20000)

    r.0.wbar.YI.t   <-  r.0.Dat$wbar.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.1.wbar.YI.t   <-  r.1.Dat$wbar.YI.t[round(r.1.Dat$YI.t,digits=4) > 0]
    r.2.wbar.YI.t   <-  r.2.Dat$wbar.YI.t[round(r.2.Dat$YI.t,digits=4) > 0]
    r.3.wbar.YI.t   <-  r.3.Dat$wbar.YI.t[round(r.3.Dat$YI.t,digits=4) > 0]   
    r.4.wbar.YI.t   <-  r.4.Dat$wbar.YI.t[round(r.4.Dat$YI.t,digits=4) > 0]   

    qHat  <-  r.0.Dat$qY.wt.t[1]
    rBar  <-  round(nTot*x*qHat)

   # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  rev(colfunc(5))

## Panel A
# Distribution of captured deleterious alleles (r) given inversion size (x)
# with benchmarks for 5th, 10th, 25th, and 40th quantiles, as well as the mean.
   pDens  <-  dpois(x=c(0:15), lambda=rBar)
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,20), ylim = c(0, 0.2), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(pDens~c(0:15), lwd=1.5, lty=1, col=transparentColor("#252525", opacity=1))
        points(pDens~c(0:15), lwd=1.5, lty=1, pch=21, col=transparentColor("#252525", opacity=1), bg=transparentColor("#252525", opacity=0.25))
        abline(v=rBar, col=1, lwd=1)
        abline(v=c(0,1,2,3,4), col=2, lwd=1, lty=2)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.5, 1.1, substitute(italic(x)==xx, list(xx=x)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.04, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(paste(Pr(italic(r)~"|"~italic(x)))), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)
        proportionalLabel( 0.5,  -0.2,  expression(italic(r)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(r)=="{0,1,2,3,4}"),
                            expression(bar(italic(r))==frac(italic(U*x),italic(h*s)))),
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
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(min(na.omit(r.4.wbar.YI.t)), max(na.omit(r.0.wbar.YI.t))), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(r.0.wbar.YI.t ~ seq_along(r.0.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1))
        lines(r.1.wbar.YI.t ~ seq_along(r.1.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1))
        lines(r.2.wbar.YI.t ~ seq_along(r.2.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.3.wbar.YI.t ~ seq_along(r.3.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r.4.wbar.YI.t ~ seq_along(r.4.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1))
        abline(h=1, lwd=1, lty=2, col=2)
        # Points depicting (arbitrary) extinction frequency of 0.00001
#        points(r.0.wbar.YI.t[length(r.0.wbar.YI.t)] ~ length(r.0.wbar.YI.t), pch=8, col='tomato')
        points(r.1.wbar.YI.t[length(r.1.wbar.YI.t)] ~ length(r.1.wbar.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        points(r.2.wbar.YI.t[length(r.2.wbar.YI.t)] ~ length(r.2.wbar.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
        points(r.3.wbar.YI.t[length(r.3.wbar.YI.t)] ~ length(r.3.wbar.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r.4.wbar.YI.t[length(r.4.wbar.YI.t)] ~ length(r.4.wbar.YI.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7))
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
                            expression(italic(r)==0),
                            expression(italic(r)==1),
                            expression(italic(r)==2),
                            expression(italic(r)==3),
                            expression(italic(r)==4)),
               lty     =  1,
               lwd     =  1.5,
               col     =  c(transparentColor(COLS[5], opacity=1),
                            transparentColor(COLS[4], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[1], opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

## Panel C
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(YI.t[1:length(r.0.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.0.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.1.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.1.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.1.Dat)
        lines(YI.t[1:length(r.2.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.2.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.2.Dat)
        lines(YI.t[1:length(r.3.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.3.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.3.Dat)
        lines(YI.t[1:length(r.4.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.4.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.4.Dat)
        # Points
#        points(YI.t[length(r.0.wbar.YI.t)] ~ length(r.0.wbar.YI.t), pch=8, col='tomato', data=r.0.Dat)
        points(YI.t[length(r.1.wbar.YI.t)] ~ seq_along(YI.t)[length(r.1.wbar.YI.t)], pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.1.Dat)
        points(YI.t[length(r.2.wbar.YI.t)] ~ seq_along(YI.t)[length(r.2.wbar.YI.t)], pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.2.Dat)
        points(YI.t[length(r.3.wbar.YI.t)] ~ seq_along(YI.t)[length(r.3.wbar.YI.t)], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.3.Dat)
        points(YI.t[length(r.4.wbar.YI.t)] ~ seq_along(YI.t)[length(r.4.wbar.YI.t)], pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.4.Dat)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(paste("Inversion Frequency (", italic(Y[I]), ")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        # Legend

    # reset Colors
     COLS  <-  c('red',
                 '#252525')
## Panel D
# deleterious allele frequencies for r = 2
     # Make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,qHat*1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(XaOv.wt.t ~ seq_along(XaOv.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(XaOv.del.t ~ seq_along(XaOv.del.t), lwd=1.5, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.1.Dat)
        lines(qI.wt.t ~ seq_along(qI.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'D', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(italic(r)==2), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.35,  0.5,   expression(paste("Deleterious Allele Frequency")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4]*0.9,
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

    # make data for plotting
    r.0.Dat  <-  makeDeterministicFigSimData(r = 0, U = U, x = x, h = h, generations = 20000)
    r.2.Dat  <-  makeDeterministicFigSimData(r = 2, U = U, x = x, h = h, generations = 20000)
    r.4.Dat  <-  makeDeterministicFigSimData(r = 4, U = U, x = x, h = h, generations = 20000)
    r.6.Dat  <-  makeDeterministicFigSimData(r = 6, U = U, x = x, h = h, generations = 20000)
    r.8.Dat  <-  makeDeterministicFigSimData(r = 8, U = U, x = x, h = h, generations = 20000)

    r.0.wbar.YI.t  <-  r.0.Dat$wbar.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.2.wbar.YI.t  <-  r.2.Dat$wbar.YI.t[round(r.2.Dat$YI.t,digits=4) > 0]
    r.4.wbar.YI.t  <-  r.4.Dat$wbar.YI.t[round(r.4.Dat$YI.t,digits=4) > 0]
    r.6.wbar.YI.t  <-  r.6.Dat$wbar.YI.t[round(r.6.Dat$YI.t,digits=4) > 0]
    r.8.wbar.YI.t  <-  r.8.Dat$wbar.YI.t[round(r.8.Dat$YI.t,digits=4) > 0]

   # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  rev(colfunc(5))

    qHat  <-  r.0.Dat$qY.wt.t[1]
    rBar  <-  round(nTot*x*qHat)

## Panel E
# Distribution of captured deleterious alleles (r) given inversion size (x)
   pDens  <-  dpois(x=c(0:25), lambda=rBar)
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,25), ylim = c(0, 0.2), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(pDens~c(0:25), lwd=1.5, lty=1, col=transparentColor("#252525", opacity=1))
        points(pDens~c(0:25), lwd=1.5, lty=1, pch=21, col=transparentColor("#252525", opacity=1), bg=transparentColor("#252525", opacity=0.25))
        abline(v=rBar, col=1, lwd=1)
        abline(v=c(0,2,4,6,8), col=2, lwd=1, lty=2)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.5, 1.1, substitute(italic(x)==xx, list(xx=x)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.04, 1.05, 'E', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.2,  expression(italic(r)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(r)=="{0,2,4,6,8}"),
                            expression(bar(italic(r))==frac(italic(U*x),italic(h*s)))),
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
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(min(na.omit(r.8.wbar.YI.t)), max(na.omit(r.0.wbar.YI.t))), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(r.0.wbar.YI.t ~ seq_along(r.0.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1))
        lines(r.2.wbar.YI.t ~ seq_along(r.2.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1))
        lines(r.4.wbar.YI.t ~ seq_along(r.4.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.6.wbar.YI.t ~ seq_along(r.6.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r.8.wbar.YI.t ~ seq_along(r.8.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1))
        abline(h=1, lwd=1, lty=2, col=2)
        # Points depicting (arbitrary) extinction frequency of 0.00001
#        points(r.0.wbar.YI.t[length(r.0.wbar.YI.t)] ~ length(r.0.wbar.YI.t), pch=8, col='tomato')
#        points(r.2.wbar.YI.t[length(r.2.wbar.YI.t)] ~ length(r.2.wbar.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
#        points(r.4.wbar.YI.t[length(r.4.wbar.YI.t)] ~ length(r.4.wbar.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
        points(r.6.wbar.YI.t[length(r.6.wbar.YI.t)] ~ length(r.6.wbar.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r.8.wbar.YI.t[length(r.8.wbar.YI.t)] ~ length(r.8.wbar.YI.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7))
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
                            expression(italic(r)==0),
                            expression(italic(r)==2),
                            expression(italic(r)==4),
                            expression(italic(r)==6),
                            expression(italic(r)==8)),
               lty     =  1,
               lwd     =  1.5,
               col     =  c(transparentColor(COLS[5], opacity=1),
                            transparentColor(COLS[4], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[1], opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )


## Panel G
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(YI.t[1:length(r.0.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.0.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.2.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.2.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.2.Dat)
        lines(YI.t[1:length(r.4.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.4.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.4.Dat)
        lines(YI.t[1:length(r.6.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.6.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.6.Dat)
        lines(YI.t[1:length(r.8.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.8.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.8.Dat)
        # Points
#        points(YI.t[length(r.0.wbar.YI.t)] ~ length(r.0.wbar.YI.t), pch=8, col='tomato', data=r.0.Dat)
#        points(YI.t[length(r.2.wbar.YI.t)] ~ seq_along(YI.t)[length(r.2.wbar.YI.t)], pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.2.Dat)
#        points(YI.t[length(r.4.wbar.YI.t)] ~ seq_along(YI.t)[length(r.4.wbar.YI.t)], pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.4.Dat)
        points(YI.t[length(r.6.wbar.YI.t)] ~ seq_along(YI.t)[length(r.6.wbar.YI.t)], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.6.Dat)
        points(YI.t[length(r.8.wbar.YI.t)] ~ seq_along(YI.t)[length(r.8.wbar.YI.t)], pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.8.Dat)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'G', cex=1.2, adj=c(0.5, 0.5), xpd=NA)


    # reset Colors
     COLS  <-  c('red',
                 '#252525')
## Panel H
# deleterious allele frequencies for r = 5
     # Make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,qHat*1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(XaOv.wt.t ~ seq_along(XaOv.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.2.Dat)
        lines(XaOv.del.t ~ seq_along(XaOv.del.t), lwd=1.5, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.2.Dat)
        lines(qI.wt.t ~ seq_along(qI.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.2.Dat)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'H', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1, expression(italic(r)==1), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        




## 3rd Column of Panels: x = 0.2
    # Set inversion size, r values
    x        <-  0.2

    # make data for plotting
    r.0.Dat   <-  makeDeterministicFigSimData(r = 0, U = U,  x = x, h = h, generations = 20000)
    r.5.Dat   <-  makeDeterministicFigSimData(r = 5, U = U,  x = x, h = h, generations = 20000)
    r.10.Dat  <-  makeDeterministicFigSimData(r = 10, U = U, x = x, h = h, generations = 20000)
    r.15.Dat  <-  makeDeterministicFigSimData(r = 15, U = U, x = x, h = h, generations = 20000)
    r.20.Dat  <-  makeDeterministicFigSimData(r = 20, U = U, x = x, h = h, generations = 20000)

    r.0.wbar.YI.t   <-  r.0.Dat$wbar.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.5.wbar.YI.t   <-  r.5.Dat$wbar.YI.t[round(r.5.Dat$YI.t,digits=4) > 0]
    r.10.wbar.YI.t  <-  r.10.Dat$wbar.YI.t[round(r.10.Dat$YI.t,digits=4) > 0]
    r.15.wbar.YI.t  <-  r.15.Dat$wbar.YI.t[round(r.15.Dat$YI.t,digits=4) > 0]
    r.20.wbar.YI.t  <-  r.20.Dat$wbar.YI.t[round(r.20.Dat$YI.t,digits=4) > 0]

   # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  rev(colfunc(5))

    qHat  <-  r.0.Dat$qY.wt.t[1]
    rBar  <-  round(nTot*x*qHat)

## Panel I
# Distribution of captured deleterious alleles (r) given inversion size (x)
     pDens  <-  dpois(x=c(0:45), lambda=rBar)
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,45), ylim = c(0, 0.4), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(pDens~c(0:45), lwd=1.5, lty=1, col=transparentColor("#252525", opacity=1))
        points(pDens~c(0:45), lwd=1.5, lty=1, pch=21, col=transparentColor("#252525", opacity=1), bg=transparentColor("#252525", opacity=0.25))
        abline(v=rBar, col=1, lwd=1)
        abline(v=c(0,5,10,15,20), col=2, lwd=1, lty=2)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.5, 1.1, substitute(italic(x)==xx, list(xx=x)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.04, 1.05, 'I', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.2,  expression(italic(r)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(r)=="{0,5,10,15,20}"),
                            expression(bar(italic(r))==frac(italic(U*x),italic(h*s)))),
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
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(min(na.omit(r.20.wbar.YI.t)), max(na.omit(r.0.wbar.YI.t))), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(r.0.wbar.YI.t ~ seq_along(r.0.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1))
        lines(r.5.wbar.YI.t ~ seq_along(r.5.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1))
        lines(r.10.wbar.YI.t ~ seq_along(r.10.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.15.wbar.YI.t ~ seq_along(r.15.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r.20.wbar.YI.t ~ seq_along(r.20.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1))
        abline(h=1, lwd=1, lty=2, col=2)
        # Points depicting (arbitrary) extinction frequency of 0.00001
        points(r.0.wbar.YI.t[length(r.0.wbar.YI.t)] ~ length(r.0.wbar.YI.t), pch=8, col='tomato')
#        points(r.5.wbar.YI.t[length(r.5.wbar.YI.t)] ~ length(r.5.wbar.YI.t), pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7))
        points(r.10.wbar.YI.t[length(r.10.wbar.YI.t)] ~ length(r.10.wbar.YI.t), pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7))
        points(r.15.wbar.YI.t[length(r.15.wbar.YI.t)] ~ length(r.15.wbar.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r.20.wbar.YI.t[length(r.20.wbar.YI.t)] ~ length(r.20.wbar.YI.t), pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7))
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
                            expression(italic(r)==0),
                            expression(italic(r)==5),
                            expression(italic(r)==10),
                            expression(italic(r)==15),
                            expression(italic(r)==20)),
               lty     =  1,
               lwd     =  1.5,
               col     =  c(transparentColor(COLS[5], opacity=1),
                            transparentColor(COLS[4], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[1], opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

## Panel K
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        lines(YI.t[1:length(r.0.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.0.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.5.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.5.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.5.Dat)
        lines(YI.t[1:length(r.10.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.10.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.10.Dat)
        lines(YI.t[1:length(r.15.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.15.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.15.Dat)
        lines(YI.t[1:length(r.20.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.20.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.20.Dat)
        # Points
        points(YI.t[length(r.0.wbar.YI.t)] ~ length(r.0.wbar.YI.t), pch=8, col='tomato', data=r.0.Dat)
#        points(YI.t[length(r.5.wbar.YI.t)] ~ seq_along(YI.t)[length(r.5.wbar.YI.t)], pch=21, col=transparentColor(COLS[4], opacity=1), bg=transparentColor(COLS[4], opacity=0.7), data=r.5.Dat)
        points(YI.t[length(r.10.wbar.YI.t)] ~ seq_along(YI.t)[length(r.10.wbar.YI.t)], pch=21, col=transparentColor(COLS[3], opacity=1), bg=transparentColor(COLS[3], opacity=0.7), data=r.10.Dat)
        points(YI.t[length(r.15.wbar.YI.t)] ~ seq_along(YI.t)[length(r.15.wbar.YI.t)], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.15.Dat)
        points(YI.t[length(r.20.wbar.YI.t)] ~ seq_along(YI.t)[length(r.20.wbar.YI.t)], pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.20.Dat)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'K', cex=1.2, adj=c(0.5, 0.5), xpd=NA)

    # reset Colors
     COLS  <-  c('red',
                 '#252525')
## Panel L
# deleterious allele frequencies for r = 5
     # Make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,qHat*1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(XaOv.wt.t ~ seq_along(XaOv.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.10.Dat)
        lines(XaOv.del.t ~ seq_along(XaOv.del.t), lwd=1.5, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.10.Dat)
        lines(qI.wt.t ~ seq_along(qI.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.10.Dat)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'L', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(italic(r)==10), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        



## 4th Column of Panels: x = 0.5
    # Set inversion size, r values
    x        <-  0.5

    # make data for plotting
    r.0.Dat  <-  makeDeterministicFigSimData(r = 0, U = U, x = x, h = h, generations = 20000)
    r.10.Dat  <-  makeDeterministicFigSimData(r = 10, U = U, x = x, h = h, generations = 20000)
    r.20.Dat  <-  makeDeterministicFigSimData(r = 20, U = U, x = x, h = h, generations = 20000)
    r.30.Dat  <-  makeDeterministicFigSimData(r = 30, U = U, x = x, h = h, generations = 20000)
    r.40.Dat  <-  makeDeterministicFigSimData(r = 40, U = U, x = x, h = h, generations = 20000)

    r.0.wbar.YI.t   <-  r.0.Dat$wbar.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.10.wbar.YI.t   <-  r.10.Dat$wbar.YI.t[round(r.10.Dat$YI.t,digits=4) < 1]
    r.20.wbar.YI.t   <-  r.20.Dat$wbar.YI.t[round(r.20.Dat$YI.t,digits=4) < 1]
    r.30.wbar.YI.t   <-  r.30.Dat$wbar.YI.t[round(r.30.Dat$YI.t,digits=4) > 0]   
    r.40.wbar.YI.t   <-  r.40.Dat$wbar.YI.t[round(r.40.Dat$YI.t,digits=4) > 0]   

   # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  rev(colfunc(5))

    qHat  <-  r.0.Dat$qY.wt.t[1]
    rBar  <-  round(nTot*x*qHat)

## Panel M
# Distribution of captured deleterious alleles (r) given inversion size (x)
     pDens  <-  dpois(x=seq(0,100, by=2), lambda=rBar)
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,100), ylim = c(0, 0.2), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(pDens~seq(0,100, by=2), lwd=1.5, lty=1, col=transparentColor("#252525", opacity=1))
        points(pDens~seq(0,100, by=2), lwd=1.5, lty=1, pch=21, col=transparentColor("#252525", opacity=1), bg=transparentColor("#252525", opacity=0.25))
        abline(v=rBar, col=1, lwd=1)
        abline(v=c(0,10,20,30,40), col=2, lwd=1, lty=2)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.5, 1.1, substitute(italic(x)==xx, list(xx=x)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.04, 1.05, 'M', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.2,  expression(italic(r)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(r)=="{0,10,20,30,40}"),
                            expression(bar(italic(r))==frac(italic(U*x),italic(h*s)))),
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
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(min(na.omit(r.40.wbar.YI.t)), max(na.omit(r.0.wbar.YI.t))), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(r.0.wbar.YI.t ~ seq_along(r.0.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1))
        lines(r.10.wbar.YI.t ~ seq_along(r.10.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1))
        lines(r.20.wbar.YI.t ~ seq_along(r.20.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.30.wbar.YI.t ~ seq_along(r.30.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r.40.wbar.YI.t ~ seq_along(r.40.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1))
        abline(h=1, lwd=1, lty=2, col=2)
        # Points depicting (arbitrary) extinction frequency of 0.00001
        points(r.0.wbar.YI.t[length(r.0.wbar.YI.t)] ~ length(r.0.wbar.YI.t), pch=8, col='tomato')
        points(r.10.wbar.YI.t[length(r.10.wbar.YI.t)] ~ length(r.10.wbar.YI.t), pch=8, col='tomato')
        points(r.20.wbar.YI.t[length(r.20.wbar.YI.t)] ~ length(r.20.wbar.YI.t), pch=8, col='tomato')
        points(r.30.wbar.YI.t[length(r.30.wbar.YI.t)] ~ length(r.30.wbar.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r.40.wbar.YI.t[length(r.40.wbar.YI.t)] ~ length(r.40.wbar.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[1], opacity=0.7))
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
                            expression(italic(r)==0),
                            expression(italic(r)==10),
                            expression(italic(r)==20),
                            expression(italic(r)==30),
                            expression(italic(r)==40)),
               lty     =  1,
               lwd     =  1.5,
               col     =  c(transparentColor(COLS[5], opacity=1),
                            transparentColor(COLS[4], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[1], opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

## Panel O
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        lines(YI.t[1:length(r.0.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.0.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.10.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.10.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.10.Dat)
        lines(YI.t[1:length(r.20.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.20.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.20.Dat)
        lines(YI.t[1:length(r.30.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.30.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.30.Dat)
        lines(YI.t[1:length(r.40.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.40.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.40.Dat)
        # Points
#        points(YI.t[length(r.0.wbar.YI.t)] ~ length(r.0.wbar.YI.t), pch=8, col='tomato', data=r.0.Dat)
        points(YI.t[length(r.10.wbar.YI.t)] ~ seq_along(YI.t)[length(r.10.wbar.YI.t)], pch=8, col='tomato', data=r.10.Dat)
        points(YI.t[length(r.20.wbar.YI.t)] ~ seq_along(YI.t)[length(r.20.wbar.YI.t)], pch=8, col='tomato', data=r.20.Dat)
        points(YI.t[length(r.30.wbar.YI.t)] ~ seq_along(YI.t)[length(r.30.wbar.YI.t)], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.30.Dat)
        points(YI.t[length(r.40.wbar.YI.t)] ~ seq_along(YI.t)[length(r.40.wbar.YI.t)], pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.40.Dat)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'O', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        # Legend

    # reset Colors
     COLS  <-  c('red',
                 '#252525')
## Panel P
# deleterious allele frequencies for r = 5
     # Make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,qHat*1.1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(XaOv.wt.t ~ seq_along(XaOv.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.30.Dat)
        lines(XaOv.del.t ~ seq_along(XaOv.del.t), lwd=1.5, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.30.Dat)
        lines(qI.wt.t ~ seq_along(qI.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.30.Dat)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'P', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(italic(r)==30), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        



## 5th Column of Panels: x = 0.8
    # Set inversion size, r values
    x    <-  0.8

    # make data for plotting
    r.0.Dat   <-  makeDeterministicFigSimData(r = 0, U = U,  x = x, h = h, generations = 20000)
    r.20.Dat  <-  makeDeterministicFigSimData(r = 20, U = U, x = x, h = h, generations = 20000)
    r.40.Dat  <-  makeDeterministicFigSimData(r = 40, U = U, x = x, h = h, generations = 20000)
    r.60.Dat  <-  makeDeterministicFigSimData(r = 60, U = U, x = x, h = h, generations = 20000)
    r.80.Dat  <-  makeDeterministicFigSimData(r = 80, U = U, x = x, h = h, generations = 20000)

    r.0.wbar.YI.t  <-  r.0.Dat$wbar.YI.t[round(r.0.Dat$YI.t,digits=4) < 1]
    r.20.wbar.YI.t  <-  r.20.Dat$wbar.YI.t[round(r.20.Dat$YI.t,digits=4) < 1]
    r.40.wbar.YI.t  <-  r.40.Dat$wbar.YI.t[round(r.40.Dat$YI.t,digits=4) < 1]
    r.60.wbar.YI.t  <-  r.60.Dat$wbar.YI.t[round(r.60.Dat$YI.t,digits=4) > 0]   
    r.80.wbar.YI.t <-  r.80.Dat$wbar.YI.t[round(r.80.Dat$YI.t,digits=4) > 0]   

   # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  rev(colfunc(5))

    qHat  <-  r.0.Dat$qY.wt.t[1]
    rBar  <-  round(nTot*x*qHat)

## Panel Q
# Distribution of captured deleterious alleles (r) given inversion size (x)
    pDens  <-  dpois(x=seq(0,160, by=4), lambda=rBar)
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,160), ylim = c(0, 0.4), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(pDens~seq(0,160, by=4), lwd=1.5, lty=1, col=transparentColor("#252525", opacity=1))
        points(pDens~seq(0,160, by=4), lwd=1.5, lty=1, pch=21, col=transparentColor("#252525", opacity=1), bg=transparentColor("#252525", opacity=0.25))
        abline(v=rBar, col=1, lwd=1)
        abline(v=c(0,20,40,60,80), col=2, lwd=1, lty=2)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.5, 1.1, substitute(italic(x)==xx, list(xx=x)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(0.04, 1.05, 'Q', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.2,  expression(italic(r)), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4],
               legend  =  c(
                            expression(italic(r)=="{0,20,40,60,80}"),
                            expression(bar(italic(r))==frac(italic(U*x),italic(h*s)))),
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
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(min(na.omit(r.60.wbar.YI.t)), max(na.omit(r.0.wbar.YI.t))), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # lines
        lines(r.0.wbar.YI.t ~ seq_along(r.0.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1))
        lines(r.20.wbar.YI.t ~ seq_along(r.20.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1))
        lines(r.40.wbar.YI.t ~ seq_along(r.40.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1))
        lines(r.60.wbar.YI.t ~ seq_along(r.60.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1))
        lines(r.80.wbar.YI.t ~ seq_along(r.80.wbar.YI.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1))
        abline(h=1, lwd=1, lty=2, col=2)
        # Points depicting (arbitrary) extinction frequency of 0.00001
        points(r.0.wbar.YI.t[length(r.0.wbar.YI.t)] ~ length(r.0.wbar.YI.t), pch=8, col='tomato')
        points(r.20.wbar.YI.t[length(r.20.wbar.YI.t)] ~ length(r.20.wbar.YI.t), pch=8, col='tomato')
        points(r.40.wbar.YI.t[length(r.40.wbar.YI.t)] ~ length(r.40.wbar.YI.t), pch=8, col='tomato')
        points(r.60.wbar.YI.t[length(r.60.wbar.YI.t)] ~ length(r.60.wbar.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7))
        points(r.80.wbar.YI.t[length(r.80.wbar.YI.t)] ~ length(r.80.wbar.YI.t), pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[1], opacity=0.7))
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
                            expression(italic(r)==0),
                            expression(italic(r)==20),
                            expression(italic(r)==40),
                            expression(italic(r)==60),
                            expression(italic(r)==80)),
               lty     =  1,
               lwd     =  1.5,
               col     =  c(transparentColor(COLS[5], opacity=1),
                            transparentColor(COLS[4], opacity=1),
                            transparentColor(COLS[3], opacity=1),
                            transparentColor(COLS[2], opacity=1),
                            transparentColor(COLS[1], opacity=1)),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )

## Panel S
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.0.Dat)), ylim = c(0,1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        lines(YI.t[1:length(r.0.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.0.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[5], opacity=1), data=r.0.Dat)
        lines(YI.t[1:length(r.20.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.20.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[4], opacity=1), data=r.20.Dat)
        lines(YI.t[1:length(r.40.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.40.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[3], opacity=1), data=r.40.Dat)
        lines(YI.t[1:length(r.60.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.60.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.60.Dat)
        lines(YI.t[1:length(r.80.wbar.YI.t)] ~ seq_along(YI.t)[1:length(r.80.wbar.YI.t)], lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.80.Dat)
        # Points
#        points(YI.t[length(r.0.wbar.YI.t)] ~ length(r.0.wbar.YI.t), pch=8, col='tomato', data=r.0.Dat)
        points(YI.t[length(r.20.wbar.YI.t)] ~ seq_along(YI.t)[length(r.20.wbar.YI.t)], pch=8, col='tomato', data=r.20.Dat)
        points(YI.t[length(r.40.wbar.YI.t)] ~ seq_along(YI.t)[length(r.40.wbar.YI.t)], pch=8, col='tomato', data=r.40.Dat)
        points(YI.t[length(r.60.wbar.YI.t)] ~ seq_along(YI.t)[length(r.60.wbar.YI.t)], pch=21, col=transparentColor(COLS[2], opacity=1), bg=transparentColor(COLS[2], opacity=0.7), data=r.60.Dat)
        points(YI.t[length(r.80.wbar.YI.t)] ~ seq_along(YI.t)[length(r.80.wbar.YI.t)], pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.7), data=r.80.Dat)
        # axes
        axis(1, las=1, labels=NA)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'S', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        # Legend

    # reset Colors
     COLS  <-  c('red',
                 '#252525')
## Panel T
# deleterious allele frequencies for r = 5
     # Make plot
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.4.Dat)), ylim = c(0,qHat*1.05), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(XaOv.wt.t ~ seq_along(XaOv.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1), data=r.60.Dat)
        lines(XaOv.del.t ~ seq_along(XaOv.del.t), lwd=1.5, lty=3, col=transparentColor(COLS[2], opacity=1), data=r.60.Dat)
        lines(qI.wt.t ~ seq_along(qI.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.60.Dat)
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'T', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.1,   expression(italic(r)==60), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.25,  expression(paste("Generations")), cex=1.2, adj=c(0.5, 0.5), xpd=NA)        

}








# FIGURE 1. Illustration of deterministic dynamics
#           for different dominance values
delMutAccumulationAutoVsSLR  <-  function() {

    # set plot layout
    layout.mat <- matrix(c(1:3), nrow=1, ncol=3, byrow=FALSE)
    layout     <- layout(layout.mat,respect=TRUE)

    # Make data for plotting 
    s        <-  0.01
    h        <-  0.25
    Ufactor  <-  2
    nTot     <-  10^4 
    U        <-  Ufactor*s
    u        <-  U/nTot
    x        <-  0.2


    # reset Colors
    ###########
    # h = 0.25

    # make data for plotting
    r.1.Dat  <-  makeDeterministicFigSimData(r = 1, x = x, h = h, generations = 2500)
    r.1.wbar.YI.t  <-  r.1.Dat$wbar.YI.t[round(r.1.Dat$YI.t,digits=4) > 0]
    autotime  <-  c(1:2500)
    q.wt.Auto  <-  p.t(u=10^-6, sdHom=0.01, h=h, t=autotime) 

   # Colors
     COLS  <-  c('red','#252525')

    # Axes limits
    relFitMin  <-  0.995
    relFitMax  <-  1.015


## Panel A
# Effect of initial del. mut. load on inversion relative fitness 
    # make plot
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')

# deleterious allele frequencies for r = 1
     # Make plot
    delMax  <-  r.1.Dat$XaOv.wt.t[1]*1.05
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.1.Dat)), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(qI.wt.t/max(qI.wt.t) ~ seq_along(qI.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
        lines(q.wt.Auto/max(q.wt.Auto) ~ autotime, lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1))
        # axes
        axis(1, las=1)
        axis(2, las=1, cex=1.2)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'A', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(paste(italic(h)==0.25,", ", italic(d)==1)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel(-0.2,  0.5,   expression(paste("Deleterious Allele Frequency/",hat(italic(q)))), cex=1.6, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.2,  expression(paste("Generations")), cex=1.6, adj=c(0.5, 0.5), xpd=NA)        
        # Legend
        legend(
               x       =  usr[2],
               y       =  usr[4]*0.9,
               legend = c(expression(italic(q)[italic(Y[I])]^italic(W)), 
                          expression(italic(q)[italic(Y[list(I,Auto)])]^italic(W))),
               lty     =  c(1,1),
               lwd     =  1.5,
               col     =  c(transparentColor(COLS[1], opacity=1),
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
    r.1.Dat  <-  makeDeterministicFigSimData(r = 1, x = x, h = h, generations = 8000)
    r.1.wbar.YI.t  <-  r.1.Dat$wbar.YI.t[round(r.1.Dat$YI.t,digits=4) > 0]
    autotime  <-  c(1:8000)
    q.wt.Auto  <-  p.t(u=10^-6, sdHom=0.01, h=h, t=autotime) 

## Panel B
# deleterious allele frequencies for r = 5
     # Make plot
     delMax  <-  r.1.Dat$XaOv.wt.t[1]*1.05
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.1.Dat)), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(qI.wt.t/max(qI.wt.t) ~ seq_along(qI.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.1.Dat)
        lines(q.wt.Auto/max(q.wt.Auto) ~ autotime, lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1))
        # axes
        axis(1, las=1, cex=1.2)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(paste(italic(h)==0.1,", ", italic(d)==1)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.2,  expression(paste("Generations")), cex=1.6, adj=c(0.5, 0.5), xpd=NA)        




    ###########
    # h = 0.01
    
    # make data for plotting
    h  <-  0.05
    r.10.Dat  <-  makeDeterministicFigSimData(r = 10, x = x, h = h, generations = 40000)
    r.10.wbar.YI.t  <-  r.10.Dat$wbar.YI.t[round(r.10.Dat$YI.t,digits=4) > 0]

    autotime  <-  c(1:40000)
    q.wt.Auto  <-  p.t(u=10^-6, sdHom=0.01, h=h, t=autotime) 

    qHat  <-  r.10.Dat$qY.wt.t[1]
    rBar  <-  round(nTot*x*qHat)

## Panel I
# deleterious allele frequencies for r = 5
     # Make plot
    delMax  <-  r.10.Dat$XaOv.wt.t[1]*1.05
     plot(NA, axes=FALSE, type='n', main='',xlim = c(0,nrow(r.10.Dat)), ylim = c(0,1), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        lines(qI.wt.t/max(qI.wt.t) ~ seq_along(qI.wt.t), lwd=1.5, lty=1, col=transparentColor(COLS[1], opacity=1), data=r.10.Dat)
        lines(q.wt.Auto/max(q.wt.Auto) ~ autotime, lwd=1.5, lty=1, col=transparentColor(COLS[2], opacity=1))
        # axes
        axis(1, las=1, cex=1.2)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel(0.04, 1.05, 'C', cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  1.05,   expression(paste(italic(h)==0.05,", ", italic(d)==10)), cex=1.5, adj=c(0.5, 0.5), xpd=NA)
        proportionalLabel( 0.5,  -0.2,  expression(paste("Generations")), cex=1.6, adj=c(0.5, 0.5), xpd=NA)        

}



###########################################
##  W-F Figs
PrFixFigStrongRecessive  <-  function() {

    # Import data for plotting
    dat  <-  read.csv(file = './data/PrFixFig_h0.01_s0.01_N1k_deterministic_q.csv', header=TRUE)
    dat$PrFixNe  <-  dat$PrFix*(dat$N/2)
    uniqueU  <-  unique(dat$Ufac)
    datN1k   <-  dat[dat$N == 1000,]
        datN1kUf1  <-  datN1k[datN1k$Ufac == uniqueU[1],]
        datN1kUf2  <-  datN1k[datN1k$Ufac == uniqueU[2],]
        datN1kUf3  <-  datN1k[datN1k$Ufac == uniqueU[3],]

    dat  <-  read.csv(file = './data/PrFixFig_h0.01_s0.01_N10k_deterministic_q.csv', header=TRUE)
    dat$PrFixNe  <-  dat$PrFix*(dat$N/2)
    uniqueU  <-  unique(dat$Ufac)
    datN10k   <-  dat[dat$N == 10000,]
        datN10kUf1  <-  datN10k[datN10k$Ufac == uniqueU[1],]
        datN10kUf2  <-  datN10k[datN10k$Ufac == uniqueU[2],]
        datN10kUf3  <-  datN10k[datN10k$Ufac == uniqueU[3],]


    # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  colfunc(3)

    # set plot layout
#    layout.mat <- matrix(c(1,2), nrow=1, ncol=2, byrow=TRUE)
#    layout     <- layout(layout.mat,respect=TRUE)

## Panel A
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
        # Legend



## Panel B
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(4,4,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c((10^-6), 2/(10^1)), log='y', ylab='', xlab='', cex.lab=1.2)
#     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c(0, 20), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80', log='y')
#        plotGrid(lineCol='grey80')
        box()
        # Points
        abline(h=(2/1000), lwd=1, lty=2, col=1)
#        abline(h=1, lwd=1, lty=1, col=1)
        points(PrFix ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[3], opacity=0.6), data=datN1kUf1)
        points(PrFix ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[3], opacity=0.6), data=datN1kUf2)
        points(PrFix ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[3], opacity=0.6), data=datN1kUf3)
        abline(h=(2/10000), lwd=1, lty=2, col=1)
        points(PrFix ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.6), data=datN10kUf1)
        points(PrFix ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.6), data=datN10kUf2)
        points(PrFix ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.6), data=datN10kUf3)
#        points(PrFix ~ x, pch=21, col=transparentColor('tomato', opacity=1), bg=transparentColor('tomato', opacity=0.6), data=IBMDat)
        # axes
        axis(1, las=1)
        axis(2, las=1, at=c(2/10^6, 2/10^5, 2/10^4, 2/10^3, 2/10^2, 2/10^1), 
                labels=c(expression(2/10^6), expression(2/10^5), expression(2/10^4), expression(2/10^3), expression(2/10^2), expression(2/10^1)))
        # Plot labels etc.
#        proportionalLabel(0.04, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.5, 1.1, 'Expand SLR', cex=1.25, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(-0.25, 0.5, expression(Fixation~Probability), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90, log='y')
        proportionalLabel( 0.5,  -0.225,  expression(paste("Inversion size (", italic(x),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        legend(
               x       =  0.25,
               y       =  2/(10^4),
               legend  =  c(
                            expression(paste(italic(N), " = ", 10^3)),
                            expression(paste(italic(N), " = ", 10^4))),
               pch     =  22,
               col     =  transparentColor(COLS[1], opacity=0),
               pt.bg   =  c(transparentColor(COLS[3], opacity=0.6),
                            transparentColor(COLS[1], opacity=0.6)),
               pt.cex  = 1.5,
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
        legend(
               x       =  0.3,
               y       =  3/(10^5),
               legend  =  c(
                            expression(paste(italic(U/s), " = ", 200)),
                            expression(paste(italic(U/s), " = ", 500)),
                            expression(paste(italic(U/s), " = ", 10^3))),
               pch     =  c(21,22,24),
               col     =  transparentColor(COLS[1], opacity=1),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
}



PrFixFigRecessive  <-  function() {

    # Import data for plotting
    dat  <-  read.csv(file = './data/PrFixFig_h0.1_s0.01_N1k_deterministic_q.csv', header=TRUE)
    dat$PrFixNe  <-  dat$PrFix*(dat$N/2)
    uniqueU  <-  unique(dat$Ufac)
    datN1k   <-  dat[dat$N == 1000,]
        datN1kUf1  <-  datN1k[datN1k$Ufac == uniqueU[1],]
        datN1kUf2  <-  datN1k[datN1k$Ufac == uniqueU[2],]
        datN1kUf3  <-  datN1k[datN1k$Ufac == uniqueU[3],]

    dat  <-  read.csv(file = './data/PrFixFig_h0.1_s0.01_N10k_deterministic_q.csv', header=TRUE)
    dat$PrFixNe  <-  dat$PrFix*(dat$N/2)
    uniqueU  <-  unique(dat$Ufac)
    datN10k   <-  dat[dat$N == 10000,]
        datN10kUf1  <-  datN10k[datN10k$Ufac == uniqueU[1],]
        datN10kUf2  <-  datN10k[datN10k$Ufac == uniqueU[2],]
        datN10kUf3  <-  datN10k[datN10k$Ufac == uniqueU[3],]


    # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  colfunc(3)



## Panel B
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(4,4,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c((10^-6), 2/(10^1)), log='y', ylab='', xlab='', cex.lab=1.2)
#     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,1), ylim = c(0, 20), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80', log='y')
#        plotGrid(lineCol='grey80')
        box()
        # Points
        abline(h=(2/1000), lwd=1, lty=2, col=1)
#        abline(h=1, lwd=1, lty=1, col=1)
        points(PrFix ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[3], opacity=0.6), data=datN1kUf1)
        points(PrFix ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[3], opacity=0.6), data=datN1kUf2)
        points(PrFix ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[3], opacity=0.6), data=datN1kUf3)
        abline(h=(2/10000), lwd=1, lty=2, col=1)
        points(PrFix ~ x, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.6), data=datN10kUf1)
        points(PrFix ~ x, pch=22, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.6), data=datN10kUf2)
        points(PrFix ~ x, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.6), data=datN10kUf3)
#        points(PrFix ~ x, pch=21, col=transparentColor('tomato', opacity=1), bg=transparentColor('tomato', opacity=0.6), data=IBMDat)
        # axes
        axis(1, las=1)
        axis(2, las=1, at=c(2/10^6, 2/10^5, 2/10^4, 2/10^3, 2/10^2, 2/10^1), 
                labels=c(expression(2/10^6), expression(2/10^5), expression(2/10^4), expression(2/10^3), expression(2/10^2), expression(2/10^1)))
        # Plot labels etc.
#        proportionalLabel(0.04, 1.05, 'B', cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(0.5, 1.1, 'Expand SLR', cex=1.25, adj=c(0.5, 0.5), xpd=NA, log='y')
        proportionalLabel(-0.25, 0.5, expression(Fixation~Probability), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90, log='y')
        proportionalLabel( 0.5,  -0.225,  expression(paste("Inversion size (", italic(x),")")), cex=1.2, adj=c(0.5, 0.5), xpd=NA, log='y')
        # Legend
        legend(
               x       =  0.925,
               y       =  2/10,
               legend  =  c(
                            expression(paste(italic(N), " = ", 10^3)),
                            expression(paste(italic(N), " = ", 10^4))),
               pch     =  22,
               col     =  transparentColor(COLS[1], opacity=0),
               pt.bg   =  c(transparentColor(COLS[3], opacity=0.6),
                            transparentColor(COLS[1], opacity=0.6)),
               pt.cex  = 1.5,
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
        legend(
               x       =  1,
               y       =  4/(10^2),
               legend  =  c(
                            expression(paste(italic(U/hs), " = ", 20)),
                            expression(paste(italic(U/hs), " = ", 50)),
                            expression(paste(italic(U/hs), " = ", 100))),
               pch     =  c(21,22,24),
               col     =  transparentColor(COLS[1], opacity=1),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
}


#########################################
## Supplement Fig.
tNeutral  <-  function(N, U, x, h, s) {
    log(N*U*x/2)/(h*s)
}
tSojourn  <-  function(N, s=0.01) {
    -2*log(2/(N*s))/(s)
}

approxAveSelectiveBenefit  <-  function(U, h, s, x){
    nTot  <-  10000
    Ufactor  <-  U/s
    u     <-  U/nTot
    qHat  <-  (U/(nTot*h*s))
    x     <-  x
    n     <-  nTot*x
    YI.0  <-  0.0001
    rVals  <-  c(0:(nTot*x*qHat))

    # Geometric Mean Variable.
    sBenGM  <-  c()
    tBen    <-  c()

    for(j in 1:length(rVals)) {
    
        # Draw random value for # del. mutations captured by inversion (r) given x
        r  <-  rVals[j]
            
        if(r == 0) {
            qHatDel  <-  0
        } else {qHatDel  <-  qHat}
    
    
        # Empty Frequency Vectors
        YI.t        <- c()
        XaOv.wt.t   <- c()
        XaSp.wt.t   <- c()
        qI.wt.t     <- c()
        qY.wt.t     <- c()
        XaOv.del.t  <- c()
        XaSp.del.t  <- c()
        qY.del.t    <- c()
        wbarYI.t    <- c()
    
        # Find equilibrium prior to inversion 
        eqs  <-  findEq(r=r, h=h, s=s, n=n, u=u, Ufactor=Ufactor, qHat=qHat, qHatDel=qHatDel)
        
        # First generation (all loci at equilibrium frequencies when inversion arises)
        YI.t[1]        <-  round(YI.prime(n=n, r=r, s=s, h=h, YI.t=YI.0, qt.I.wt=0, qt.Y.wt=eqs$qY.wt.init, qt.Y.del=eqs$qY.del.init, XaOv.t.wt=eqs$XaOv.wt.init, XaOv.t.del=eqs$XaOv.del.init), digits=9)
        XaOv.wt.t[1]   <-  round(XaOv.wt.prime(u=u, s=s, h=h, XaOv.t.wt=eqs$XaOv.wt.init, XaSp.t.wt=eqs$XaSp.wt.init), digits=9)
        XaSp.wt.t[1]   <-  round(XaSp.wt.prime(u=u, s=s, h=h, YI.t=YI.0, qt.I.wt=0, qt.Y.wt=eqs$qY.wt.init, XaOv.t.wt=eqs$XaOv.wt.init), digits=9)
        qI.wt.t[1]     <-  round(qI.wt.prime(n=n, r=r, u=u, s=s, h=h, YI.t=YI.0, qt.I.wt=0, XaOv.t.wt=eqs$XaOv.wt.init, XaOv.t.del=eqs$XaOv.del.init), digits=9)
        qY.wt.t[1]     <-  round(qY.wt.prime(u=u, s=s, h=h, qt.Y.wt=eqs$qY.wt.init, XaOv.t.wt=eqs$XaOv.wt.init), digits=9)
        XaOv.del.t[1]  <-  round(XaOv.del.prime(u=u, s=s, h=h, XaOv.t.del=eqs$XaOv.del.init, XaSp.t.del=eqs$XaSp.del.init), digits=9)
        XaSp.del.t[1]  <-  round(XaSp.del.prime(u=u, s=s, h=h, YI.t=YI.0, qt.Y.del=eqs$qY.del.init, XaOv.t.del=eqs$XaOv.del.init), digits=9)
        qY.del.t[1]    <-  round(qY.del.prime(u=u, s=s, h=h, qt.Y.del=eqs$qY.del.init, XaOv.t.del=eqs$XaOv.del.init), digits=9)
        wbarYI.t[1]    <-  round((YI.t[1]/YI.0), digits=9)
        if(wbarYI.t[1] < 1) {
            next
        }
    
        # Subsequent generations
        i=2
        while(wbarYI.t[i-1] > 1) {
            YI.t[i]        <-  round(YI.prime(n=n, r=r, s=s, h=h, YI.t=YI.t[i-1], qt.I.wt=qI.wt.t[i-1], qt.Y.wt=qY.wt.t[i-1], qt.Y.del=qY.del.t[i-1], XaOv.t.wt=XaOv.wt.t[i-1], XaOv.t.del=XaOv.del.t[i-1]), digits=9)
            XaOv.wt.t[i]   <-  round(XaOv.wt.prime(u=u, s=s, h=h, XaOv.t.wt=XaOv.wt.t[i-1], XaSp.t.wt=XaSp.wt.t[i-1]), digits=9)
            XaSp.wt.t[i]   <-  round(XaSp.wt.prime(u=u, s=s, h=h, YI.t=YI.t[i-1], qt.I.wt=qI.wt.t[i-1], qt.Y.wt=qY.wt.t[i-1], XaOv.t.wt=XaOv.wt.t[i-1]), digits=9)
            XaOv.del.t[i]  <-  round(XaOv.del.prime(u=u, s=s, h=h, XaOv.t.del=XaOv.del.t[i-1], XaSp.t.del=XaSp.del.t[i-1]), digits=9)
            XaSp.del.t[i]  <-  round(XaSp.del.prime(u=u, s=s, h=h, YI.t=YI.t[i-1], qt.Y.del=qY.del.t[i-1], XaOv.t.del=XaOv.del.t[i-1]), digits=9)
            qI.wt.t[i]     <-  round(qI.wt.prime(n=n, r=r, u=u, s=s, h=h, YI.t=YI.t[i-1], qt.I.wt=qI.wt.t[i-1], XaOv.t.wt=XaOv.wt.t[i-1], XaOv.t.del=XaOv.del.t[i-1]), digits=9)
            qY.wt.t[i]     <-  round(qY.wt.prime(u=u, s=s, h=h, qt.Y.wt=qY.wt.t[i-1], XaOv.t.wt=XaOv.wt.t[i-1]), digits=9)
            qY.del.t[i]    <-  round(qY.del.prime(u=u, s=s, h=h, qt.Y.del=qY.del.t[i-1], XaOv.t.del=XaOv.del.t[i-1]), digits=9)
            wbarYI.t[i]    <-  round((YI.t[i]/ YI.t[i-1]), digits=9)
            i  <-  i + 1
        }
        
        sBenGM[j]  <-  exp(mean(log(wbarYI.t))) - 1
        tBen[j]    <-  i
    }
   
   cbind(sBenGM,tBen)
}

suppTneutral_TsojournFig  <-  function(h = 0.1, s = 0.01, U = c(0.02, 0.05, 0.1), x=0.1) {

    # Generate Plotting Data
    N  <-  seq(1000,20000, by=50)

    # Colors
    # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  c(colfunc(5), "dodgerblue")

    # set plot layout
    layout.mat <- matrix(c(1,2,3), nrow=1, ncol=3, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)


## Panel A
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
    # Panel A: x = 0.1
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')    
        sBen  <-  simAveSelectiveBenefit(U=U[1], h=h, s=s, x=x)
        tSojourn0  <-  tSojourn(N=N, s=sBen[1,1])
        tSojourn1  <-  tSojourn(N=N, s=sBen[2,1])
        tNeut      <-  tNeutral(N=N, U=U[1], x=x, h=h, s=s)
     plot(NA, axes=FALSE, type='n', main='', xlim = c(min(N),max(N)), ylim = c(0,max(tSojourn1)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
#        abline(h=1, lwd=1.5, lty=1, col=1)
        lines(tSojourn0 ~ N, lwd=1.5, col=transparentColor(COLS[2], opacity=1))
        lines(tSojourn1 ~ N, lwd=1.5, col=transparentColor(COLS[3], opacity=1))
        lines(tNeut ~ N, lwd=1.5, lty=2, col=transparentColor(COLS[6], opacity=1))
        # benchmarks
        filter1  <-  tSojourn1 > tNeut
        abline(v=c(N[filter1][1]))
        # axes
        axis(1, las=1)
        axis(2, las=1)
        # Plot labels etc.
        proportionalLabel(-0.25,  0.5,   expression(Generations), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.25,  expression(paste(italic(N))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)


    # Panel B: x = 0.15
        sBen  <-  simAveSelectiveBenefit(U=U[2], h=h, s=s, x=x)
        tSojourn0  <-  tSojourn(N=N, s=sBen[1])
        tSojourn1  <-  tSojourn(N=N, s=sBen[2])
        tSojourn2  <-  tSojourn(N=N, s=sBen[3])
        tNeut  <-  tNeutral(N=N, U=U[2], x=x, h=h, s=s)
     plot(NA, axes=FALSE, type='n', main='', xlim = c(min(N),max(N)), ylim = c(0,max(tSojourn2)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        lines(tSojourn0 ~ N, lwd=1.5, col=transparentColor(COLS[2], opacity=1))
        lines(tSojourn1 ~ N, lwd=1.5, col=transparentColor(COLS[3], opacity=1))
        lines(tSojourn2 ~ N, lwd=1.5, col=transparentColor(COLS[4], opacity=1))
        lines(tNeut ~ N, lwd=1.5, lty=2, col=transparentColor(COLS[6], opacity=1))
        # benchmarks
        filter1  <-  tSojourn1 > tNeut
        filter2  <-  tSojourn2 > tNeut
        abline(v=c(N[filter1][1], N[filter2][1]))
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.5,  -0.25,  expression(paste(italic(N))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        abline(h=c(sBen[1,2], sBen[2,2], sBen[3,2]))

    # Panel C: x = 0.7
        sBen  <-  simAveSelectiveBenefit(U=U[3], h=h, s=s, x=x)
        tSojourn0  <-  tSojourn(N=N, s=sBen[1])
        tSojourn1  <-  tSojourn(N=N, s=sBen[2])
        tSojourn2  <-  tSojourn(N=N, s=sBen[3])
        tSojourn3  <-  tSojourn(N=N, s=sBen[4])
        tNeut  <-  tNeutral(N=N, U=U, x=x[3], h=h, s=s)
     plot(NA, axes=FALSE, type='n', main='', xlim = c(min(N),max(N)), ylim = c(0,max(tSojourn3)), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Lines
        lines(tSojourn0 ~ N, lwd=1.5, col=transparentColor(COLS[2], opacity=1))
        lines(tSojourn1 ~ N, lwd=1.5, col=transparentColor(COLS[3], opacity=1))
        lines(tSojourn2 ~ N, lwd=1.5, col=transparentColor(COLS[4], opacity=1))
        lines(tSojourn3 ~ N, lwd=1.5, col=transparentColor(COLS[5], opacity=1))
        lines(tNeut ~ N, lwd=1.5, lty=2, col=transparentColor(COLS[6], opacity=1))
        # benchmarks
        filter1  <-  tSojourn1 > tNeut
        filter2  <-  tSojourn2 > tNeut
        filter3  <-  tSojourn3 > tNeut
        abline(v=c(N[filter1][1], N[filter2][1], N[filter3][1]))
        # axes
        axis(1, las=1)
        axis(2, las=1, labels=NA)
        # Plot labels etc.
        proportionalLabel( 0.5,  -0.25,  expression(paste(italic(N))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)


        # Legend
        legend(
               x       =  usr[2]*0.965,
               y       =  usr[4],
               legend  =  c(
                            expression(paste(italic(N), " = ", 10^4)),
                            expression(paste(italic(N), " = ", 10^3)),
                            expression(paste(italic(N), " = ", 10^2))),
               pch     =  22,
               col     =  transparentColor(COLS[1], opacity=0),
               pt.bg   =  c(transparentColor(COLS[1], opacity=0.6),
                            transparentColor(COLS[2], opacity=0.6),
                            transparentColor(COLS[3], opacity=0.6)),
               pt.cex  = 1.5,
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
        legend(
               x       =  usr[2],
               y       =  usr[4]*0.775
               ,
               legend  =  c(
                            expression(paste(italic(U/s), " = ", 2)),
                            expression(paste(italic(U/s), " = ", 5)),
                            expression(paste(italic(U/s), " = ", 10))),
               pch     =  c(21,22,23),
               col     =  transparentColor(COLS[1], opacity=1),
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
}





sim_Tben_Tfix_Fig  <-  function() {

    # Import data for plotting
    dat1  <-  read.csv(file = './data/timeBen_Fix_h0.1_s0.01_U0.02_x0.2_N10k.csv', header=TRUE)
    dat2  <-  read.csv(file = './data/timeBen_Fix_h0.1_s0.01_U0.02_x0.5_N10k.csv', header=TRUE)
    dat3  <-  read.csv(file = './data/timeBen_Fix_h0.1_s0.01_U0.02_x0.8_N10k.csv', header=TRUE)

    # Colors
    colfunc  <-  colorRampPalette(c('#252525', 'grey70'))
    COLS     <-  colfunc(3)
#    COLS     <-  c("#252525", "grey70")

    # set plot layout
    layout.mat <- matrix(c(1), nrow=1, ncol=1, byrow=TRUE)
    layout     <- layout(layout.mat,respect=TRUE)

## Panel A
# Effect of initial del. mut. load on inversion frequency dynamics 
    # make plot
    par(omi=c(0.5, 0.5, 0.75, 0.5), mar = c(3,3,1,1), bty='o', xaxt='s', yaxt='s')    
     plot(NA, axes=FALSE, type='n', main='', xlim = c(0,10000), ylim = c(0,3000), ylab='', xlab='', cex.lab=1.2)
        usr  <-  par('usr')
        rect(usr[1], usr[3], usr[2], usr[4], col='white', border=NA)
        plotGrid(lineCol='grey80')
        box()
        # Points
        points(tBen ~ N, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.6), data=dat1)
        points(tFix ~ N, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[1], opacity=0.6), data=dat1)
        points(tBen ~ N, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=dat2)
        points(tFix ~ N, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[2], opacity=0.6), data=dat2)
        points(tBen ~ N, pch=21, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[3], opacity=0.6), data=dat3)
        points(tFix ~ N, pch=24, col=transparentColor(COLS[1], opacity=1), bg=transparentColor(COLS[3], opacity=0.6), data=dat3)
        # axes
        axis(1, las=1)
        axis(2, las=1)
#        axis(2, las=1, at=c((1/2)*10^-4, (1/2)*10^-3, (1/2)*10^-2, (1/2)*10^-1), 
#               labels=c(expression(1/2%*%10^-4), expression(1/2%*%10^-3), expression(1/2%*%10^-2), expression(2%*%10^-1)))
        # Plot labels etc.
        proportionalLabel(-0.25,  0.5,   expression(Generations), cex=1.2, adj=c(0.5, 0.5), xpd=NA, srt=90)        
        proportionalLabel( 0.5,  -0.25,  expression(paste(italic(N))), cex=1.2, adj=c(0.5, 0.5), xpd=NA)
        # Legend
        legend(
               x       =  usr[2]*0.16685,
               y       =  usr[4]*0.775,
               legend  =  c(
                            expression(paste(italic(bar(t)[fix]))),
                            expression(paste(italic(bar(t)[ben])))),
               pch     =  c(24,21),
               col     =  transparentColor(COLS[1], opacity=1),
               pt.bg   =  NA,
               pt.cex  = 1.2,
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )    
        legend(
               x       =  usr[2]*0.25,
               y       =  usr[4],
               legend  =  c(
                            expression(italic(x)~'='~0.8),
                            expression(italic(x)~'='~0.5),
                            expression(italic(x)~'='~0.2)),
               pch     =  22,
               col     =  NA,
               pt.bg   =  c(transparentColor(COLS[1], opacity=0.6),
                            transparentColor(COLS[2], opacity=0.6),
                            transparentColor(COLS[3], opacity=0.6)),
               pt.cex  = 1.5,
               cex     =  1,
               xjust   =  1,
               yjust   =  1,
               bty     =  'n',
               border  =  NA
               )
}