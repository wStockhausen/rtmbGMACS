#'
#'@title Calculate size transition matrices by year, sex, shell condition
#'
#'@description Function to calculate size transition matrices by year, sex, shell condition.
#'
#'@param mc - model configuration object
#'
#'@return list with elements:
#'mnZAM_cxz: mean size after molt, by time block and sex 
#'T_cxzz: 4d array with size transition matrices by time block/sex
#'T_yxszz: 5d array with size transition matrices by year/sex/shell condition
#'
#'@import ggplot2
#'@import reshape2
#'
#'@details None.
#'
#'@export
#'
calcZTM<-function(mc,showPlot=TRUE){
    d<-mc$dims;
    p<-mc$params$growth;
    
    mnZAM_cxz   <- dimArray(mc,'pc_growth.x.z')
    prZAM_cxzz  <- dimArray(mc,'pc_growth.x.z.zp')
    mnZAM_yxsz  <- dimArray(mc,'y.x.s.z')
    prZAM_yxszz <- dimArray(mc,'y.x.s.z.zp')
    mdfr.pr<-NULL;
    mdfr.mn<-NULL;
    for (t in names(p$blocks)){
        tb<-p$blocks[[t]];
        yrs<-as.character(tb$years);
        prZAM_xzz   <- dimArray(mc,'x.z.zp');#size transition matrix
        for (x in d$x$nms){
            mnZAM_cxz[t,x,]<-exp(tb$a_x[x])*(d$z$vls^tb$b_x[x]);
            #mnZs = grA*mfexp(log(grB/grA)/log(zGrB/zGrA)*log(zBs/zGrA));
            grA<-tb$a_x[x]; zGrA<-tb$za_x[x];
            grB<-tb$b_x[x]; zGrB<-tb$zb_x[x];
            mnZAM_cxz[t,x,]<-grA*exp(log(grB/grA)/log(zGrB/zGrA)*log(d$z$vls/zGrA));
            for (z in 1:(d$z$n-1)){ #looping over pre-molt size bins
                idx<-z:d$z$n;#index over CUTPOINTS up to start of final bin
                cumZAM<-pgamma(d$zc$vls[idx],shape=mnZAM_cxz[t,x,z]/tb$s_x[x],scale=tb$s_x[x]);#integrated up each cutpoint
                prZAM<-c(first_difference(cumZAM),1.0-cumZAM[length(cumZAM)]);#contribution to each bin
                #TODO: no truncation here!!
                prZAM<-prZAM/sum(prZAM);
                prZAM_xzz[x,z,z:d$z$n]<-prZAM;#note that ROWS (z) here are pre-molt, COLUMNS (zp) are post-molt
            }#z
            prZAM_xzz[x,d$z$n,d$z$n]<-1.0;
            prZAM_cxzz[t,x,,]<-t(prZAM_xzz[x,,]);#columns (zp) now pre-molt, rows (z) post-molt
            for (s in d$s$nms){
                for (y in yrs) {
                    #indep of shell condition
                    mnZAM_yxsz[y,x,s,] <-mnZAM_cxz[t,x,];
                    prZAM_yxszz[y,x,s,,]<-prZAM_cxzz[t,x,,];
                }#y
            }#s
        }#x
    }#t
    if (showPlot){
        mdfr<-melt(mnZAM_cxz,value.name='val');
        pl <- ggplot(aes(x=z,y=val,color=x,shape=x),data=mdfr);
        pl <- pl + geom_abline(intercept=0,slope=1,linetype=3,color='black')
        pl <- pl + geom_line(size=2);
        pl <- pl + geom_point(size=5);
        pl <- pl + guides(color=guide_legend(''),shape=guide_legend(''));
        pl <- pl + labs(x='pre-molt size (mm)',y='mean post-molt size (mm)',title='Mean Growth');
        pl <- pl + facet_grid(pc~.);
        print(pl);
        mdfr<-melt(prZAM_cxzz,value.name='val');
        pl <- ggplot(aes(x=zp,y=z,fill=val,size=val),data=mdfr);
        pl <- pl + geom_point(alpha=0.6,shape=21);
        pl <- pl + scale_size_area(max_size=10);
        pl <- pl + scale_fill_gradient()
        pl <- pl + geom_abline(intercept=0,slope=1,linetype=3,color='black')
        pl <- pl + labs(x='pre-molt size (mm)',y='post-molt size (mm)',title='Growth Transition Matrices');
        pl <- pl + guides(fill=guide_colorbar(expression(pr*bgroup("(",paste(Z[post],"|",Z[pre]),")")),order=1,alpha=1),
                        size=guide_legend('',order=2));
        pl <- pl + facet_grid(pc~x);
        print(pl);
    }
    return(list(mnZAM_cxz=mnZAM_cxz,T_cxzz=prZAM_cxzz,mnZAM_yxsz=mnZAM_yxsz,T_yxszz=prZAM_yxszz))
}