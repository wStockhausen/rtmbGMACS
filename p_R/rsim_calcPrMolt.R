#'
#'@title Calculate probability at size of molting for immature crab by sex, shell condition
#'
#'@description Function to calculate probability at size of molting for immature crab by sex, shell condition.
#'
#'@param mc - model configuration object
#'@param showPlot - flag to show plots
#'
#'@return prMolt_yxsz
#'
#'@details None.
#'
#'@import reshape2
#'@import ggplot2
#'
#'@export
#'
calcPrMolt<-function(mc,showPlot=TRUE){
    if (mc$type!='TC'){
        throwModelTypeError(mc$type,'TC','calcPrMolt()');
    }
    
    d<-mc$dims;
    p<-mc$params$molting;
    
    prMolt_yxsz <- dimArray(mc,'y.x.s.z',val=0);    
    mdfr<-NULL;
    for (t in names(p$blocks)){
        tb<-p$blocks[[t]];
        yrs<-as.character(tb$years);
        for (x in d$x$nms){
            for (s in d$s$nms) {
                z50 <- tb$z50_xs[x,s];
                sdv <- tb$sdv_xs[x,s];
                mp_z<-dimArray(mc,'z');
                mp_z[] <- 1.0 - plogis(d$z$vls,z50,sdv);
                mdfrp<-melt(mp_z,value.name='val');
                mdfrp$x<-x;
                mdfrp$s<-s;
                mdfrp$fac<-paste(x,s,sep=', ');
                mdfrp$t<-t;
                mdfr<-rbind(mdfr,mdfrp);
                for (y in yrs) prMolt_yxsz[y,x,s,]<-mp_z;
            }#s
        }#x
    }#t
    
    if (showPlot){
        pl <- ggplot(aes(x=z,y=val,color=x),data=mdfr)
        pl <- pl + geom_line()
        pl <- pl + labs(x='size (mm)',y='pr(molt|size)')
        pl <- pl + guides(color=guide_legend('sex'));
        pl <- pl + facet_grid(t~s);
        print(pl);
    }
    
    return(prMolt_yxsz);
}