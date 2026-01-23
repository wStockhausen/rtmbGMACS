#'
#'@title Calculate probability-at-size (pre-molt) for immature crab that a molt is to maturity, by sex and shell condition
#'
#'@description Function to calculate probability-at-size (pre-molt) for immature crab that a molt is to maturity, by sex and shell condition.
#'
#'@param mc - model configuration object
#'@param showPlot - flag to show plots
#'
#'@return list with elements:
#'prM2M_cxz
#'prM2M_yxsz
#'
#'@details None.
#'
#'@import reshape2
#'@import ggplot2
#'
#'@export
#'
calcPrM2M<-function(mc,showPlot=TRUE){
    if (mc$type!='TC'){
        throwModelTypeError(mc$type,'calcPrMoltToMaturity()');
    }
    
    d<-mc$dims;
    p<-mc$params$prM2M;
    
    prM2M_cxz  <- dimArray(mc,'pc_prM2M.x.z',val=0);    
    prM2M_yxsz <- dimArray(mc,'y.x.s.z',val=0);    
    for (t in names(p$blocks)){
        tb<-p$blocks[[t]];
        yrs<-as.character(tb$years);
        for (x in d$x$nms){
            z50 <- tb$z50_x[x];
            sdv <- tb$sdv_x[x];
            prM2M_cxz[t,x,]<-plogis(d$z$vls,z50,sdv);
            for (m in d$m$nms){
                for (s in d$s$nms) {
                    for (y in yrs) prM2M_yxsz[y,x,s,]<-prM2M_cxz[t,x,];
                }#s
            }#m      
        }#x
    }#t
    
    if (showPlot){
        mdfr<-melt(prM2M_cxz,value.name='val')
        pl <- ggplot(aes(x=z,y=val,color=x),data=mdfr)
        pl <- pl + geom_line()
        pl <- pl + labs(x='size (mm)',y='pr(molt-to-maturity|size, molt)')
        pl <- pl + guides(color=guide_legend('sex'));
        pl <- pl + facet_grid(pc~.);
        print(pl);
    }
    
    return(list(prM2M_cxz=prM2M_cxz,prM2M_yxsz=prM2M_yxsz));
}