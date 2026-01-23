#'
#'@title Calculate selectivity/retention functions.
#'
#'@description Function to calculate selectivity/retention functions.
#'
#'@param mc - model configuration object
#'@param showPlot - flag to show plots
#'
#'@return array sel_cz with dimensions 'pc', 'z'. The label associated with
#'each pc can be obtained as mc$dims$selfcns$lbls[pc].
#'
#'@import reshape2
#'@import ggplot2
#'
#'@details None.
#'
#'@export
#'
calcSelFcns<-function(mc,showPlot=TRUE){
    d<-mc$dims;
    ps<-mc$params$selfcns;
    
    sel_cz<-dimArray(mc,'pc_selfcns.z',val=NA);
    for (c in 1:d$selfcns$n){
        si<-ps[[c]];
        sel_cz[c,]<-calcSelectivity(si$type,d$z$vls,si$params);
    }
    if (showPlot){
        mdfr<-melt(sel_cz,value.name='val');
        mdfr$selfcns<-d$selfcns$lbls[mdfr$pc];
        p <- ggplot(aes(x=z,y=val,color=selfcns,shape=selfcns),data=mdfr);
        p <- p + geom_point(size=5,alpha=0.5);
        p <- p + geom_line();
        p <- p + labs(x='size (mm)',y='selectivity/retention')
        p <- p + guides(color=guide_legend('',order=1),
                        shape=guide_legend('',order=1));
        print(p)
    }
    return(sel_cz);
}