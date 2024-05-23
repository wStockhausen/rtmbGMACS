plotResults<-function(dfr,showPlots=FALSE){
  hasWTSPLOTS = requireNamespace("wtsPlots");
  #--plot results
  ps = list();
  for (v in unique(dfr$name)){
    #--v=unique(dfr$name)[1];
    cat("plotting",v,"\n")
    dfrp = dfr |> dplyr::filter(name==v);
    nD = dfrp$nDim[1];#--number of dimensions
    cat("\tnD =",nD,"\n")
    if (nD==1){
      if (nrow(dfrp)==length(unique(dfrp$model))){
        #--plot scalar-valued quantities
        p = ggplot(dfrp,aes(x=name,y=Est,ymin=Est-Std,ymax=Est+Std,color=model,fill=model)) +
              #geom_boxplot(aes(lower=Est-Std,middle=Est,upper=Est+Std),stat="identity") +
              geom_point(size=6) +
              geom_errorbar() +
              labs(x=v);
        if (hasWTSPLOTS) p = p + wtsPlots::getStdTheme();
        if (showPlots) print(p);
        ps[[v]] = p;
      } else {
        #--plot vectors-valued quantities
        p = ggplot(dfrp,aes(x=i1,y=Est,ymin=Est-Std,ymax=Est+Std,color=model,fill=model)) +
              geom_ribbon(alpha=0.3,colour=NA) + geom_line() +
              labs(x="index",y=v);
        if (hasWTSPLOTS) p = p + wtsPlots::getStdTheme();
        if (showPlots) print(p);
        ps[[v]] = p;
      }
    } else if (nD==2) {
      #--plot matrices
      p = ggplot(dfrp,aes(x=i1,y=Est,ymin=Est-Std,ymax=Est+Std,color=model,fill=model)) +
            geom_ribbon(alpha=0.3,colour=NA) + geom_line() +
            facet_wrap(facets=paste0("i",2:nD)) +
            labs(x="index",y=v);
      if (hasWTSPLOTS) p = p + wtsPlots::getStdTheme();
      if (showPlots) print(p);
      ps[[v]] = p;
    }
  }
  return(ps);
}
