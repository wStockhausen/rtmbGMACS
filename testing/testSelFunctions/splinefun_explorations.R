#--data values
x = rep(c(1:10),3);
y = sin(2*pi*x/10);           #--true y's
g = rep.int(1:3,c(10,10,10));
yo = y+rnorm(length(x),0,0.1);#--observed y's
dfro = tibble::tibble(x=x,true=y,obs=yo,grp=g) |>
         tidyr::pivot_longer(c(true,obs));
ggplot(dfro,aes(x=x,y=value,colour=name,group=grp)) +
  geom_point(data=dfro |> dplyr::filter(name=="obs")) +
  geom_line(data=dfro |> dplyr::filter(name=="true"));


#--parameter values
xs = c(1,4,7,10);#--first and last values are dumm
ys = sin(2*pi*xs/10);
dx = 0.1;
xs = c(xs[1]-dx*rev(1:3),xs,xs[length(xs)]+dx*(1:3));
ys = c(rep(ys[1],3),ys,rep(ys[length(ys)],3));

dfrs  = tibble::tibble(x=xs,value=ys,name="knots")
ggplot(dfrs,aes(x=x,y=value,colour=name)) +
  geom_point() +
  geom_line();

#--add parameter values to function
spline = splinefun(xs,ys,method="natural");

#--calculate spline at values that include the x's data
xp = seq(from=-3,to=14,by=0.1);
dfrsv = tibble::tibble(x=xp,y=spline(xp,0),d1ydx1=spline(xp,1),d2ydx2=spline(xp,2),d3ydx3=spline(xp,3)) |>
        tidyr::pivot_longer(!x) |>
        dplyr::mutate(name=factor(name,levels=c("y","d1ydx1","d2ydx2","d3ydx3")));


ggplot(dfrsv,aes(x=x,y=value,colour=name)) +
  geom_line() + geom_point(data=dfrs) +
  geom_point(data=dfro)
