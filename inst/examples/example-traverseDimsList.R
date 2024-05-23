#--define nested dimensions list
#----nested dimensions represent sex ("x"), maturity ("m"),
#----shell condition ("s"), age ("a"), and size ("z")
vZs = list( MALE=  list(IMMATURE=list(`NEW SHELL`=list(`NO AGEs`=seq(40,60,5))),
                          MATURE=list(`NEW SHELL`=list(`PTM 0`=seq(80,120,10)),
                                      `OLD SHELL`=list(`PTM 1`=seq(80,120,10),
                                                       `PTM 2`=seq(80,120,10)))),
          FEMALE=list(IMMATURE=list(`NEW SHELL`=list(`NO AGEs`=seq(25,50,5))),
                        MATURE=list(`NEW SHELL`=list(`PTM 0`=seq(50,100,10)),
                                    `OLD SHELL`=list(`PTM 1`=seq(50,100,10),
                                                     `PTM 2`=seq(50,100,10)))));
attr(vZs,"dmnms")<-c("x","m","s","a","z");#--define dmnms attribute
require(rtmbGMACS);
dfr = traverseDimsList(vZs,"z");
