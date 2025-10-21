library(dplyr)
library(ggplot2)
library(lhs)
library(dplyr) #need for bind_rows for incidences_sum
library(reshape2) #need for melt for total.inc_long

# runSEIR takes the SEIR model parameters and initial condition 
# and returns a time series for each group
SEIIRfunct <- function(beta, omega, sigma, v.o, v.a, v.m, v.sh, v.abx,
                       epsilon.a, epsilon.a.y,
                       epsilon.s, epsilon.s.y, epsilon.m.T, epsilon.m.T.y, 
                       epsilon.s.T, epsilon.s.T.y,
                       alpha.m, alpha.s,
                       tau, q.o, q, delta, delta.y,
                       theta.m, theta.s,
                       gamma.a, gamma.m, gamma.s,
                       gamma.m.abx, gamma.s.abx,
                       mu.m, mu.m.y, mu.s, mu.s.y,
                       initial.state,
                       step.size = 1,
                       freq.dependent=TRUE,
                       final.only=FALSE){
  
  #If the model is frequency dependent we modify beta based on the total population size
  beta.divisor <- ifelse(freq.dependent,
                         initial.state["Sy"]+initial.state["S"]+
                           initial.state["Ey"]+initial.state["E"]+
                           initial.state["Iy_a_sh"]+ initial.state["Ia_sh"]+ 
                           initial.state["Iy_m_syU"]+initial.state["Im_syU"]+
                           initial.state["Iy_m_sh"]+initial.state["Im_sh"]+
                           initial.state["Iy_m_syT"]+initial.state["Im_syT"]+
                           initial.state["Iy_m_abx"]+initial.state["Im_abx"]+
                           initial.state["Iy_s_syU"]+initial.state["Is_syU"]+
                           initial.state["Iy_s_sh"]+initial.state["Is_sh"]+
                           initial.state["Iy_s_syT"]+initial.state["Is_syT"]+
                           initial.state["Iy_s_abx"]+initial.state["Is_abx"]+
                           initial.state["R"]+initial.state["R_abx"]+initial.state["D"],
                         1)
  ## Define the number of columns to go in the final dataset 
  col.count <-37 ## change this every time you add a new arrow
  
  #create the parameter vector.
  param <- c(beta.scaled=beta/beta.divisor, omega=omega, sigma=sigma, 
             v.o=v.o, v.a=v.a, v.m=v.m, v.sh=v.sh, v.abx=v.abx,
             epsilon.a=epsilon.a, epsilon.a.y=epsilon.a.y, epsilon.s=epsilon.s, epsilon.s.y=epsilon.s.y,
             epsilon.m.T=epsilon.m.T, epsilon.m.T.y=epsilon.m.T.y, epsilon.s.T=epsilon.s.T, epsilon.s.T.y=epsilon.s.T.y,
             alpha.m=alpha.m, alpha.s=alpha.s,
             tau=tau,q.o=q.o, q=q, delta=delta, delta.y=delta.y,
             theta.m=theta.m, theta.s=theta.s,
             gamma.a=gamma.a, gamma.m=gamma.m, gamma.s=gamma.s,
             gamma.m.abx=gamma.m.abx, gamma.s.abx=gamma.s.abx,
             mu.m=mu.m, mu.m.y=mu.m.y, mu.s=mu.s, mu.s.y=mu.s.y)
  
  #Since we are not using a fancy solver we will need to run this on our own.
  #note that the simulation ends once there are 0 people in groups E or I
  t <- 0
  y <- initial.state
  
  SEIIR.output <- matrix(ncol=length(initial.state)+1, nrow=1)
  colnames(SEIIR.output) <- c("time", "Sy", "S", "Ey","E",
                              "Iy_a_sh","Ia_sh",
                              "Iy_m_syU", "Im_syU", "Iy_m_sh", "Im_sh", 
                              "Iy_m_syT", "Im_syT", "Iy_m_abx", "Im_abx",
                              "Iy_s_syU", "Is_syU", "Iy_s_sh", "Is_sh", 
                              "Iy_s_syT", "Is_syT", "Iy_s_abx", "Is_abx",
                              "R", "D", "R_abx")
  SEIIR.output[1,] <- c(t,y)
  SEIIR.output
  incidences <- matrix(ncol=col.count+1, nrow=1)
  colnames(incidences) <- c("time",
                            "incident.aging",
                            "incident.exposed.y",
                            "incident.exposed",
                            "incident.cases.a.y", 
                            "incident.cases.a",
                            "incident.cases.symp.mU.y", "incident.cases.symp.mU",
                            "incident.cases.symp.mT.y", "incident.cases.symp.mT",
                            "incident.cases.symp.sU.y", "incident.cases.symp.sU",
                            "incident.cases.symp.sT.y", "incident.cases.symp.sT",
                            "shedding.cases.mU.y", "shedding.cases.mU", "shedding.cases.mT.y", "shedding.cases.mT",
                            "abx.cases.m.y", "abx.cases.m",
                            "shedding.cases.s.y", "shedding.cases.s", "abx.cases.s.y", "abx.cases.s",
                            "nat.recovered.a.y", "nat.recovered.a",
                            "nat.recovered.m.y", "nat.recovered.m", 
                            "died.m.y", "died.m", "abx.recovered.m.y", "abx.recovered.m",
                            "nat.recovered.s.y", "nat.recovered.s", 
                            "died.s.y", "died.s", "abx.recovered.s.y", "abx.recovered.s")
  incidences[1,] <- rep(0,col.count+1)
  
  while(
    y["Ey"]>0 |
    y["E"]>0 |
    y["Iy_a_sh"]>0 |
    y["Ia_sh"]>0 |
    y["Iy_m_syU"]>0 | y["Im_syU"]>0 | y["Iy_m_sh"]>0 | y["Im_sh"]>0 | 
    y["Iy_m_syT"]>0 | y["Im_syT"]>0 | y["Iy_m_abx"]>0 | y["Im_abx"]>0 |
    y["Iy_s_syU"]>0 | y["Is_syU"]>0 | y["Iy_s_sh"]>0 | y["Is_sh"]>0 | 
    y["Iy_s_syT"]>0 | y["Is_syT"]>0 | y["Iy_s_abx"]>0 | y["Is_abx"]>0) {
    
    t <- t+step.size 
    print(y)
    if(y["R"] > 500000){
      stop("ERROR: generating people from thin air")
    }
    
    lambda <- (param["beta.scaled"]*(y["Iy_s_syU"] + y["Iy_s_syT"])) + (v.o*param["beta.scaled"]*(y["Is_syU"] + y["Is_syT"])) +
      (v.sh*param["beta.scaled"]*y["Iy_s_sh"]) + (v.o*v.sh*param["beta.scaled"]*y["Is_sh"]) +
      (v.a*param["beta.scaled"]*y["Iy_a_sh"]) + (v.o*v.a*param["beta.scaled"]*y["Ia_sh"]) +
      ((v.m*param["beta.scaled"])*(y["Iy_m_syU"]+y["Iy_m_syT"])) + ((v.o*v.m*param["beta.scaled"])*(y["Im_syU"]+y["Im_syT"])) + 
      (v.sh*v.m*param["beta.scaled"]*y["Iy_m_sh"]) + (v.o*v.sh*v.m*param["beta.scaled"]*y["Im_sh"]) + 
      (v.m*v.abx*param["beta.scaled"]*y["Iy_m_abx"]) + (v.o*v.m*v.abx*param["beta.scaled"]*y["Im_abx"])
    lambda
    #at this point, lambda is basically beta, looks complicated since is so many different infected types
    #is a number, should look like a number. print to look at
    
    #calculate the probability of infection and recovery in this time step
    #probability that thing happens in a given time step, updates every time step. range 0-1
    #need one of these per arrow, moderate/severe specific etc
    Pr.expose <- 1-exp(-step.size*lambda)
    
    Pr.aging <- 1-exp(-step.size*param["omega"])
    
    Pr.infect.all <- 1-exp(-step.size*param["sigma"])

    Pr.mU.shed <- 1-exp(-step.size*param["alpha.m"]) #the rate at which moderates who don't seek care move to shed
    Pr.mT.shed <- 1-exp(-step.size*param["alpha.m"]) #this is only rate of move to shed
    Pr.mT.abx <- 1-exp(-step.size*param["delta"]*param["tau"]) #this is only rate of moderates who seek care move to abx
    # Pr.mT.shed <- 1-exp(-step.size*(1-param["q"])*param["alpha.m"])
    # Pr.mT.abx <- 1-exp(-step.size*param["q"]*param["delta"]*param["tau"])
    Pr.sU.shed <- 1-exp(-step.size*param["alpha.s"]) #the rate at which severes who don't seek care move to shed
    Pr.sT.abx <- 1-exp(-step.size*param["tau"]) #the rate at which severes who do seek care move to abx
    
    Pr.a.recover <- 1-exp(-step.size*param["gamma.a"]) #rate at which asymptomatics recover
    Pr.m.shed.recover <- 1-exp(-step.size*param["gamma.m"]) #this only incorporates the rate at which shedding moderates recover
    Pr.m.shed.die <- 1-exp(-step.size*param["mu.m"]) #this only incorporates the rate at which shedding moderates die
    # Pr.m.shed.recover <- 1-exp(-step.size*(1-param["theta.m"])*param["gamma.m"])
    # Pr.m.shed.die <- 1-exp(-step.size*param["theta.m"]*param["mu.m"])
    Pr.m.abx.recover <- 1-exp(-step.size*param["gamma.m.abx"]) #rate at which moderates who receive abx recover
    Pr.s.shed.recover <- 1-exp(-step.size*param["gamma.s"]) #this only incorporates the rate of recover among shedding severes
    Pr.s.shed.die <- 1-exp(-step.size*param["mu.s"]) #this only incorporates the rate of die among shedding severes
    # Pr.s.shed.recover <- 1-exp(-step.size*(1-param["theta.s"])*param["gamma.s"])
    # Pr.s.shed.die <- 1-exp(-step.size*param["theta.s"]*param["mu.s"])
    Pr.s.abx.recover <- 1-exp(-step.size*param["gamma.s.abx"]) #rate of recover among severes who receive abx
    
    ## Hailey's additions
    Pr.expose.y <- 1-exp(-step.size*lambda)
    Pr.infect.all.y <- 1-exp(-step.size*param["sigma"])
    Pr.mU.shed.y <- 1-exp(-step.size*param["alpha.m"]) #the rate at which moderates who don't seek care move to shed
    Pr.mT.shed.y <- 1-exp(-step.size*param["alpha.m"]) #this is only rate of move to shed
    Pr.mT.abx.y <- 1-exp(-step.size*param["delta.y"]*param["tau"]) #this is only rate of moderates who seek care move to abx
    Pr.sU.shed.y <- 1-exp(-step.size*param["alpha.s"]) #the rate at which severes who don't seek care move to shed
    Pr.sT.abx.y <- 1-exp(-step.size*param["tau"]) #the rate at which severes who do seek care move to abx
    
    Pr.a.recover.y <- 1-exp(-step.size*param["gamma.a"]) #rate at which children asymptomatics recover
    Pr.m.shed.recover.y <- 1-exp(-step.size*param["gamma.m"]) #this only incorporates the rate at which shedding moderates recover
    Pr.m.shed.die.y <- 1-exp(-step.size*param["mu.m.y"]) #this only incorporates the rate at which shedding moderates die
    Pr.m.abx.recover.y <- 1-exp(-step.size*param["gamma.m.abx"]) #rate at which moderates who receive abx recover
    Pr.s.shed.recover.y <- 1-exp(-step.size*param["gamma.s"]) #this only incorporates the rate of recover among shedding severes
    Pr.s.shed.die.y <- 1-exp(-step.size*param["mu.s.y"]) #this only incorporates the rate of die among shedding severes
    Pr.s.abx.recover.y <- 1-exp(-step.size*param["gamma.s.abx"]) #rate of recover among severes who receive abx
    
    #draw random variable from binomial distribution for new number of
    #using prob from above, is the actual number based on prob and denom (num ppl in that box)
    #should be whole numbers. if aren't something's wrong
    #rbinom(num observations, num trials, prob success each trial)
    incident.exposed <- rbinom(1, y["S"], Pr.expose) #the number who go So to Eo
    
    incident.cases.all <- rbinom(1, y["E"], Pr.infect.all) #this is the number that leave E
    incident.cases.a <- round(incident.cases.all*param["epsilon.a"],0) #this the number that enter I_a
    #now only have (number that leave)-(number that go to A) left to go to symptomatics
    incident.cases.symp <- (incident.cases.all-incident.cases.a) #this is the number that go to symptomatic I's
    incident.cases.symp.s <- round(incident.cases.symp*param["epsilon.s"],0) #this is the number that go to severe symptomatic I's
    incident.cases.symp.m <- (incident.cases.symp-incident.cases.symp.s) #this is the number that goes to moderate symptomatic I's
    incident.cases.symp.mT <- round(incident.cases.symp.m*param["epsilon.m.T"],0) #this is the number that go to moderate seek care
    incident.cases.symp.mU <- (incident.cases.symp.m-incident.cases.symp.mT) #this is the number that go to moderate don't seek care
    incident.cases.symp.sT <- round(incident.cases.symp.s*param["epsilon.s.T"],0) #this is the number that go to severe seek care
    incident.cases.symp.sU <- (incident.cases.symp.s-incident.cases.symp.sT) #this is the number that go to severe don't seek care
    # incident.cases.all
    # sum(incident.cases.a,incident.cases.symp.mU,incident.cases.symp.mT,incident.cases.symp.sU,incident.cases.symp.sT)
    
    shedding.cases.mU <- rbinom(1, y["Im_syU"], Pr.mU.shed) #the number who move to moderate shedding
    m.who.will.abx <- round(y["Im_syT"]*param["q"],0) #the number of moderates who seek care who will eventually receive abx#are treating q like a proportion, though this not actually correct. the math works though and conceptually is ok via prop.m.abx equation
    shedding.cases.mT <- rbinom(1, y["Im_syT"]-m.who.will.abx, Pr.mT.shed) #the number of moderates who seek care who don't receive abx and move to shed
    abx.cases.m <- rbinom(1, m.who.will.abx, Pr.mT.abx) #the number of moderates who receive abx at this step
    # shedding.cases.mT <- rbinom(1, y["Im_syT"], Pr.mT.shed)
    # abx.cases.m <- rbinom(1, y["Im_syT"], Pr.mT.abx)
    shedding.cases.s <- rbinom(1, y["Is_syU"], Pr.sU.shed) #the number of severes who don't seek care and move to shed
    abx.cases.s <- rbinom(1, y["Is_syT"], Pr.sT.abx) #the number of severes who seek care and receive abx
    
    nat.recovered.a <- rbinom(1, y["Ia_sh"], Pr.a.recover) #the number of asymptomatics who recover
    m.who.will.die <- round(y["Im_sh"]*param["theta.m"],0) #this is the whole number of moderates who will eventually die
    nat.recovered.m <- rbinom(1, y["Im_sh"]-m.who.will.die, Pr.m.shed.recover) #this is the actual number of I_sh who naturally recover. this incorporates the stochastic prob that will actually move in this step #rbinom requires whole numbers
    died.m <- rbinom(1, m.who.will.die, Pr.m.shed.die) #the number of moderates who move from I_sh to dead
    # nat.recovered.m <- rbinom(1, y["Im_sh"], Pr.m.shed.recover)
    # died.m <- rbinom(1, y["Im_sh"], Pr.m.shed.die)
    abx.recovered.m <- rbinom(1, y["Im_abx"], Pr.m.abx.recover) #number of moderates who move from abx to recover in this step
    s.who.will.die <- round(y["Is_sh"]*param["theta.s"],0) #this is the whole number of severes who will die
    nat.recovered.s <- rbinom(1, y["Is_sh"]-s.who.will.die, Pr.s.shed.recover) #number of severes who move from Is_sh to recover w/o abx
    died.s <- rbinom(1, s.who.will.die, Pr.s.shed.die) #number of Is_sh who die in this time step
    abx.recovered.s <- rbinom(1, y["Is_abx"], Pr.s.abx.recover) #number of severes who recover after abx in this time step
    
    ## Hailey's additions
    incident.aging <- rbinom(1, y["Sy"], Pr.aging) #the number who go Sy to So
    incident.exposed.y <- rbinom(1, y["Sy"], Pr.expose.y) #the number who go Sy to Ey
    incident.cases.all.y <- rbinom(1, y["Ey"], Pr.infect.all.y) #this is the number that leave Ey
    incident.cases.a.y <- round(incident.cases.all.y*param["epsilon.a.y"],0) #this the number that enter Iy_a
    incident.cases.symp.y <- (incident.cases.all.y-incident.cases.a.y) #this is the number that go to symptomatic Iy's
    incident.cases.symp.s.y <- round(incident.cases.symp.y*param["epsilon.s.y"],0)
    incident.cases.symp.m.y <- (incident.cases.symp.y-incident.cases.symp.s.y) #this is the number that goes to moderate symptomatic I's
    incident.cases.symp.mT.y <- round(incident.cases.symp.m.y*param["epsilon.m.T.y"],0) #this is the number that go to moderate seek care
    incident.cases.symp.mU.y <- (incident.cases.symp.m.y-incident.cases.symp.mT.y) #this is the number that go to moderate don't seek care
    incident.cases.symp.sT.y <- round(incident.cases.symp.s.y*param["epsilon.s.T.y"],0) #this is the number that go to severe seek care
    incident.cases.symp.sU.y <- (incident.cases.symp.s.y-incident.cases.symp.sT.y) #this is the number that go to severe don't seek care
    
    shedding.cases.mU.y <- rbinom(1, y["Iy_m_syU"], Pr.mU.shed.y) #the number who move to moderate shedding
    m.who.will.abx.y <- round(y["Iy_m_syT"]*param["q"],0) #the number of moderates who seek care who will eventually receive abx#are treating q like a proportion, though this not actually correct. the math works though and conceptually is ok via prop.m.abx equation
    shedding.cases.mT.y <- rbinom(1, y["Iy_m_syT"]-m.who.will.abx.y, Pr.mT.shed.y) #the number of moderates who seek care who don't receive abx and move to shed
    abx.cases.m.y <- rbinom(1, m.who.will.abx.y, Pr.mT.abx.y) #the number of moderates who receive abx at this step
    shedding.cases.s.y <- rbinom(1, y["Iy_s_syU"], Pr.sU.shed.y) #the number of young severes who don't seek care and move to shed
    abx.cases.s.y <- rbinom(1, y["Iy_s_syT"], Pr.sT.abx.y) #the number of young severes who seek care and receive abx
    
    nat.recovered.a.y <- rbinom(1, y["Iy_a_sh"], Pr.a.recover.y) #the number of children asymptomatics who recover
    m.who.will.die.y <- round(y["Iy_m_sh"]*param["theta.m"],0) #this is the whole number of moderates who will eventually die
    nat.recovered.m.y <- rbinom(1, y["Iy_m_sh"]-m.who.will.die.y, Pr.m.shed.recover.y) #this is the actual number of I_sh who naturally recover. this incorporates the stochastic prob that will actually move in this step #rbinom requires whole numbers
    died.m.y <- rbinom(1, m.who.will.die.y, Pr.m.shed.die.y) #the number of moderates who move from I_sh to dead
    abx.recovered.m.y <- rbinom(1, y["Iy_m_abx"], Pr.m.abx.recover.y) #number of moderates who move from abx to recover in this step
    s.who.will.die.y <- round(y["Iy_s_sh"]*param["theta.s"],0) #this is the whole number of severes who will die
    nat.recovered.s.y <- rbinom(1, y["Iy_s_sh"]-s.who.will.die.y, Pr.s.shed.recover.y) #number of severes who move from Is_sh to recover w/o abx
    died.s.y <- rbinom(1, s.who.will.die.y, Pr.s.shed.die.y) #number of Is_sh who die in this time step
    abx.recovered.s.y <- rbinom(1, y["Iy_s_abx"], Pr.s.abx.recover.y) #number of severes who recover after abx in this time step
    
    #Find the deltas for each compartment
    ## Children
    dSy <- -incident.aging - incident.exposed.y
    dEy <- incident.exposed.y - incident.cases.a.y
    dIa.sh.y <-incident.cases.a.y - nat.recovered.a.y
    
    dIm.syU.y <- incident.cases.symp.mU.y - shedding.cases.mU.y
    dIm.sh.y <- shedding.cases.mU.y + shedding.cases.mT.y - nat.recovered.m.y - died.m.y
    dIm.syT.y <- incident.cases.symp.mT.y - shedding.cases.mT.y - abx.cases.m.y
    dIm.abx.y <- abx.cases.m.y - abx.recovered.m.y
    
    dIs.syU.y <- incident.cases.symp.sU.y - shedding.cases.s.y
    dIs.sh.y <- shedding.cases.s.y - nat.recovered.s.y - died.s.y
    dIs.syT.y <- incident.cases.symp.sT.y - abx.cases.s.y
    dIs.abx.y <- abx.cases.s.y - abx.recovered.s.y
    
    ## Adults
    dS <- incident.aging - incident.exposed
    dE <- incident.exposed - incident.cases.a - 
      incident.cases.symp.mU - incident.cases.symp.mT - 
      incident.cases.symp.sU - incident.cases.symp.sT
    dIa.sh <- incident.cases.a - nat.recovered.a
    
    dIm.syU <- incident.cases.symp.mU - shedding.cases.mU
    dIm.sh <- shedding.cases.mU + shedding.cases.mT - nat.recovered.m - died.m
    dIm.syT <- incident.cases.symp.mT - shedding.cases.mT - abx.cases.m
    dIm.abx <- abx.cases.m - abx.recovered.m
    
    dIs.syU <- incident.cases.symp.sU - shedding.cases.s
    dIs.sh <- shedding.cases.s - nat.recovered.s - died.s
    dIs.syT <- incident.cases.symp.sT - abx.cases.s
    dIs.abx <- abx.cases.s - abx.recovered.s
    
    ##Recovery/Death
    dR <- nat.recovered.a + nat.recovered.a.y + nat.recovered.m + nat.recovered.m.y + nat.recovered.s + nat.recovered.s.y
    dD <- died.m + died.m.y + died.s + died.s.y
    dR.abx <- abx.recovered.m + abx.recovered.m.y + abx.recovered.s + abx.recovered.s.y
    
    sum(dIa.sh.y,dIa.sh,dIm.syU.y,dIm.syU,dIm.syT.y,dIm.syT,dIs.syU.y,dIs.syU,dIs.syT.y,dIs.syT)
    
    deltas <- c(dSy,dS,dEy,dE,dIa.sh.y,dIa.sh,dIm.syU.y,dIm.syU,dIm.sh.y,dIm.sh,dIm.syT.y,dIm.syT,dIm.abx.y,dIm.abx,dIs.syU.y,dIs.syU,dIs.sh.y,dIs.sh,dIs.syT.y,dIs.syT,dIs.abx.y,dIs.abx,dR,dD,dR.abx) # calculate step sizes. this is net change for each box
    z <- c(incident.aging,
           incident.exposed, incident.exposed.y,
           incident.cases.a, incident.cases.a.y,
           incident.cases.symp.mU, incident.cases.symp.mU.y, incident.cases.symp.mT, incident.cases.symp.mT.y,
           incident.cases.symp.sU, incident.cases.symp.sU.y, incident.cases.symp.sT, incident.cases.symp.sT.y,
           shedding.cases.mU, shedding.cases.mU.y, shedding.cases.mT, shedding.cases.mT.y, abx.cases.m, abx.cases.m.y,
           shedding.cases.s, shedding.cases.s.y, abx.cases.s, abx.cases.s.y,
           nat.recovered.a, nat.recovered.a.y,
           nat.recovered.m, nat.recovered.m.y, died.m, died.m.y, abx.recovered.m, abx.recovered.m.y,
           nat.recovered.s, nat.recovered.s.y, died.s, died.s.y, abx.recovered.s, abx.recovered.s.y)
    #this is the incident (new) number of ppl along each arrow at each step
    
    y <- y+deltas 
    # print(c(t,y))
    
    #trick to speed up the code
    if (!final.only){
      SEIIR.output<-rbind(SEIIR.output, c(t,y))
      incidences<-rbind(incidences, c(t,z))
    }
  }
  
  if(final.only){
    SEIIR.output[1,]<-c(t,y)
    incidences[1,]<-c(t,z)
  }
  
  combo<-list(SEIIR.output=as.data.frame(SEIIR.output),incidences=as.data.frame(incidences))
  return(combo)
}

