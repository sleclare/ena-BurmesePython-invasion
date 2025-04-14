setwd("./data_files/data")
####set up####
#clear environment
rm(list = ls())

#load packages#
#archive web source to install enaR: https://cran.r-project.org/src/contrib/Archive/enaR/
#to install enaR, first need to install packages 'stringr', 'sna', 'network', 'gdata', & 'limSolve'
#Code: install.packages("enaR_3.0.0.tar.gz", repos=NULL, type="source")
#When new version of gdata is installed: 
library(enaR)
library(cheddar)

#Read in original Networks
origram <- read.enam(file = "GRAMDRY.xlsx")
oricyp <- read.enam(file = "CYPDRY.xlsx")
orimang <- read.enam(file = "MANGDRY.xlsx")

#unpack original network
namemat <- c('Flow','inp','resp','exp','out','biom', 'living') #create names for the matrices in the order that they appear in the network data
for(i in 1:7){
  dat <- unpack(origram)[[i]] #pulls the relevant data from the network data. Note the double square bracket here, needed because of the list of lists. 
  assign(x = namemat[i], dat) 
}

#Read in Python Diet Data
diet <- read.csv("pythondiet.csv", header = T, stringsAsFactors = F)
#Read in average mass data by compartment
mass_by_compartment <- read.csv("compartment_mass.csv")

#Community data format for cheddar package derived metrics
nodes <- read.csv("nodes.csv", header = T)
properties <- read.csv("properties.csv")
properties <- as.list(properties)

####Functions####
#Estimates biomass in grams of carbon
carbon_biomass <- function(df, density, d, p) {
  avg_mass <- mean(df)
  carbon_mass <- avg_mass*density*d*p
  return(carbon_mass)
}
#Estimates annual consumption rate in grams of carbon
est_acr <- function(a, b, mass, fmr, density, p, conversion = 1){
  daily <- (a*mass^b)/fmr #gdw/day/ind.
  acr <- daily*density*p*c #annual consumption rate in gC/yr/m^2
  return(acr)
} 
#calculates disparity proportion
disprop <- function(Flow, inp, resp, exp){
  input <- apply(Flow, 2, sum) + inp
  output <- apply(Flow, 1, sum) + resp + exp
  return(input/output)
}
#Integrates new compartment into existing network
adjnet <- function(net, c.mass, acr, prey.prop, predator.prop, production, egestion,
                   fmr, pmr, inflows, outflows, newbiom = NULL, biom.est = NULL, st.adj){
  #original network
  namemat <- c('Flow','inp','resp','exp','out','biom', 'living')  #create names for the matrices in the order that they appear in the network data
  for(i in 1:7){
    dat <- unpack(origram)[[i]]  #pulls the relevant data from the network data. Note the double square bracket here, needed because of the list of lists. 
    assign(x = namemat[i], dat) 
  }
  
  ###import new compartment####
  
  ##Enter in a 'empty' scenario compartment
  python <- vector(mode = "numeric", length = length(biom)) #length = original number of compartments
  Flow <- cbind(python, Flow)
  python <- vector(mode = "numeric", length = (length(biom)+1)) #length = original number of compartments + scenario compartment
  Flow <- rbind(python, Flow)
  ##enter in new value to existing component vectors for scenario compartment
  inp <- c(0, inp)
  resp <- c(0, resp)
  exp <- c(0, exp)
  out <- c(0, out)
  biom <- c(0, biom)
  living <- c(TRUE, living)
  
  #original component metrics prior to adjustments
  origflow <- Flow
  origbiom <- biom
  originp <- inp
  origresp <- resp
  origexp <- exp
  
  #disparity proportion
  x <- disprop(Flow, inp, resp, exp)
  lowlim <- min(x, na.rm = T)
  uplim <- max(x, na.rm = T)
  
  ##Enter in value estimates to scenario compartment####
  ##biomass estimate
  biom[1] <- c.mass
  
  ##flow matrix estimates
  ##Input flows - prey proportions
  prey <- prey.prop*acr
  
  for(i in 1:length(inflows)){
    Flow[inflows[i],1] <- prey[i]
  }
  
  if(length(inflows) != length(prey)) {
    inp[1] <- prey[length(prey)]
  }
  ##output flows - predation and detritus
  predator <- acr*production*pmr*predator.prop
  
  for(i in 1:length(outflows)){
    Flow[1, outflows[i]] <- predator[i]
  }
  
  if(length(outflows) != length(predator)){
    exp[1] <- predator[length(predator)]
  }
  
  #Detritus outflows
  Flow[1,ncol(Flow)-2] <- acr*egestion
  Flow[1,ncol(Flow)] <- acr*production*(1-pmr)
  
  #respiration estimates
  resp[1] <- acr*fmr
  out[1] <- resp[1]+exp[1]
  
  C_1 <- sum(Flow[,1], inp[1])/biom[1] #consumption/biomass ratio derived from the original network
  C_1_flow <- sum(Flow[,1])/sum(Flow[,1], inp[1])
  C_1_in <- inp[1]/sum(Flow[,1], inp[1])
  
  P_1 <- sum(Flow[1,],exp[1])/biom[1] #production/biomass ratio derived from the original network
  P_1_flow <- sum(Flow[1,])/sum(Flow[1,],exp[1])
  P_1_ex <- exp[1]/sum(Flow[1,],exp[1])
  
  R_1 <- resp[1]/biom[1] #respiration/biomass ratio derived from the original network
  
  ###Adjusting biomass compartments for novel outputs & recalculating flows####
  
  if(is.null(newbiom) == FALSE){
    newbiom <- newbiom + 1
    for (i in 1:length(newbiom)) {
      #obtain flow proportions and ecoystem ratios
      C <- sum(origflow[,newbiom[i]], originp[newbiom[i]])/origbiom[newbiom[i]] #consumption/biomass ratio derived from the original network
      C.flow <- sum(origflow[,newbiom[i]])/sum(origflow[,newbiom[i]], originp[newbiom[i]])
      C.in <- 1-C.flow 
      
      P <- sum(origflow[newbiom[i],],origexp[newbiom[i]])/origbiom[newbiom[i]] #production/biomass ratio derived from the original network
      P.flow <- sum(origflow[newbiom[i],])/sum(origflow[newbiom[i],],origexp[newbiom[i]])
      P.ex <- 1-P.flow
      
      R <- origresp[newbiom[i]]/origbiom[newbiom[i]] #respiration/biomass ratio derived from the original network
      
      biom[newbiom[i]] <- biom.est[i]
      resources <- which(Flow[,newbiom[i]]>0)
      biomass <- biom[resources]
      p.ratio <- vector(length = length(resources))
      if(length(p.ratio) > 0){
        for (j in 1:length(p.ratio)) {
          if(resources[j]==1){
            p.ratio[j] <- P_1
          } else {
            p.ratio[j] <- sum(origflow[resources[j],],origexp[resources[j]])/origbiom[resources[j]] 
          }
        }
        inprop <- (biomass*p.ratio)/sum(biomass*p.ratio)  
      } else{
        inprop <- 0
      }
      
      #adjusting consumer proportions
      consumers <- which(Flow[newbiom[i],]>0)
      det.flow.na <- vector( length = length(consumers))
      for (j in 1:length(det.flow.na)) {
        if(consumers[j] == 65){
          det.flow.na[j] <- Flow[newbiom[i],consumers[j]]/sum(Flow[i,])
        } else{
          if(consumers[j] == 66){
            det.flow.na[j] <- Flow[newbiom[i],consumers[j]]/sum(Flow[i,])
          } else{
            if(consumers[j] == 67) {
              det.flow.na[j] <- Flow[newbiom[i],consumers[j]]/sum(Flow[i,])
            } else{
              det.flow.na[j] <- NA
            }
          }
        }
      }
      det.flow <- na.omit(det.flow.na)
      pred.flow <- 1 - sum(Flow[newbiom[i], (ncol(Flow)-2):(ncol(Flow))]/sum(Flow[newbiom[i],]))
      if(pred.flow > 1e-15){
        predators <- consumers[1:(length(det.flow.na) - length(det.flow))]
        biomass <- biom[predators]
        c.ratio <- vector(length = length(predators))
        for (j in 1:length(c.ratio)) {
          if(predators[j]==1){
            c.ratio[j] <- C_1
          }else{
            c.ratio[j] <- sum(origflow[,predators[j]], originp[predators[j]])/origbiom[predators[j]]
          }
        }
        outprop <- c((pred.flow*((biomass*c.ratio)/sum(biomass*c.ratio))),det.flow)
      } else{
        outprop <- det.flow 
      }
      
      #Flow Recalculation
      Flow[resources,newbiom[i]] <- inprop*C*biom[newbiom[i]]*C.flow #new ingoing flows
      inp[newbiom[i]] <- C*biom[newbiom[i]]*C.in #new input flows from outside system
      Flow[newbiom[i],consumers] <- outprop*P*biom[newbiom[i]]*P.flow #new outgoing flows
      exp[newbiom[i]] <- P*biom[newbiom[i]]*P.ex #new export estimates
      resp[newbiom[i]] <- R*biom[newbiom[i]] #new respiration estimates
    }
  }
  
  #disparity proportion
  x <- disprop(Flow, inp, resp, exp); x
  if(any(x>uplim | x<lowlim)){
    for (n in 1:2) {
      if(is.na(any(x < lowlim)) == FALSE){
        if(any(x<lowlim)){
          newbiom <- which(x<lowlim)
          for (i in newbiom){
            #Ecosystem ratios
            if(i==1){
              C <- C_1
              C.flow <- C_1_flow
              C.in <- C_1_in
              
              P <- P_1
              P.flow <- P_1_flow
              P.ex <- P_1_ex
              
              R <- R_1
            }else{
              C <- sum(origflow[,i], originp[i])/origbiom[i] #consumption/biomass ratio derived from the original network
              C.flow <- sum(origflow[,i])/sum(origflow[,i], originp[i])
              C.in <- 1-C.flow 
              
              P <- sum(origflow[i,],origexp[i])/origbiom[i] #production/biomass ratio derived from the original network
              P.flow <- sum(origflow[i,])/sum(origflow[i,],origexp[i])
              P.ex <- 1-P.flow
              
              R <- origresp[i]/origbiom[i] #respiration/biomass ratio derived from the original network
            }
            #adjusted biomass
            biom[i] <- biom[i]-(biom[i]*st.adj)
            #adjusting resources proportions
            resources <- which(Flow[,i]>0)
            biomass <- biom[resources]
            p.ratio <- vector(length = length(resources))
            if(length(p.ratio) > 0){
              for (j in 1:length(p.ratio)) {
                if(resources[j]==1){
                  p.ratio[j] <- P_1
                } else {
                  p.ratio[j] <- sum(origflow[resources[j],],origexp[resources[j]])/origbiom[resources[j]] 
                }
              }
              inprop <- (biomass*p.ratio)/sum(biomass*p.ratio)  
            } else{
              inprop <- 0
            }
            
            #adjusting consumer proportions
            consumers <- which(Flow[i,]>0)
            det.flow.na <- vector( length = length(consumers))
            for (j in 1:length(det.flow.na)) {
              if(consumers[j] == 65){
                det.flow.na[j] <- Flow[i,consumers[j]]/sum(Flow[i,])
              } else{
                if(consumers[j] == 66){
                  det.flow.na[j] <- Flow[i,consumers[j]]/sum(Flow[i,])
                } else{
                  if(consumers[j] == 67) {
                    det.flow.na[j] <- Flow[i,consumers[j]]/sum(Flow[i,])
                  } else{
                    det.flow.na[j] <- NA
                  }
                }
              }
            }
            det.flow <- na.omit(det.flow.na)
            pred.flow <- 1 - sum(Flow[i, (ncol(Flow)-2):(ncol(Flow))]/sum(Flow[i,]))
            if(pred.flow > 1e-15){
              predators <- consumers[1:(length(det.flow.na) - length(det.flow))]
              biomass <- biom[predators]
              c.ratio <- vector(length = length(predators))
              for (j in 1:length(c.ratio)) {
                if(predators[j]==1){
                  c.ratio[j] <- C_1
                }else{
                  c.ratio[j] <- sum(origflow[,predators[j]], originp[predators[j]])/origbiom[predators[j]]
                }
              }
              outprop <- c((pred.flow*((biomass*c.ratio)/sum(biomass*c.ratio))),det.flow)
            } else{
              outprop <- det.flow 
            }
            
            #Flow Recalculation
            Flow[resources,i] <- inprop*C*biom[i]*C.flow #new ingoing flows
            inp[i] <- C*biom[i]*C.in #new input flows from outside system
            Flow[i,consumers] <- outprop*P*biom[i]*P.flow #new outgoing flows
            exp[i] <- P*biom[i]*P.ex #new export estimates
            resp[i] <- R*biom[i] #new respiration estimates
          }
        } 
      }
      if(any(x>uplim)){
        newbiom <- which(x>uplim)
        for (i in newbiom){
          #Ecosystem ratios
          if(i==1){
            C <- C_1
            C.flow <- C_1_flow
            C.in <- C_1_in
            
            P <- P_1
            P.flow <- P_1_flow
            P.ex <- P_1_ex
            
            R <- R_1
          }else{
            C <- sum(origflow[,i], originp[i])/origbiom[i] #consumption/biomass ratio derived from the original network
            C.flow <- sum(origflow[,i])/sum(origflow[,i], originp[i])
            C.in <- 1-C.flow 
            
            P <- sum(origflow[i,],origexp[i])/origbiom[i] #production/biomass ratio derived from the original network
            P.flow <- sum(origflow[i,])/sum(origflow[i,],origexp[i])
            P.ex <- 1-P.flow
            
            R <- origresp[i]/origbiom[i] #respiration/biomass ratio derived from the original network
          }
          #adjusted biomass
          biom[i] <- biom[i]+(biom[i]*st.adj)
          #adjusting resources proportions
          resources <- which(Flow[,i]>0)
          biomass <- biom[resources]
          p.ratio <- vector(length = length(resources))
          if(length(p.ratio) > 0){
            for (j in 1:length(p.ratio)) {
              if(resources[j]==1){
                p.ratio[j] <- P_1
              } else {
                p.ratio[j] <- sum(origflow[resources[j],],origexp[resources[j]])/origbiom[resources[j]] 
              }
            }
            inprop <- (biomass*p.ratio)/sum(biomass*p.ratio)  
          } else{
            inprop <- 0
          }
          
          #adjusting consumer proportions
          consumers <- which(Flow[i,]>0)
          det.flow.na <- vector( length = length(consumers))
          for (j in 1:length(det.flow.na)) {
            if(consumers[j] == 65){
              det.flow.na[j] <- Flow[i,consumers[j]]/sum(Flow[i,])
            } else{
              if(consumers[j] == 66){
                det.flow.na[j] <- Flow[i,consumers[j]]/sum(Flow[i,])
              } else{
                if(consumers[j] == 67) {
                  det.flow.na[j] <- Flow[i,consumers[j]]/sum(Flow[i,])
                } else{
                  det.flow.na[j] <- NA
                }
              }
            }
          }
          det.flow <- na.omit(det.flow.na)
          pred.flow <- 1 - sum(Flow[i, (ncol(Flow)-2):(ncol(Flow))]/sum(Flow[i,]))
          if(pred.flow > 1e-15){
            predators <- consumers[1:(length(det.flow.na) - length(det.flow))]
            biomass <- biom[predators]
            c.ratio <- vector(length = length(predators))
            for (j in 1:length(c.ratio)) {
              if(predators[j]==1){
                c.ratio[j] <- C_1
              }else{
                c.ratio[j] <- sum(origflow[,predators[j]], originp[predators[j]])/origbiom[predators[j]]
              }
            }
            outprop <- c((pred.flow*((biomass*c.ratio)/sum(biomass*c.ratio))),det.flow)
          } else{
            outprop <- det.flow 
          }
          
          #Flow Recalculation
          Flow[resources,i] <- inprop*C*biom[i]*C.flow #new ingoing flows
          inp[i] <- C*biom[i]*C.in #new input flows from outside system
          Flow[i,consumers] <- outprop*P*biom[i]*P.flow #new outgoing flows
          exp[i] <- P*biom[i]*P.ex #new export estimates
          resp[i] <- R*biom[i] #new respiration estimates
        }
      }
      x <- disprop(Flow, inp, resp, exp); x
    }
    
    adjmodel <- pack(flow = Flow, input = inp, respiration = resp, 
                     export = exp, output = out, 
                     storage = biom, living = living)
    adjmodel <- balance(adjmodel, method = "O")
    
    return(adjmodel)
  } else {
    adjmodel <- pack(flow = Flow, input = inp, respiration = resp, 
                     export = exp, output = out, 
                     storage = biom, living = living)
    adjmodel <- balance(adjmodel, method = "O")
    
    return(adjmodel)
  }
}
#Calculates individual ascendency coefficients for each compartment
asc_coef <- function(net){
  T.ulan <- as.extended(net)
  N <- ncol(T.ulan) # set up N
  ami <- mat.or.vec(N,N) # initialize ascendency matrix
  cap <- mat.or.vec(N,N) # initialize capacity matrix
  
  #Calculate Total System Throughput
  TSTp <- sum(T.ulan)
  
  ## H = Total Flow Diversity
  h <- T.ulan/sum(T.ulan)
  h2 <- log2(h)
  h2[!is.finite(h2)] <- 0
  H = - sum(h * h2)   # Total Flow Diversity
  
  CAP <- H * TSTp     # Capactity
  
  # loop through T.ulan to calculate AMI
  for (i in 1:N){
    for (j in 1:N){
      if (T.ulan[i,j] == 0){
        ami[i,j] <- 0
      }else{
        ami[i,j] <- T.ulan[i,j]/TSTp * log2((T.ulan[i,j]*TSTp)/(sum(T.ulan[,j])*sum(T.ulan[i,])))
      }
    }
  }
  #individual contributions sensitivities
  A <- apply(ami, 2, sum) #contribution of each taxon to overall system performance
  Tp <- apply(T.ulan, 1, sum) #Trhoughput for each taxon
  
  Sens. <- (A/Tp)*TSTp #Ascendency Sensitivity Coefficients; Constribution that each taxon makes to the total system
  Sens. <- Sens.[-c((ncol(T.ulan)-2):ncol(T.ulan))]
  return(Sens.)
}
#Calculates throughflow and closeness centrality for each compartment
centrality <- function(net) {
  namemat <- c('Flow','inp','resp','exp','out','biom', 'living')  #create names for the matrices in the order that they appear in the network data
  for(i in 1:7){  #loops over the unpack function
    dat <- unpack(net)[[i]]  #pulls the relevant data from the network data. Note the double square bracket here, needed because of the list of lists. 
    assign(x = namemat[i], dat)  #assign the data to its own object. 
  }
  
  ##Centrality####
  #weighted betweenness centrality
  ccentrality <- closeness_w(Flow, directed = TRUE, gconly = FALSE)
  ccentrality <- ccentrality[,2]
  #Throughflow Centrality
  T. <- apply(Flow,1,sum) + inp;   # input throughflow (assuming steady state)
  TST <- sum(T.)  # total system throughflow
  
  T.out <- apply(Flow,2,sum) + out #output throughflow
  T.Central <- T.out/TST #Throughflow Centrality
  return(list('C.Centrality' = ccentrality, 'T.Centrality' = T.Central))
} 
#modified NodeQuantitativeDescriptor function from cheddar package, removing trophic chain metrics to reduce computational load
NodeDescript <- function(community, weight) {
  # Bersier et al (2002) Ecology
  
  if(!is.Community(community)) stop('Not a Community')
  
  # A common operation
  vlog2v <- function(v)   v*log2(v)
  
  b <- PredationMatrix(community, weight)
  bIn <- colSums(b)
  bOut <- rowSums(b)
  
  # Diversity of inflows and outflows
  HN <- -rowSums(vlog2v(t(b)/bIn), na.rm=TRUE)    # p 2397, eq 5
  HP <- -rowSums(vlog2v(b/bOut), na.rm=TRUE)      # p 2397, eq 6
  
  # Equivalent numbers of prey and predators
  nN <-2^HN          # p 2397, eq 7
  nN[0==bIn] <- 0
  nP <-2^HP          # p 2397, eq 8
  nP[0==bOut] <- 0
  
  d.prime <- nN/(nN + nP)         # p 2397, eq 9
  d <- bIn*nN/(bIn*nN + bOut*nP)  # p 2397, eq 10
  
  # Standardised quantitative G and V, p 2400, eq 28-29
  g.prime <- nN*NumberOfNodes(community) / sum(nN)
  g <- bIn*nN*NumberOfNodes(community) / sum(bIn*nN)
  v.prime <- nP*NumberOfNodes(community) / sum(nP)
  v <- bOut*nP*NumberOfNodes(community) / sum(bOut*nP)
  
  # Bersier et al table 1
  res <- cbind(NResources=NumberOfResources(community), 
               NConsumers=NumberOfConsumers(community), 
               bIn,
               bOut,
               nN, 
               nP, 
               d.prime, 
               d,
               g.prime, 
               g, 
               v.prime, 
               v)
  rownames(res) <- unname(NP(community, 'node'))
  return (res)
}
#modified QuantitativeDescriptor function from cheddar package, removing trophic chain metrics to reduce computational load
QuantitativeRevised <- function(community, weight, top.level.threshold=0.99){
  # A common operation
  vlog2v <- function(v)   v*log2(v)
  
  # Community-level descriptors are derived from the node-level decsriptors
  
  # TODO prevent double-computation of chains
  
  np <- NodeDescript(community, weight)
  
  b <- PredationMatrix(community, weight)
  sumb <- sum(b)
  
  # Diversity of inflows and outflows
  HN <- -rowSums(vlog2v(t(b)/np[,'bIn']), na.rm=TRUE)    # p 2397, eq 5
  HP <- -rowSums(vlog2v(b/np[,'bOut']), na.rm=TRUE)      # p 2397, eq 6
  
  # See comment at bottom of 2397 about about top-level threshold
  fracT.q.prime <- mean(np[,'d.prime']>=top.level.threshold)
  fracI.q.prime <- mean(0<np[,'d.prime'] & np[,'d.prime']<top.level.threshold)
  fracB.q.prime <- mean(0==np[,'d.prime'])
  
  fracT.q <- mean(np[,'d']>=top.level.threshold)
  fracI.q <- mean(0<np[,'d'] & np[,'d']<top.level.threshold)
  fracB.q <- mean(0==np[,'d'])
  
  # Ratios of resources to consumers - p 2398, eq 11 and 12
  NP.q.prime <- 2^(-sum(vlog2v(np[,'nP']/sum(np[,'nP'])), na.rm=TRUE)) / 
    2^(-sum(vlog2v(np[,'nN']/sum(np[,'nN'])), na.rm=TRUE))
  
  NP.q <- 2^(-sum(vlog2v((np[,'bOut']*np[,'nP'])/sum(np[,'bOut']*np[,'nP'])), na.rm=TRUE)) / 
    2^(-sum(vlog2v((np[,'bIn']*np[,'nN'])/sum(np[,'bIn']*np[,'nN'])), na.rm=TRUE))
  
  # Link properties
  # p 2398, eq 13
  LD.q.prime <- (sum(np[,'nP']) + sum(np[,'nN'])) / (2*NumberOfNodes(community))
  # p 2398, eq 14
  LD.q <- (sum(np[,'bOut']*np[,'nP']/sumb, na.rm=TRUE) + 
             sum(np[,'bIn']*np[,'nN']/sumb, na.rm=TRUE))/2
  
  # p 2398, col 2
  C.q.prime <- LD.q.prime/NumberOfNodes(community)
  C.q <- LD.q/NumberOfNodes(community)
  
  # ...the sum of the diversity of outflows weighted by the total outflows, 
  # of the diversity of inflows weighted by the total inflows. Phi can be 
  # thought of as the average amount of choice in trophic pathways (Ulanowicz 
  # and Wolff 1991)
  Phi <- sum(HP*np[,'bOut']/sumb) + sum(HN*np[,'bIn']/sumb) # p 2398, eq 16
  m <- 2^(Phi/2) # Effective connectance per node p 2398, eq 15
  
  # Page 2398, eq 18
  PhiAB <- function(A, B)
  {
    # A and B should be functions that take a community and return node 
    # names or indices.
    A <- A(community)
    B <- B(community)
    bOut <- rowSums(b)[A]
    bIn <- colSums(b)[B]
    return (sum((bOut/sumb) * -vlog2v(b[A,B]/bOut),  na.rm=TRUE) + 
              sum((bIn/sumb)  * -vlog2v(t(b[A,B])/bIn), na.rm=TRUE))
  }
  
  fracTI.q <- PhiAB(IntermediateNodes, TopLevelNodes) / Phi
  fracTB.q <- PhiAB(BasalNodes, TopLevelNodes) / Phi
  fracII.q <- PhiAB(IntermediateNodes, IntermediateNodes) / Phi
  fracIB.q <- PhiAB(BasalNodes, IntermediateNodes) / Phi
  
  # Page 2399, top left
  PhiPrime <- sum(HP/NumberOfNodes(community)) + sum(HN/NumberOfNodes(community))
  PhiABPrime <- function(A, B)
  {
    # A and B should be functions that take a community and return node 
    # names or indices.
    A <- A(community)
    B <- B(community)
    bOut <- rowSums(b)[A]
    bIn <- colSums(b)[B]
    s <- NumberOfNodes(community)
    return (sum((1/s) * -vlog2v(b[A,B]/bOut),  na.rm=TRUE) + 
              sum((1/s) * -vlog2v(t(b[A,B])/bIn), na.rm=TRUE))
  }
  
  fracTI.q.prime <- PhiABPrime(IntermediateNodes, TopLevelNodes) / PhiPrime
  fracTB.q.prime <- PhiABPrime(BasalNodes, TopLevelNodes) / PhiPrime
  fracII.q.prime <- PhiABPrime(IntermediateNodes, IntermediateNodes) / PhiPrime
  fracIB.q.prime <- PhiABPrime(BasalNodes, IntermediateNodes) / PhiPrime
  
  # Quantitative generality and vulnerability, p 2400, eq 24-27
  nT <- length(TopLevelNodes(community))
  nI <- length(IntermediateNodes(community))
  nB <- length(BasalNodes(community))
  G.q.prime <- sum(np[,'nN'])/(nT+nI)
  G.q <- sum(np[,'nN']*np[,'bIn']/sumb, na.rm=TRUE)
  V.q.prime <- sum(np[,'nP'])/(nI+nB)
  V.q <- sum(np[,'nP']*np[,'bOut']/sumb, na.rm=TRUE)
  
  # Bersier et al table 2
  tlps <- TLPS(community, node.properties=c('IsTopLevelNode', 'IsIntermediateNode', 'IsBasalNode'))
  fracTI <- with(tlps, sum(resource.IsIntermediateNode & consumer.IsTopLevelNode))/nrow(tlps)
  fracTB <- with(tlps, sum(resource.IsBasalNode & consumer.IsTopLevelNode))/nrow(tlps)
  fracII <- with(tlps, sum(resource.IsIntermediateNode & consumer.IsIntermediateNode))/nrow(tlps)
  fracIB <- with(tlps, sum(resource.IsBasalNode & consumer.IsIntermediateNode))/nrow(tlps)
  # This is far faster than ChainLengths
  Qualitative <- c(FractionTopLevelNodes(community), 
                   FractionIntermediateNodes(community), 
                   FractionBasalNodes(community), 
                   sum(NumberOfConsumers(community)>0) / sum(NumberOfResources(community)>0), 
                   LinkageDensity(community), 
                   DirectedConnectance(community), 
                   fracTI, 
                   fracTB, 
                   fracII, 
                   fracIB,
                   mean(TrophicGenerality(community)[NumberOfResources(community)>0]), 
                   mean(TrophicVulnerability(community)[NumberOfConsumers(community)>0]), 
                   sd(NormalisedTrophicGenerality(community)), 
                   sd(NormalisedTrophicVulnerability(community)))
  
  Unweighted <- c(fracT.q.prime, 
                  fracI.q.prime, 
                  fracB.q.prime, 
                  NP.q.prime, 
                  LD.q.prime,
                  C.q.prime, 
                  fracTI.q.prime, 
                  fracTB.q.prime, 
                  fracII.q.prime, 
                  fracIB.q.prime,
                  G.q.prime, 
                  V.q.prime, 
                  sd(np[,'g.prime']), 
                  sd(np[,'v.prime']))
  
  Weighted <- c(fracT.q, 
                fracI.q, 
                fracB.q, 
                NP.q, 
                LD.q,
                C.q, 
                fracTI.q, 
                fracTB.q, 
                fracII.q, 
                fracIB.q, 
                G.q, 
                V.q, 
                sd(np[,'g']), 
                sd(np[,'v']))
  res <- cbind(Qualitative, Unweighted, Weighted)
  
  rownames(res) <- c('Fraction top level', 
                     'Fraction intermediate', 'Fraction basal', 'Ratio resources:consumers', 
                     'Link density', 'Connectance', 'Fraction links top:intermediate', 
                     'Fraction links top:basal', 'Fraction links intermediate:intermediate', 
                     'Fraction links intermediate:basal', 'Generality', 'Vulnerability', 
                     'SD standardised generality', 'SD standardised vulnerability')
  return (res)
} 

####Estimating Python Parameters####
#Load Packages
library(plyr)
library(ggplot2)

#bp = Burmese python

###Biomass##
bp_mass_data <- diet[!is.na(diet$MassG), 'MassG']
bp_density <- 5.00e-06 #ind/m^2 (Snow et al, 2007)
bp_d <- 0.25 #conversion from gww to gdw for snakes (jorgensen et al, 1991)
bp_p <- 0.45 #conversion from gdw to gC for snakes (jorgensen et al, 1991)

#Biomass Estimate
bp_c_mass <- carbon_biomass(bp_mass_data, density = bp_density, d = bp_d, p = bp_p)

###Flow Matrix##
##Annual Consumption Rate

bp_slope <- 0.00865 #allometric slope for carnivorous reptiles (Nagy 2022)
bp_intercept <- 0.963 #allometric intercept for carnivorous reptiles (Naggy 2022)
bp_fmr <- 0.45 #Proportion of python energy budget that is allocated to field metabolic rate (Cox & Secor, 2007)
c <- 365 #conversion factor from day to yr.

#Annual Consumption Rate Estimate
bp_acr <- est_acr(a = bp_slope, b = bp_intercept, mass = mean(bp_mass_data), fmr = bp_fmr,
                  density = bp_density, p = bp_p, conversion = c) 

##proportion of prey items within diet

#sub-setting ID'd prey species from the data set and acquiring diet frequency
id_diet <- diet[diet$SciName != "Unidentifiable" & diet$SciName != "NotYetID",]
diet_freq <- count(id_diet, 'CommonName')
diet_observation <- sum(diet_freq$freq)
diet_proportion <- diet_freq$freq/diet_observation
#delineating each prey item into their respective network compartment
prey.name <- c("alligator","muskrat", "mice&rats", "rabbit", "raccoon", "oppossum", "otter", "W-T deer", "bobcat", "grebe", "bitterns",
               "ducks", "gruiformes", "passerines", "input")
pe.j <- c(diet_proportion[1],diet_proportion[64], 
          sum(diet_proportion[7],diet_proportion[13],diet_proportion[18],diet_proportion[26],
              diet_proportion[37],diet_proportion[38],diet_proportion[40],diet_proportion[41],
              diet_proportion[49],diet_proportion[72]),
          sum(diet_proportion[24], diet_proportion[48],diet_proportion[71]),
          diet_proportion[59],diet_proportion[74],diet_proportion[62],diet_proportion[77],diet_proportion[11], 
          sum(diet_proportion[57],diet_proportion[70]),sum(diet_proportion[2],diet_proportion[43]), 
          sum(diet_proportion[9],diet_proportion[22],diet_proportion[33],diet_proportion[47],
              diet_proportion[51], diet_proportion[56], diet_proportion[81]),
          sum(diet_proportion[3],diet_proportion[15],diet_proportion[16], diet_proportion[42], 
              diet_proportion[44], diet_proportion[58], diet_proportion[75]),
          sum(diet_proportion[10], diet_proportion[17], diet_proportion[25], diet_proportion[28], 
              diet_proportion[39], diet_proportion[50], diet_proportion[55], diet_proportion[61], diet_proportion[65]))
gdw.conv <- 0.35 #conversion from grams of wet weight to grams of dry weight derived from ATLSS
gdw.j <- c(mass_by_compartment$avg.mass[1],mass_by_compartment$avg.mass[2],mass_by_compartment$avg.mass[3], 
           mass_by_compartment$avg.mass[4], mass_by_compartment$avg.mass[5],mass_by_compartment$avg.mass[6],
           mass_by_compartment$avg.mass[7],mass_by_compartment$avg.mass[8], mass_by_compartment$avg.mass[9],
           mass_by_compartment$avg.mass[10],mass_by_compartment$avg.mass[11],mass_by_compartment$avg.mass[12],
           mass_by_compartment$avg.mass[13],mass_by_compartment$avg.mass[14])*gdw.conv
w.inflow.prop <- pe.j*gdw.j
#proportion of prey outside of network consolidated to input compartment
input.pe <- c(diet_proportion[4],diet_proportion[5],sum(diet_proportion[6],diet_proportion[34]),
              diet_proportion[8],diet_proportion[12],diet_proportion[14],diet_proportion[21],diet_proportion[23],
              sum(diet_proportion[27],diet_proportion[30],diet_proportion[73]),diet_proportion[29],
              diet_proportion[31],diet_proportion[32],diet_proportion[36],diet_proportion[45],diet_proportion[46],
              diet_proportion[52],diet_proportion[53],sum(diet_proportion[54],diet_proportion[68]),diet_proportion[60],
              diet_proportion[63],diet_proportion[66],diet_proportion[67],diet_proportion[69],diet_proportion[76],
              diet_proportion[78],diet_proportion[79],diet_proportion[80],diet_proportion[82])
gdw.input <- c(mass_by_compartment$avg.mass[15],mass_by_compartment$avg.mass[16],mass_by_compartment$avg.mass[17],
               mass_by_compartment$avg.mass[18],mass_by_compartment$avg.mass[19],mass_by_compartment$avg.mass[20],
               mass_by_compartment$avg.mass[23],mass_by_compartment$avg.mass[24],mass_by_compartment$avg.mass[25],
               mass_by_compartment$avg.mass[26],mass_by_compartment$avg.mass[27],mass_by_compartment$avg.mass[28],
               mass_by_compartment$avg.mass[30],mass_by_compartment$avg.mass[31],mass_by_compartment$avg.mass[32],
               mass_by_compartment$avg.mass[33],mass_by_compartment$avg.mass[34],mass_by_compartment$avg.mass[35],
               mass_by_compartment$avg.mass[36],mass_by_compartment$avg.mass[37],mass_by_compartment$avg.mass[38],
               mass_by_compartment$avg.mass[39],mass_by_compartment$avg.mass[40],mass_by_compartment$avg.mass[41],
               mass_by_compartment$avg.mass[42],mass_by_compartment$avg.mass[43],mass_by_compartment$avg.mass[44],
               mass_by_compartment$avg.mass[45])*gdw.conv

w.input.prop <- sum(input.pe*gdw.input)
#proportion by Network Compartment
w.prop <- c(w.inflow.prop,w.input.prop)

bp_prey_prop <- w.prop/sum(w.prop)

#Visualization of proportion by Network Compartment
Bp.prey.comprop <- data.frame(prey.name, bp_prey_prop)

Bp.prey.comprop$prey.name <- factor(Bp.prey.comprop$prey.name, levels = Bp.prey.comprop$prey.name[order(Bp.prey.comprop$bp_prey_prop)])
prop_graph <- Bp.prey.comprop
prop_graph$inflows <- Bp.prey.comprop$bp_prey_prop*bp_acr
ggplot(prop_graph, aes(x=inflows, y=prey.name)) + theme_bw() + geom_bar(stat = "identity") +
  ggtitle("Distribution of carbon inflows for the Burmese Python by Prey Node") +
  xlab("Carbon Inflows (gC/m2/yr)") + ylab("Nodes")+
  theme(axis.text.y = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 15))+
  theme(axis.text=element_text(size=20))

##proportion of predation/outflows for the Burmese python compartment

#Predation rate estimates
bp_production <- 0.35 #Production rate (Cox & Secor, 2007) - Biomass Turnover that flows into other compartments based on balance assumption of network model
bp_pmr <- 0.45 #Proportion of mortality attributed to predation (Pittman & Bartoszek, 2021)
bp_egestion <- 0.20 #Egestion rate (Cox & Secor, 2007)

#Calculating predation proportions
predator.name <- c("snakes", "alligators", "bobcat","export")

#Outflows
within.predator.biomass <- c(biom[43], biom[45], biom[54])
within.predator.consumption <- c(sum(Flow[,43])/biom[43],sum(Flow[,45])/biom[45], sum(Flow[,54])/biom[54]) #Consumption/biomass ratio of each compartment taken from original network

#Exports from predators in neighboring ecosystems
#cypress
cyp.namemat <- c('c.Flow','c.inp','c.resp','c.exp','c.out','c.biom', 'c.living')
for(i in 1:7){
  dat <- unpack(oricyp)[[i]]
  assign(x = cyp.namemat[i], dat) 
}
#mangrove
mang.namemat <- c('m.Flow','m.inp','m.resp','m.exp','m.out','m.biom', 'm.living')
for(i in 1:7){
  dat <- unpack(orimang)[[i]]
  assign(x = mang.namemat[i], dat) 
}

export.predation <- sum(c.biom[39]*(sum(c.Flow[,39])/c.biom[39]), m.biom[63]*(sum(m.Flow[,63])/m.biom[63]), m.biom[70]*(sum(m.Flow[,70])/m.biom[70]), 
                        m.biom[76]*(sum(m.Flow[,76])/m.biom[76]))
predator.consumption <- c(within.predator.biomass*within.predator.consumption, export.predation)
#proportion of predator outflows based on biomass consumption values of determined predators
bp_predation_prop <- predator.consumption/sum(predator.consumption)

#Visualization of proporiton by Network Compartment
Bp.predator.comprop <- data.frame(predator.name, bp_predation_prop)

Bp.predator.comprop$predator.name <- factor(Bp.predator.comprop$predator.name, levels = Bp.predator.comprop$predator.name[order(Bp.predator.comprop$bp_predation_prop)])
Bp.predator.comprop$outflows <- Bp.predator.comprop$bp_predation_prop*bp_pmr*bp_acr*bp_production
Bp.sed <- bp_egestion*bp_acr
Bp.det <- (1 - bp_pmr)*bp_acr*bp_production
Bp.resp <- bp_fmr*bp_acr
outflow.name <- c("snakes", "alligators", "bobcat","export", "sediment carbon", "refractory detritus", "respiration")
outflow_dist <- c(Bp.predator.comprop$outflows,Bp.sed,Bp.det,Bp.resp)
Bp.outflows <- as.data.frame(outflow_dist)
Bp.outflows$name <- outflow.name

ggplot(Bp.outflows, aes(x=outflow_dist, y=reorder(name, outflow_dist))) + theme_bw() + geom_bar(stat = "identity") +
  xlab("Carbon Outflows (gC/m2/yr)") + ylab("Nodes")+
  theme(axis.text.y = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 15))+
  theme(axis.text=element_text(size=20))

#Position of each flow compartment
bp_inflows <- c(46:52, 54:55, 57:59, 62, 64) #position of each prey compartment after integration of additional compartment
bp_outflows <- c(44, 46, 55) #position of each predator compartment after integration of additional compartment

newbiom <- c(48,49,50,53,54) #Meso-Mammal compartments
biom.est <- c(length(newbiom))
for (i in 1:length(newbiom)) {
  if(newbiom[i]==48){
    biom.est[i] <- biom[newbiom[i]]-(biom[newbiom[i]]*0.999)
  } else{
    if(newbiom[i]==49){
      biom.est[i] <- biom[newbiom[i]]-(biom[newbiom[i]]*0.993)
    } else{
      if(newbiom[i]==50){
        biom.est[i] <- biom[newbiom[i]]-(biom[newbiom[i]]*0.989)
      } else{
        if(newbiom[i]==53){
          biom.est[i] <- biom[newbiom[i]]-(biom[newbiom[i]]*0.941)
        } else{
          if(newbiom[i]==54){
            biom.est[i] <- biom[newbiom[i]]-(biom[newbiom[i]]*0.875)
          }
        }
      }
    }
  }
}

invgram <- adjnet(net = origram, c.mass = bp_c_mass, acr = bp_acr, prey.prop = bp_prey_prop,
                  predator.prop = bp_predation_prop , production = bp_production, 
                  egestion = bp_egestion, fmr = bp_fmr, pmr = bp_pmr, inflows = bp_inflows, 
                  outflows = bp_outflows, newbiom = newbiom, biom.est = biom.est,st.adj=0.1)

#unpack invaded network for analyses
namemat <- c('post_Flow','post_inp','post_resp','post_exp','post_out','post_biom', 'post_living')
for(i in 1:7){
  dat <- unpack(invgram)[[i]]
  assign(x = namemat[i], dat)
}

####Analyses########################
library(network)
library(dplyr)
library(reshape2)
library(gt)
library(tnet)
library(mvnormtest)
library(MVN)
library(raster)
library(cluster)
library(vegan)
library(factoextra)
library(circlize)
library(MASS)
library(janitor)

#Trophic Analysis####

#Lindeman Trophic Analysis
#Function (Get source)
plot.lindeman <- function(x = 'model', enatroagg="troagg", primprod, type = 1){
  if (class(x) != 'network'){warning('x is not a network class object')}
  
  # define primary producers & nonliving (by name)
  u <- unpack(x)
  vn <- x%v%'vertex.names'
  nonliving = which(u$living == FALSE)
  
  ## apply Trophic Aggregation
  if (enatroagg == "troagg"){
    enatroagg <- enaTroAgg(x)
  }else{}
  
  ns <- as.data.frame(enatroagg$ns)
  ## primprod?
  if (exists("primprod")){
    warning("Please supply a vector of primary producers as 'primprod'.")
  }else{}
  
  ## count number of compartments for LS
  ntl = length(enatroagg$GC[enatroagg$GC])
  
  # text scaling  -  to make nicer plots when NTL is large
  if(ntl <= 5){
    lvl.1 = 1.3
    lvl.2 = 1
    lvl.3 = 0.95
    lvl.4 = 1
  } else {
    lvl.1 = 1
    lvl.2 = 0.95
    lvl.3 = 0.85
    lvl.4 = 0.85
    
  }
  
  
  # define corners of first 2 rectangles
  x1  = 10 ; y1  = 15 ; x2  = 15 ; y2  = 20 # compartment I
  x1d = 10 ; y1d =  5 ; x2d = 15 ; y2d = 10 # detritus pool
  
  #    x1  = 10 ; y1  = 25 ; x2  = 15 ; y2  = 30 # compartment I
  #    x1d = 10 ; y1d = 15 ; x2d = 15 ; y2d = 20 # detritus pool
  
  x1c = 10 ; y1c = 5 ;  x2c = 15 ; y2c = 10 # new lindeman I+D
  
  #draw an empty plot
  opar <- par(oma = c(1,1,1,1), mar = c(0,0,0,0))
  plot(c(0,(ntl+1)*10),c(0,30),type='n',xlab='',ylab='',axes=FALSE)
  
  if(type == 1){
    
    # Detritus rectangle and its label
    rect(x1d,y1d,x2d,y2d,lwd=3, border = "grey75")
    text(x1d+2.4,y2d-2,labels='D',cex=lvl.1)
    arrows(x1d-5,y1d+2.5,x1d,y1d+2.5,lwd=2,length = 0.2,angle=0)#arrow Detrital input
    polygon(c(x1d,x1d-1,x1d-1),c(y1d+2.5,y1d+2.9,y1d+2.1),col='black',lwd=2)
    text(x1d-3,y1d+3.5,labels=as.factor(round(ns$DetritalInput, 2)),cex=lvl.3)#value Detrital Input
    arrows(x2d,y2d,x2+5,y2d+5,lwd=2,angle=0,length = 0.2)# arrow for Detritivory from D to II
    polygon(c(x2d+5,x2d+4,x2d+4.5),c(y2d+5,y2d+4.5,y2d+4),col='black',lwd=2)
    text(x2d+2.5,y1-3,labels=as.factor(round(ns$DetritalInput, 2)),cex=lvl.3,pos=4)# value Detritivory
    
    
    #draw all other rectangles
    for(comp in 1:ntl) {
      rect(x1,y1,x2,y2,lwd=3, border = "grey75")#draw rectables
      text(x1+2.4,y2-2, labels = as.roman(comp), cex=lvl.1)#label compartments
      
      #add efficiency to all boxes but last compartment
      if(comp <= ntl-1){
        if(comp == 2){
          ef = (enatroagg$GC[comp+1]*100)/sum(enatroagg$GC[comp]+ns$Detritivory)#eff trophic level 2
        }else{
          ef = enatroagg$GC[comp+1]*100/enatroagg$GC[comp]
        } #efficiency other compartments
        if(ef >= 0.01){
          text(x1+2.4,y2-3.5,labels=paste(round(ef,2),'%',sep=''), cex = lvl.4)
        }else{
          text(x1+2.4,y2-3.5,labels=paste(sprintf('%1.1e',ef),'%',sep=''), cex = lvl.4)
        }
      }
      
      
      #exogenous input
      if(comp == 1){
        arrows(x1d,y1d+2.5,x1d-5,y1d+2.5,lwd=2,length = 0.1,angle=140)
        arrows(x1,y1+2.5,x1-5,y1+2.5,lwd=2,length = 0.1,angle=140)
      }
      
      #Import
      if(!is.null(enatroagg$CI)){
        if(comp > 1){
          ci = enatroagg$CI[comp]
          if(ci > 0){
            arrows(x1-2,y2+4,x1,y2,lwd=2,length = 0.2,angle=0)
            arrows(x1,y2,x1-2,y2+4,lwd=2,length = 0.1,angle=140)
            polygon(c(x1,x1-0.75,x1),c(y2,y2+0.5,y2+1),col='black',lwd=2)
            if(ci >= 0.01){text(x1-1.5,y2+4.9,labels=round(enatroagg$CI[comp],2),cex = lvl.3,pos=2)
            }else{text(x1-1.5,y2+4.9,labels=sprintf('%1.1e',ci),pos=2,cex = lvl.3)}
          }}}
      
      
      #Grazing chain
      arrows(x1-5,y1+2.5,x1,y1+2.5,lwd=2,length = 0.2,angle=0)#draw arrows flow between rectangles
      polygon(c(x1,x1-1,x1-1),c(y1+2.5,y1+2.8,y1+2.2),col='black',lwd=2)
      gc = enatroagg$GC[comp]#back to detrital pool
      if(gc >= 0.01){
        text(x1-3,y1+3.5,labels=format(round(gc,2),nsmall=2),cex=lvl.4)
      }else{
        text(x1-3,y1+3.5,labels=sprintf('%1.1e',gc),cex=lvl.4)
      }
      
      #Export
      arrows(x2,y2,x2+1.5,y2+1.5,lwd=2,angle=0,length = 0.2)# oblique arrow for exports
      polygon(c(x2+1.5,x2+0.5,x2+1),c(y2+1.5,y2+1,y2+0.5),col='white',lwd=2)
      ex = enatroagg$CE[comp]  #back to detrital pool
      if(ex >= 0.01){
        text(x2+2.6,y2+2.2,labels=format(round(ex,2),nsmall=2),cex = lvl.3,pos=2)
      } else{
        text(x2+2.6,y2+2.2,labels=sprintf('%1.1e',ex),pos=2,cex = lvl.3)
      }
      
      #Detrital Pool
      if( comp ==1){
        arrows(x1+2.5,y1,x1+2.5,y1-5,lwd=2,angle=0)
        polygon(c(x1+2.5,x1+2,x1+3),c(y1-5,y1-4,y1-4),col='black',lwd=2)
      } #arrow from I to D
      else{
        arrows(x1+2.5,y1,x1+2.5,y1-7.5,lwd=2,length = 0.2,angle=0)
        polygon(c(x1+2.5,x1+2,x1+3),c(y1-7.5,y1-6.5,y1-6.5),col='black',lwd=2)
      }#other arrows back to D
      dp = enatroagg$RDP[comp]#back to detrital pool
      if(comp == 1){
        if(dp >= 0.01){
          text(x1+2.55,y1-3,labels=format(round(dp,2),nsmall=2), pos=2, cex = lvl.3)
        } else {
          text(x1+2.55,y1-3,labels=sprintf('%1.1e',dp), pos=4, cex = lvl.3)
        }
      } else {
        if(dp >= 0.01){
          text(x1+2.55,y1-3,labels=format(round(dp,2),nsmall=2), pos=4, cex = lvl.3)
        }else{
          text(x1+2.55,y1-3,labels=sprintf('%1.1e',dp), pos=4, cex = lvl.3)
        }
      }
      
      # Respiration
      arrows(x1+3.5,y1,x1+3.5,y1-1.2,lwd=2,angle=90,length=0.1)
      arrows(x1+3.5,y1,x1+3.5,y1-1.5,lwd=2,angle=90,length=0.05)
      re = enatroagg$CR[comp]#respiration
      if(re >= 0.01){text(x1+3.8,y1-1,labels=format(round(re,2),nsmall=2),pos=4,cex = lvl.3)
      }else{text(x1+3.8,y1-1,labels=sprintf('%1.1e',re),pos=4,cex = lvl.3)}
      
      #move coordinates for next rectangle
      x1=x1+10 ; x2=x2+10
    }
    
    arrows(x1-7.5,y1d+2.5,x2d,y1d+2.5,lwd=2,angle=0)#arrow connecting horizontally to Detrital pool
    polygon(c(x2d,x2d+1,x2d+1),c(y1d+2.5,y1d+2.9,y1d+2.1),col='black',lwd=2)
    text(x2d+1.5,y1d+3.5,lab=format(round(sum(enatroagg$RDP[2:ntl]),2),nsmall=2),pos=4,cex = lvl.3)#sum of II-... back to detritus pool
    
  } else {
    
    # UNCHANGED (Oct. 4, 2017)
    
    #################################################################################################
    # add LS merged I+D
    #draw all other rectangles
    ntl = length(enatroagg$GC[enatroagg$LS >0])
    
    for(comp in 1:ntl) {
      rect(x1c, y1c, x2c, y2c, lwd=3) #draw rectables
      if(comp==1){
        text(x1c+2.4,y2c-2,labels='I + D',cex=1.3)
      }else{
        text(x1c+2.4,y2c-2,labels=as.roman(comp),cex=1.3)
      }#label compartments
      #add efficiency to all boxes but last compartment
      if(comp <= ntl-1){
        ef = (enatroagg$TE[comp]*100)#efficiency
        if(ef >= 0.01){text(x1c+2.4,y2c-3.5,labels=paste(round(ef,2),'%',sep=''))
        }else{text(x1c+2.4,y2c-3.5,labels=paste(sprintf('%1.1e',ef),'%',sep=''))}
      }
      
      #exogenous input 'haken'
      if(comp == 1){
        arrows(x1c,y1c+2.5,x1c-5,y1c+2.5,lwd=2,length = 0.1,angle=140)
      }
      
      #Import
      if(!is.null(enatroagg$CI)){
        if(comp > 1){
          ci = enatroagg$CI[comp]
          if(ci > 0){
            arrows(x1c-2,y2c+4,x1c,y2c,lwd=2,length = 0.2,angle=0)
            arrows(x1c,y2c,x1c-2,y2c+4,lwd=2,length = 0.1,angle=140)
            polygon(c(x1c,x1c-0.75,x1c),c(y2c,y2c+0.5,y2c+1),col='black',lwd=2)
            if(ci >= 0.01){text(x1c-1.5,y2c+4.9,labels=round(enatroagg$CI[comp],2),cex = lvl.3,pos=2)
            }else{text(x1c-1.5,y2c+4.9,labels=sprintf('%1.1e',ci),cex = lvl.3,pos=2)}
          }}}
      
      
      
      #Grazing chain
      arrows(x1c-5,y1c+2.5,x1c,y1c+2.5,lwd=2,length = 0.2,angle=0)#draw arrows flow between rectangles
      polygon(c(x1c,x1c-1,x1c-1),c(y1c+2.5,y1c+2.9,y1c+2.1),col='black',lwd=2)
      gc = enatroagg$LS[comp]# GC value
      
      if(comp==1){
        text(x1c-5,y1c+4,labels=format(round(gc,2),nsmall=2),font=2,cex=1.0,pos=2)#total input to I+D, bold
        text(x1c-3,y1c+3.5,labels=round(ns$DetritalInput + enatroagg$GC[comp],2),cex=0.75)
      }else{
        if(gc >= 0.01){
          text(x1c-3,y1c+3.5,labels=format(round(gc,2),nsmall=2),cex=0.75)
        } else {
          text(x1c-3,y1c+3.5,labels=sprintf('%1.1e',gc),cex=0.75)
        }
      }
      #Export
      arrows(x2c,y2c,x2c+1.5,y2c+1.5,lwd=2,angle=0,length = 0.2)# oblique arrow for exports
      polygon(c(x2c+1.5,x2c+0.5,x2c+1),c(y2c+1.5,y2c+1,y2c+0.5),col='white',lwd=2)
      ex = enatroagg$CE[comp]#back to detrital pool
      
      if(comp==1){
        text(x2c+2.6,y2c+2.2,
             labels=round(sum(u$output[c(primprod,nonliving)], na.rm = T), 2),
             cex = lvl.3, pos=2)
        arrows(x1c+2.5,y1c,x1c+2.5,y1c-2,lwd=2,angle=90,length=0.0)
        arrows(x1c+2.5,y1c-2,x1c+0.5,y1c-2,lwd=2,angle=90,length=0.0)
        arrows(x1c+0.5,y1c-2,x1c+0.5,y1c,lwd=2,angle=90,length=0.0)
        polygon(c(x1c+0.5,x1c,x1c+1),c(y1c,y1c-1,y1c-1),col='black',lwd=2)
        
      }else{
        if(ex >= 0.01){
          text(x2c+2.6,y2c+2.2,labels=format(round(ex,2),nsmall=2),cex = lvl.3,pos=2)
        }else{
          text(x2c+2.6,y2c+2.2,labels=sprintf('%1.1e',ex),pos=2,cex = lvl.3)
        }
      }
      #Detrital Pool
      dp = enatroagg$RDP[comp]#back to detrital pool
      if( comp ==1){arrows(5,y1c+1.7,x1c,y1c+1.7,lwd=2,angle=0)
        polygon(c(x1c,x1c-1,x1c-1),c(y1c+1.7,y1c+2,y1c+1.4),col='black',lwd=2)
        text(x1c,y1c-3,labels=format(round(dp,2),nsmall=2),pos=4,cex = lvl.3)
      }else{arrows(x1c+2.5,y1c,x1c+2.5,y1c-5,lwd=2,length = 0.2,angle=0)
        polygon(c(x1c+2.5,x1c+2,x1c+3),c(y1c-5,y1c-4,y1c-4),col='black',lwd=2)#other arrows back to D
        if(dp >= 0.01){text(x1c+2.55,y1c-3,labels=format(round(dp,2),nsmall=2),pos=4,cex = lvl.3)
        }else{text(x1c+2.55,y1c-3,labels=sprintf('%1.1e',dp),pos=4,cex = lvl.3)}}
      
      # Respiration
      arrows(x1c+3.5,y1c,x1c+3.5,y1c-1.2,lwd=2,angle=90,length=0.1)
      arrows(x1c+3.5,y1c,x1c+3.5,y1c-1.5,lwd=2,angle=90,length=0.05)
      re = enatroagg$CR[comp]#respiration
      if(re >= 0.01){text(x1c+3.8,y1c-1,labels=format(round(re,2),nsmall=2),pos=4,cex = lvl.3)
      }else{text(x1c+3.8,y1c-1,labels=sprintf('%1.1e',re),pos=4,cex = lvl.3)}
      
      x1c=x1c+10 ; x2c=x2c+10
    }
    arrows(x1c-7.5,y1c-5,5,y1c-5,lwd=2,angle=0)#arrow connecting horizontally to Detrital pool
    arrows(5,y1c-5,5,y1c+1.7,lwd=2,angle=0)#arrow connecting horizontally to Detrital pool
    #polygon(c(x2d,x2d+1,x2d+1),c(y1d+1.5,y1d+2.9,y1d+2.1),col='black',lwd=2)
    text(5,y1c+0.8,lab=format(round(sum(enatroagg$RDP[2:ntl]),2),nsmall=2),pos=4,cex = lvl.3)#sum of II-... back to detritus pool
  }
  rm(opar)
}

#Trophic Spine Visualization
primprod <- c(3:6)
plot.lindeman(origram, enatroagg = "troagg", primprod = primprod, type = 2)
plot.lindeman(invgram, enatroagg = "troagg", primprod = primprod, type = 2)

#Pre invaded Trophic proportions by taxonomy
orig_troagg <- enaTroAgg(origram)

orig_TL_prop <- matrix(data = NA, nrow = 63, ncol = 63)
for (i in 1:length(orig_troagg$A[1,])) {
  for (j in 1:length(orig_troagg$A[1,])) {
    orig_TL_prop[i,j] <- orig_troagg$A[i,j]/sum(orig_troagg$A[i,])
  }
}

colnames(orig_TL_prop) <- colnames(orig_troagg$A)
orig_TL_prop <- orig_TL_prop[1:6,]

orig_TL_prop_df <- as.data.frame(t(orig_TL_prop))
colnames(orig_TL_prop_df) <- c("TL1", "TL2", "TL3", "TL4", "TL5", "TL6")
# tax.num.delin <- as.factor(c(2,2,1,1,1,1,3,3,3,3,3,3,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,7,7,7,7,
#                    10,10,11,11,11,12,12,9,12,12,8,8,8,8,8,8,8,8))
tax.num.delin <- as.factor(c(2,2,1,1,1,1,3,3,3,3,3,3,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,7,7,7,7,
                             11,11,10,12,12,14,14,10,14,13,9,9,9,9,9,9,9,9))
orig_TL_prop_df$tax.num.delin <- tax.num.delin
orig_TL_prop_df <- orig_TL_prop_df %>% replace(is.na(.), 0)
orig_TL_prop_by_taxonomy <- aggregate(. ~ tax.num.delin, orig_TL_prop_df, sum)
orig_TL_prop_by_taxonomy$taxonomy <- c("Primary Producer","Microfauna","Aquatic Invertebrate","Terrestrial Invertebrate",
                                       "Fish","Amphibean","Reptile","Bird","Mammalian Herbivore","Rodent","Meso-Mammal", "Panther", "Other Mammalian Carnivores")

orig_TL_prop_by_taxonomy <- orig_TL_prop_by_taxonomy %>% 
  tidyr::pivot_longer(starts_with("TL"), names_to = "Trophic.Level", values_to = "Proportion") %>% 
  arrange(tax.num.delin)

orig_TL_prop_by_taxonomy$Trophic.Level <- gsub("TL", "", orig_TL_prop_by_taxonomy$Trophic.Level)

#lock in data frame order
orig_TL_prop_by_taxonomy$Trophic.Level <- factor(orig_TL_prop_by_taxonomy$Trophic.Level, 
                                                 levels = unique(orig_TL_prop_by_taxonomy$Trophic.Level))
#color vector
mycol <- c("green4","greenyellow","tan3","tan4","yellow","orange","orangered","maroon","lightblue","skyblue2","skyblue3", "blue", "skyblue4")
label.order <- c("Primary Producer","Microfauna","Aquatic Invertebrate","Terrestrial Invertebrate",
                 "Fish","Amphibean","Reptile","Bird","Mammalian Herbivore","Rodent","Meso-Mammal", "Panther","Other Mammalian Carnivores")

#graph
pre_invaded_TL <- ggplot(orig_TL_prop_by_taxonomy, aes(x = Trophic.Level, y = Proportion))+
  geom_bar(aes(fill = tax.num.delin),position = "stack", stat = "identity")+
  xlab("Trophic Level")+
  labs(fill = "Taxonomy")+
  ggtitle("Pre-Invaded")+
  scale_fill_manual(values = mycol, labels = label.order)+
  theme(legend.title = element_text(18), legend.text = element_text(15))

#post invaded Trophic proportions by taxonomy
inv_troagg <- enaTroAgg(invgram)

inv_TL_prop <- matrix(data = NA, nrow = 64, ncol = 64)
for (i in 1:length(inv_troagg$A[1,])) {
  for (j in 1:length(inv_troagg$A[1,])) {
    inv_TL_prop[i,j] <- inv_troagg$A[i,j]/sum(inv_troagg$A[i,])
  }
}
colnames(inv_TL_prop) <- colnames(inv_troagg$A)
inv_TL_prop <- inv_TL_prop[1:6,]

inv_TL_prop_df <- as.data.frame(t(inv_TL_prop))
colnames(inv_TL_prop_df) <- c("TL1", "TL2", "TL3", "TL4", "TL5", "TL6")

tax.num.delin <- as.factor(c(8,2,2,1,1,1,1,3,3,3,3,3,3,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,6,6,6,6,6,6,7,7,7,7,
                             11,11,10,12,12,14,14,10,14,13,9,9,9,9,9,9,9,9))
inv_TL_prop_df$tax.num.delin <- tax.num.delin
inv_TL_prop_by_taxonomy <- aggregate(. ~ tax.num.delin, inv_TL_prop_df, sum)
inv_TL_prop_by_taxonomy$taxonomy <- c("Primary Producer","Microfauna","Aquatic Invertebrate","Terrestrial Invertebrate",
                                      "Fish","Amphibian","Reptile", "Python", "Bird","Mammalian Herbivore","Rodent","Meso-Mammal", "Panter","Other Mammalian Carnivores")
inv_TL_prop_by_taxonomy <- inv_TL_prop_by_taxonomy %>% 
  tidyr::pivot_longer(starts_with("TL"), names_to = "Trophic.Level", values_to = "Proportion") %>% 
  arrange(tax.num.delin)

inv_TL_prop_by_taxonomy$Trophic.Level <- gsub("TL", "", inv_TL_prop_by_taxonomy$Trophic.Level)
#Visualization
#lock in data frame order
inv_TL_prop_by_taxonomy$Trophic.Level <- factor(inv_TL_prop_by_taxonomy$Trophic.Level, 
                                                levels = unique(inv_TL_prop_by_taxonomy$Trophic.Level))

#color vector
mycol <- c("green4","greenyellow","tan3","tan4","yellow","orange","orangered", "black", "maroon","lightblue","skyblue2","skyblue3","blue","skyblue4")
label.order <- c("Primary Producer","Microfauna","Aquatic Invertebrate","Terrestrial Invertebrate",
                 "Fish","Amphibian","Reptile", "Python","Bird","Mammalian Herbivore","Rodent","Meso-Mammal", "Panther", "Other Mammalian Carnivores")

#graph
post_invaded_TL <- ggplot(inv_TL_prop_by_taxonomy, aes(x = Trophic.Level, y = Proportion, fill = tax.num.delin))+
  geom_bar(position = "stack", stat = "identity")+
  xlab("Trophic Level")+
  ggtitle("Post-Invaded")+
  scale_fill_manual(values = mycol, labels = label.order)+
  labs(fill = "Taxonomy")+
  theme(legend.title = element_text(18), legend.text = element_text(15))


gridExtra::grid.arrange(pre_invaded_TL,post_invaded_TL, nrow = 1, ncol =2)

#Compositional Shift
pyth.prop <- inv_TL_prop_by_taxonomy[which(inv_TL_prop_by_taxonomy$tax.num.delin == 8),]
pyth.prop$prop.diff <- pyth.prop$Proportion
test.con <- inv_TL_prop_by_taxonomy %>% filter(tax.num.delin != 8)
test.con$prop.diff <- test.con$Proportion - orig_TL_prop_by_taxonomy$Proportion
test.con <- rbind(test.con, pyth.prop)
test.con <- test.con[order(test.con$tax.num.delin, decreasing = FALSE), ]
tmycol <- c("green4","greenyellow","tan3","tan4","yellow","orange","orangered", "black", "maroon","lightblue","skyblue2","skyblue3","blue","skyblue4")
tlabel.order <- c("Primary Producer","Microfauna","Aquatic Invertebrate","Terrestrial Invertebrate",
                  "Fish","Amphibean","Reptile", "Python","Bird","Mammalian Herbivore","Rodent","Meso-Mammal", "Panther", "Other Mammalian Carnivores")

#graph
ggplot(test.con, aes(x = Trophic.Level, y = prop.diff, fill = tax.num.delin))+
  geom_bar(position = "stack", stat = "identity")+
  scale_fill_manual(values = tmycol, labels = tlabel.order)+
  xlab("Trophic Level")+
  ylab("Difference")+
  ggtitle("Trophic Compostional Shift Post Invasion")+
  labs(fill = "Taxonomy")+
  theme(legend.title = element_text(20), legend.text = element_text(17))+
  geom_hline(yintercept = 0, size = 0.8)+
  geom_segment(aes(x = 1, y = 0.05, xend = 1, yend = 0.35),
               arrow = arrow(length = unit(0.5, "cm")), size = 1.5)+
  geom_segment(aes(x = 1, y = -0.05, xend = 1, yend = -0.35),
               arrow = arrow(length = unit(0.5, "cm")), size = 1.5)+
  annotate("text", x=2, y=-0.4, label= "decreasing composition") +
  annotate("text", x=2, y=0.4, label= "increasing composition")

#total proportional shift of biomass composition post invasion across trophic levels
test.con %>% filter(prop.diff > 0) %>% group_by(Trophic.Level) %>% summarize(shift = sum(prop.diff))
#table of proportional differences different than zero
df <- test.con[which(test.con$prop.diff != 0),]
df <- df[order(df$Trophic.Level),]
#save table as csv
write.csv(df, file = "tl.comp.shift.csv")

#Mixed Trophic Impacts####
invaded_mti <- enaMTI(invgram)
longData <- melt(invaded_mti$M)
ggplot(data = longData, mapping = aes(x = Var2, y = Var1))+
  geom_raster(aes(fill=value))+
  scale_fill_gradient2(low = "red", mid = "white", high = " dark green", limits = c(-1,1))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Impact on")+
  ylab("Impact from")

python_mti <- as.data.frame(invaded_mti$M[1,])
colnames(python_mti) <- "value"
python_mti$name <- rownames(python_mti)

python_mti_peffect <- python_mti %>% filter(value >= 2e-03); python_mti_peffect
python_mti_neffect <- python_mti %>% filter(value <= 0); python_mti_neffect
python_mti_effect <- rbind(python_mti_neffect, python_mti_peffect)
python_mti_effect <- python_mti_effect[-c(2:3),]

ggplot(data = python_mti_effect, mapping = aes(x = reorder(name, value), y = value))+
  geom_bar(stat = "identity", aes(fill = ifelse(value<0, "red", "blue")))+
  theme(axis.text.x = element_text(angle =100, hjust = 1, vjust = 0.1, size = 15))+
  theme(axis.text.y = element_text(size = 15))+
  ylab("Trophic Impact")+
  xlab("Species")+
  theme(axis.ticks = element_blank(), 
        legend.position = "none")

#Biomass comparison pre and post invasion
biomass_diff <- as.data.frame((post_biom[-1] - biom)/biom)
colnames(biomass_diff) <- "Change.in.Biomass"
biomass_diff$name <- rownames(Flow)
biomass_diff <- biomass_diff[which(biomass_diff$Change.in.Biomass != 0),]

ggplot(data = biomass_diff, mapping = aes(x = reorder(name, Change.in.Biomass), y = Change.in.Biomass))+
  geom_bar(stat = "identity", aes(fill = ifelse(Change.in.Biomass<0, "red", "blue")))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = -0.1, size = 7))+
  xlab("Node")+
  ylab("Proportional change in Biomass")+
  theme(axis.text.y = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 15))+
  theme(axis.text=element_text(size=20))+
  theme(axis.ticks = element_blank(), 
        legend.position = "none")

#Network Indicators####

#Functional Properties
origFlow <- enaFlow(origram)
invFlow <- enaFlow(invgram)
production_to_biomass <- c(origFlow$ns[,3]/sum(biom[1:63]), invFlow$ns[,3]/sum(post_biom[1:64]))
origAsc <- enaAscendency(origram)
invAsc <- enaAscendency(invgram)

#Table 
Flow_indices <- as.data.frame(rbind(origFlow$ns[,3:5],invFlow$ns[,3:5]))
Flow_indices <- t(Flow_indices)
colnames(Flow_indices) <- c("Pre.Burmese.Python", "Post.Burmese.Python")
Flow_indices <- as.data.frame(rbind(Flow_indices, production_to_biomass))

#Table
asc_indices <- as.data.frame(rbind(origAsc[,c(2,4:8)], invAsc[,c(2,4:8)]))
asc_indices <- t(asc_indices)
colnames(asc_indices) <- c("Pre.Burmese.Python", "Post.Burmese.Python")

#Ecosystem Indices Table
ecosystem_indicies <- as.data.frame(rbind(Flow_indices, asc_indices))
ecosystem_indicies <- ecosystem_indicies %>% 
  mutate(Percent.Difference = ((ecosystem_indicies[,2] - ecosystem_indicies[,1])/(ecosystem_indicies[,1]))*100)
rownames(ecosystem_indicies) <- c("Total System Troughput", "Average Path Length", "Finn Cycling Index",
                                  "Production:Biomass", "Average Mutual Information","Capacity", "Ascendancy", "Overhead", 
                                  "Relative Ascendancy", "Relative Overhead")
write.csv(ecosystem_indicies, "./Results/ecosystem_indicies.csv")

#Structural properties
pre_community <- LoadCommunity('C:/Users/biogirl92/OneDrive - University of Florida/Documents/Graduate Assistanship/Thesis Project/Code/github/data_files/cheddar_format_pre')

resource <- rep(rownames(post_Flow), each = length(post_biom))
consumer <- rep(colnames(post_Flow), times = length(post_biom))
weight <- vector(length = length(post_biom))
for (i in 1:length(post_biom)) {
  if(i == 1){
    weight <- post_Flow[i,]
  } else{
    weight <- c(weight, post_Flow[i,])
  }
}
trophic.links <- data.frame(cbind(resource, consumer, weight))
rownames(trophic.links) <- NULL
trophic.links$weight <- as.numeric(trophic.links$weight)
trophic.links <- trophic.links %>% filter(trophic.links[,3] != 0)
trophic.links <- trophic.links[-which(trophic.links[,2] == "Utricularia" | trophic.links[,2] == "Sediment Carbon" 
                                      | trophic.links[,2] == "Labile Detritus" | trophic.links[,2] == "Refractory Detritus"),]

post_community <- Community(nodes = nodes, properties = properties, trophic.links = trophic.links)

origdesc <- QuantitativeRevised(pre_community, 'weight', top.level.threshold=0.99)
invdesc <- QuantitativeRevised(post_community, 'weight', top.level.threshold=0.99)

topology_indices <- cbind(origdesc[,3], invdesc[,3], (((invdesc[,3] - origdesc[,3])/origdesc[,3])*100))
colnames(topology_indices) <- c("Pre.Burmese.Python", "Post.Burmese.Python", "Percent.Difference")
write.csv(topology_indices, file = "./Results/topology_indices.csv")

###Data Visulaization for system level metrics
sys_tbl <- rbind(ecosystem_indicies, topology_indices)
sys_tbl <- sys_tbl[-c(17:20,23:24),]
sys_tbl <- round(sys_tbl, 2)
sys_tbl <- cbind(row.names(sys_tbl), sys_tbl[,1:3])
colnames(sys_tbl) <- c("Metric", "Pre.Invasion", "Post.Invasion", "Percent.Difference")

net_table <- gt(sys_tbl)
net_table <- net_table |>
  tab_header(
    title = md("**Pre- and Post-Invasion**"),
    subtitle = md("Functional Metrics")
  ) |>
  tab_row_group(
    label = md("*Structural*"),
    rows = c(11:18)
  ) |>
  tab_row_group(
    label = md("*Functional*"),
    rows = 1:10
  ) |>
  tab_style(
    style = cell_text(weight = "bold"), 
    locations = cells_column_labels(columns=c("Metric", "Pre.Invasion", "Post.Invasion", "Percent.Difference"))
  )

#Post_Invasion PCoA####

invCentral <- centrality(invgram)

ETL <- inv_troagg$ETL
TDC <- apply(invFlow$TDC, 1, sum)
ASC <- asc_coef(invgram)
C.Centrality <- invCentral$C.Centrality
T.Centrality <- invCentral$T.Centrality

#compartment by metric matrix
inv_functional <- cbind(ETL, TDC, ASC, T.Centrality, C.Centrality)

#Data transformation and scaling
par(mfrow = c(2,3))
mapply(hist, as.data.frame(inv_functional), main = colnames(inv_functional))

log.inv <- inv_functional
log.inv[, c(2,4:5)] <- log10(inv_functional[,c(2,4:5)])

par(mfrow = c(2,3))
mapply(hist, as.data.frame(log.inv), main = colnames(log.inv))

#test for multivariate normality - Hard NO
mshapiro.test(t(inv_functional))
mvn(inv_functional, mvnTest = "mardia")
#standardization
rsum <- rowSums(log.inv)
csum <- colSums(log.inv)

cv(csum)
cv(rsum)
Zinv <- scale(log.inv); Zinv

##PCA
inv_functional_pca <- princomp(Zinv, cor = F)
#eigen values
eigenVal<- (inv_functional_pca$sdev*sqrt(41/40))^2
propVar<-eigenVal/sum(eigenVal)
cumVar<-cumsum(propVar)
pca_Table<-t(rbind(eigenVal,propVar,cumVar))
pca_Table

##cluster analysis
#determinging number of groups
wss <- rep(0, 20)

#Run a loop for 1 to 20 clusters:
for (i in 1:20) # sets the number of times the loop will be run i.e., the number of clusters in this case)
  
  wss[i] <- sum(kmeans(Zinv, centers = i,nstart=25)$withinss) # run the kmeans function for each number of clusters (i) and extract the within sum of squares for each.

#Vector of within group sum of squares
wss 

par(mfrow = c(1,1))
plot(1:20, wss, type = "b", xlab = "Number of groups", ylab = "Within groups sum of squares") 

#silouette widths
sil <- rep(0,20)
for (i in 2:20)
  sil[i] <- summary(silhouette(kmeans(Zinv, centers=i, iter.max=100, nstart=25)$cluster, dist(Zinv)))$avg.width
plot(2:20, sil[2:20], type = "b", xlab = "Number of groups", ylab = "average silhouette width ")

#4 clusters
inv.kop <- kmeans(Zinv, centers= 6, iter.max=10, nstart=25)
pairs(Zinv, panel=function(x,y,z) text(x,y,inv.kop$cluster))
#testing groups
groups <- inv.kop$cluster
set.seed(11)
inv_MRPP<-mrpp(Zinv, groups, permutations = 1000) ##From Here

hist(inv_MRPP$boot.deltas, xlim = c(1.0,3.2), main = "Histogram of pre-functional MRPP deltas" )
points(inv_MRPP$delta,0, pch=19,col="red", bg="red", cex=2)

loadings(inv_functional_pca)
group.num.delin <- as.factor(c(5,1,2,1,1,2,1,2,2,2,2,2,3,2,5,5,3,3,3,2,3,2,2,2,2,3,5,3,3,5,3,3,5,5,5,4,3,4,
                               5,5,5,4,3,5,3,5,2,3,4,4,3,5,6,3,6,6,4,6,4,6,6,4,4,3,1,1,1))

fviz_pca_biplot(inv_functional_pca, label ="ind",geom.ind="point", pointshape = 19, habillage = group.num.delin, title = "Post-Invaded PCA",
                xlim = c(-5,5), ylim = c(-5,5),addEllipses=TRUE, ellipse.level=0.95)+
  xlab("PC1 (72.8%)")+
  ylab("PC2 (19.1%)")+
  geom_label(label = "periphyton", x = -4.09, y = 0.27, label.padding = unit(0.2, "lines"), color = "black")+
  geom_label(label = "refractory detritus", x = -3.91, y = -0.15, label.padding = unit(0.2, "lines"), color = "black")+
  geom_label(label = "rats&mice", x = -0.7, y = -0.05, label.padding = unit(0.2, "lines"), color = "black")+
  geom_label(label = "turtles", x = -0.5, y = -0.89, label.padding = unit(0.2, "lines"), color = "black")+
  geom_label(label = "dollar sunfish", x = -0.05, y = 0.69, label.padding = unit(0.2, "lines"), color = "black")+
  geom_label(label = "mesoinverts", x = -2.01, y = 0.1, label.padding = unit(0.2, "lines"), color = "black")+
  geom_label(label = "killifishes", x = -1.03, y = -0.57, label.padding = unit(0.2, "lines"), color = "black")+
  geom_label(label = "python", x = 0.69, y = 1.86, label.padding = unit(0.2, "lines"), color = "red")+
  geom_label(label = "bass", x = 0.72, y = 1.48, label.padding = unit(0.2, "lines"), color = "black")+
  geom_label(label = "alligator", x = 1.24, y = 1.3, label.padding = unit(0.2, "lines"), color = "black")+
  geom_label(label = "panther", x = 3.62, y = -0.09, label.padding = unit(0.2, "lines"), color = "black")+
  geom_label(label = "bobcat", x = 2.82, y = 0.46, label.padding = unit(0.2, "lines"), color = "black")+
  geom_label(label = "Gruiformes", x = 1.27, y = -1.49, label.padding = unit(0.2, "lines"), color = "black")+
  geom_label(label = "tadpoles", x = 1.53, y = -2.26, label.padding = unit(0.2, "lines"), color = "black")+
  geom_label(label = "raccoons", x = 1.15, y = -0.81, label.padding = unit(0.2, "lines"), color = "black")
#Trophic Similarity
sp.similarity <- TrophicSimilarity(post_community)
pyth.similarity <- sp.similarity[,1]

