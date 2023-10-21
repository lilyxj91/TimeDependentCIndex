##################################################################
##  Name:           TimeDependentCIndex.R  
##
##  Description:    The original code developed for evaluating dynamic 
##                  and predictive discrimination for recurrent event models
##                  using a time-dependent C-index. 
##
## Functions:       TimeDependentCIndex
##                  my.reReg
##                  objective_ee
##                  objective_pl
##                  My.counts
##                  alpha.estimate
##
## Author:          Jian Wang
## Completion date: May 2023                                          
## Modified:         
## Date             Initials                  Change                           
##################################################################


##################################################################
## Name :              TimeDependentCIndex
##
## Description:        Time-dependent C-index for recurrent events
##
## Usage:              TimeDependentCIndex(data,outcome.type,casecontrol,
##                                      prev,printout)
##
## Required arguments: train.data:    Input training dataset file in the .csv format.
##                              
##                     validate.data: Input validatation dataset file in the .csv format.
##                               
##                     n.var:         Number of predictors. 
##                                  
##                     B:             Number of perburtations.
##
## Author:             Jian Wang                        
## Completion date:    May 2023                            
## Modified:                                       
## Date                Initials                   Change          
##################################################################

TimeDependentCIndex <- function(train.data,
                                validate.data,
                                n.var,
                                B=200)
{
  
  ## Use training data (train.data) to calculate the recurrent event regression coefficients
  ## For estimation
  n<-length(unique(train.data$id)) ## sample size
  init.value<-matrix(c(0.5,0.5,0.5),ncol=1)
  fit<-my.reReg(train.data,init.value,X=paste('X',1:n.var,sep = ""))
  ## For variance estimation
  ## perturbation estimation
  fit.pb.tmp<-data.frame(matrix(NA,B,n.var))
  set.seed(1)  ## for reproducing the results
  for (ii in 1:B){
    w.list<-rep(rexp(n),table(train.data$id))
    fit.pb.tmp[ii,]<-my.reReg.pb(w.list,train.data,init.value,X=paste('X',1:n.var,sep = ""))
  }
  results_train<-fit
      
  ## Use validation data (validate.data) to calculate the C-index
  data<-validate.data
  ## Calculate the risk prediction scores
  ## For estimation
  data$S_z<-t(as.numeric(results_train) %*% t(data[,paste('X',1:n.var,sep = "")]))
  ## For variance estimation
  pb.list<-list()
  for (ii in 1:B){
    tmp.list<-list()
    tmp.list$S_z<-t(as.numeric(fit.pb.tmp[ii,]) %*% t(data[,paste('X',1:n.var,sep = "")]))
    pb.list[[ii]]<-tmp.list
  }
    
  ## Calculate the alpha parameters
  s_h<-sort(unique(data$t))
  max.fu<-max(validate.data$t.obs)
  ## The parameter estimation
  data.risk<-data[,-which(names(data) %in% c("indi",paste('X',1:n.var,sep = "")))];
  ## Calculate the counts
  tableij = lapply(s_h[2:length(s_h)],My.counts,data.risk)
  ## alpha.estimate
  results_alpha<-alpha.estimate(tableij,s_h,max.fu)
      
  ## The variance estimation: perturbation
  n.validate<-length(unique(data$id))
  alpha.tmp<-data.frame(matrix(NA,B,3))
  set.seed(2) # for reproducing the results
  for (ii in 1:B){
    data.risk$w<-rep(rexp(n.validate),table(data.risk$id))
    data.risk$S_z<-pb.list[[ii]]$S_z
    tableij = lapply(s_h[2:length(s_h)],My.counts.pb,data.risk)
    alpha.tmp[ii,]<-alpha.estimate(tableij,s_h,max.fu)
  }
  results_alpha.pb.all<-alpha.tmp
  
  ## Create the C-index plot
  
  ## For the plot, because the curve becomes increasingly variable at 
  ## the end of the study due to less information available, we will
  ## not plot the entire study time period. In stead, we plot the curve
  ## up to the 95% quantile of the maximum followup of the data.
  cutoff<-sort(data$t)[round(length(data$t)[1]*0.95)]
  
  t<-seq(0,cutoff,0.1)
  tt<-(t-max.fu/2)/(max.fu/2)
  C.t<-data.frame(matrix(NA,1,length(tt)))
  a0<-results_alpha[1];a1<-results_alpha[2];a2<-results_alpha[3]
  C.t<-exp(a0+a1*tt+a2*tt^2)/(1+exp(a0+a1*tt+a2*tt^2))
  ## Calculate the confidence interval
  C.t.sd<-C.t.low<-C.t.high<-data.frame(matrix(NA,1,length(tt)))
  alpha.tmp<-results_alpha.pb.all
  C.tmp<-data.frame(matrix(NA,dim(alpha.tmp)[1],length(tt)))
  for (kk in 1:dim(alpha.tmp)[1]){
    a0<-alpha.tmp[kk,1];a1<-alpha.tmp[kk,2];a2<-alpha.tmp[kk,3]
    C.tmp[kk,]<-exp(a0+a1*tt+a2*tt^2)/(1+exp(a0+a1*tt+a2*tt^2))
  }
  C.t.sd<-apply(C.tmp,2,sd)  
  C.t.low<-C.t-1.96*C.t.sd
  C.t.high<-C.t+1.96*C.t.sd
  
  data.used<-data.frame(x=t,y<-C.t,
                        y.low<-C.t.low,
                        y.high<-C.t.high)
  names(data.used)<-c('Time','y','y.low','y.high')
  
  maintitle <-paste('Time-dependent C-index')
  
  f<-ggplot(data = data.used, aes(x = Time))+
    geom_line(aes(y=y,colour="y",linetype="y"), size=1.5)+
    geom_line(aes(y=`y.low`,colour="y",linetype="95% CI"), size=0.8)+
    geom_line(aes(y=`y.high`,colour="y",linetype="95% CI"), size=0.8)+
    scale_colour_manual("",
                        breaks = c("y"),
                        values = c("y"="blue")) +
    scale_linetype_manual("",values=c("y"="solid",
                                      "95% CI"="dotdash"))+
    guides(linetype="none")+
    xlab('Time (years)') + ylab(expression(hat(C)~(t))) +
    labs(color="Legend")+
    scale_y_continuous(limits=c(0, 0.8),breaks=seq(0,0.8,0.2))+
    scale_x_continuous(limits=c(0, max(data.used$Time)),breaks=seq(0,15,1))+
    labs(title =maintitle)+ 
    theme_bw() +
    theme(axis.text.x = element_text(colour = "black", size = 23),
          axis.title.x = element_text(colour = "black", size = 27),
          axis.text.y = element_text(colour = "black", size = 23),
          axis.title.y = element_text(colour = "black", size = 27),         
          panel.grid.major =element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size =27, colour = "black", vjust = 1)) +
    theme(plot.margin = unit(c(1,1,1,1), "cm"))+
    theme(legend.position = "none")+
    guides(colour = guide_legend(nrow = 1))
  print(f)

  ## Return the alpha estimations
  return(list(results_alpha,results_alpha.pb.all))
 
}

############################################
## recurrent event regression ##############
############################################
## proportional means model
my.reReg<-function(data,beta1,X){
  data0j=data%>%group_by(id)%>%slice(1)%>%dplyr::select(-c("indi","t"))
  data0i=data%>%ungroup()%>%filter(Num!=0)%>%dplyr::select(-c("indi","t.obs"))
  
  input_list = list()
  for (i in 1:length(data0i$id)){
    this_input=list()
    this_input$id = data0i$id[i]
    this_input$t = data0i$t[i]
    this_input$xi =  as.matrix(data0i[i,X])
    selection=which(data0j$t.obs>=this_input$t)
    this_input$xij =  as.matrix(data0j[selection,X]) 
    input_list[[i]] = this_input
  }
  
  #par1 = neldermead(beta1, fn =objective_ee,input_list=input_list,control =list(xtol_rel = 1e-10, maxeval = 10000))$par
  par2 = neldermead(beta1, fn =objective_pl,input_list=input_list,control =list(xtol_rel = 1e-10, maxeval = 10000))$par
  
  return(par2)
}

## For variance estimation: using perturbation
my.reReg.pb<-function(w,data,beta1,X){
  data$w<-w
  data0j=data%>%group_by(id)%>%slice(1)%>%dplyr::select(-c("indi","t"))
  data0i=data%>%ungroup()%>%filter(Num!=0)%>%dplyr::select(-c("indi","t.obs"))
  
  input_list = list()
  for (i in 1:length(data0i$id)){
    this_input=list()
    this_input$id = data0i$id[i]
    this_input$t = data0i$t[i]
    this_input$xi =  as.matrix(data0i[i,X])
    this_input$wi = data0i$w[i]
    selection=which(data0j$t.obs>=this_input$t)
    this_input$xij =  as.matrix(data0j[selection,X]) #cbind(data0j$X1[selection],data0j$X2[selection],data0j$X3[selection])
    input_list[[i]] = this_input
  }
  
  # par1.pb = neldermead(beta1, fn =objective_ee_pb,input_list=input_list,control =list(xtol_rel = 1e-10, maxeval = 10000))$par
  par2.pb = neldermead(beta1, fn =objective_pl_pb,input_list=input_list,control =list(xtol_rel = 1e-10, maxeval = 10000))$par
  
  return(par2.pb)
}

## Objective functions
## Method A
## If use estimating equation

objective_ee = function(beta,input_list){
  uu = matrix(c(0,0,0),ncol=3) 
  n = length(input_list)
  for (i in 1:n){
    xi = input_list[[i]]$xi
    xij = input_list[[i]]$xij
    g=c(xij%*%beta)
    uu = uu + xi - exp(g)%*%xij/sum(exp(g)) 
  }
  return (sum(abs(uu)^2))
}

objective_ee_pb = function(beta,input_list){
  uu = matrix(c(0,0,0),ncol=3) 
  n = length(input_list)
  for (i in 1:n){
    xi = input_list[[i]]$xi
    xij = input_list[[i]]$xij
    wi = input_list[[i]]$wi
    g=c(xij%*%beta)
    uu = uu + wi*xi - wi*exp(g)%*%xij/sum(exp(g)) 
  }
  return (sum(abs(uu)^2))
}


## Method B
## If mimic log partial likelihood

objective_pl = function(beta,input_list){
  nloglikeli = 0
  n = length(input_list)
  for (i in 1:n){
    xi = input_list[[i]]$xi
    xij = input_list[[i]]$xij
    nloglikeli = nloglikeli - xi%*%beta + log(sum(exp(xij%*%beta))) 
  }
  return (nloglikeli)
}

objective_pl_pb = function(beta,input_list){
  nloglikeli = 0
  n = length(input_list)
  for (i in 1:n){
    xi = input_list[[i]]$xi
    xij = input_list[[i]]$xij
    wi = input_list[[i]]$wi
    nloglikeli = nloglikeli - wi*xi%*%beta + wi*log(sum(exp(xij%*%beta))) 
  }
  return (nloglikeli)
}


######################################################
## calculate the counts of concordance ###############
######################################################
## For estimation
My.counts<-function(ss,data.risk){
  data.risk<-data.risk[data.risk$t.obs>ss,]
  data.risk$s_h<-ss
  data.risk$indicator<-(data.risk$t<=apply(data.risk[,c('s_h','t.obs')],1,min) & data.risk$Num!=0)
  data.risk.agg<-data.risk %>% group_by(id) %>% mutate(N.counts = sum(indicator))
  data.risk.agg<-aggregate(data.risk.agg[,-which(names(data.risk.agg) %in% c('t','indicator','Num'))], by=list(ID=data.risk.agg$id), FUN=unique)
  
  data.risk.agg<-data.risk.agg[order(data.risk.agg$N.counts),]
  id<-data.risk.agg$id
  s_h<-data.risk.agg$s_h
  vt<-data.risk.agg$N.counts
  vm<-data.risk.agg$S_z
  n0<-length(vt)
  index<-length(vt)
  nk<-length(index)
  listij<-c()
  cnt=0
  cnt1=0
  for(k in 1:index){
    for(j in k:n0){
      if (vt[j]>vt[k]){
        cnt=cnt+1
        cnt1=cnt1+as.numeric(vm[j]>vm[k])
      }
    }
  }
  return(list(n_h=cnt,n.s_h=cnt1))
}

## For variance estimation: using perturbation
My.counts.pb<-function(ss,data.risk){
  data.risk<-data.risk[data.risk$t.obs>ss,]
  data.risk$s_h<-ss
  data.risk$indicator<-(data.risk$t<=apply(data.risk[,c('s_h','t.obs')],1,min) & data.risk$Num!=0)
  data.risk.agg<-data.risk %>% group_by(id) %>% mutate(N.counts = sum(indicator))
  data.risk.agg<-aggregate(data.risk.agg[,-which(names(data.risk.agg) %in% c('t','indicator','Num'))], by=list(ID=data.risk.agg$id), FUN=unique)
  
  data.risk.agg<-data.risk.agg[order(data.risk.agg$N.counts),]
  id<-data.risk.agg$id
  s_h<-data.risk.agg$s_h
  vt<-data.risk.agg$N.counts
  vm<-data.risk.agg$S_z
  w<-data.risk.agg$w
  n0<-length(vt)
  index<-length(vt)
  nk<-length(index)
  listij<-c()
  cnt=0
  cnt1=0
  for(k in 1:index){
    for(j in k:n0){
      if (vt[j]>vt[k]){
        cnt=cnt+1*w[j]*w[k]
        cnt1=cnt1+as.numeric(vm[j]>vm[k])*w[j]*w[k]
      }
    }
  }
  return(list(n_h=cnt,n.s_h=cnt1))
}

#############################################
## process to likelihood ####################
#############################################
alpha.estimate<-function(tableij,s_h,max.fu){
  n_h<-lapply(tableij, function(x) x$n_h)
  n.s_h<-lapply(tableij, function(x) x$n.s_h)
  score = cbind(1,(s_h[2:length(s_h)]-max.fu/2)/(max.fu/2),((s_h[2:length(s_h)]-max.fu/2)/(max.fu/2))^2)
  totalnum = unlist(n_h)
  eventsnum = unlist(n.s_h)
  opt.model <- neldermead(c(0,0,0), fn =My.loglikehood,score=score,totalnum=totalnum,eventsnum=eventsnum)
  return(opt.model$par)
}

#############################################
## pseudo partial-likelihood ################
#############################################
My.loglikehood<-function(beta,score,totalnum,eventsnum){
  g <- score%*%beta
  lnL <- -sum(t(eventsnum)%*%g- t(totalnum)%*%log(1+exp(g)))
  return(c(lnL))
}


