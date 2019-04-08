


#Set varibles to make time points easier to spot
time0<-1
time45<-2
time90<-3

#Input raw data
qpcr<-list()
qpcr[[time0]] <- read.table("qpcr-0.txt", header=T,sep="\t",as.is=T)
qpcr[[time45]] <- read.table("qpcr-45.txt", header=T,sep="\t",as.is=T)
qpcr[[time90]] <- read.table("qpcr-90.txt", header=T,sep="\t",as.is=T)

for (time_point in time0:time90)
{
  qpcr[[time_point]][qpcr[[time_point]]=="No Ct"]<-NA
}


#Prepare lists
ct_values<-list()
ct_values_means<-list()
d_ct_values_means<-list()
input_values<-list()
input_mean<-list()
d_input_mean<-list()
dd_ct_values_means<-list()


for (time_point in time0:time90)
{
  
  
  #Create object
  ct_values[[time_point]]<- matrix(NA,0,2)
  
  #Add reps 1:8 to object
  for (rep in 1:8) {
    y_coord <- rep*12-11
    ct_values[[time_point]] <-rbind(
      ct_values[[time_point]],
      cbind(as.numeric(qpcr[[time_point]]$Ct..dR[y_coord:(y_coord+4)]),
            as.numeric(qpcr[[time_point]]$Ct..dR[(y_coord+5):(y_coord+9)]
            ))
    )
  }
  
  
  
  #Average techical reps from same isogenic replicate
  ct_values_means[[time_point]] <- matrix(NA,0,2)
  for(rep in 1:8){
    y_coord <- rep*5-4
    ct_values_means[[time_point]] <- rbind(ct_values_means[[time_point]], colMeans(ct_values[[time_point]][y_coord:(y_coord+4),1:2]))
  }
  
  #Import Input from time_point mins

  
  input_values[[time_point]] <- cbind(
                                      as.numeric(qpcr[[time_point]]$Ct..dR[ 12*(0:4)+11]),
                                      as.numeric(qpcr[[time_point]]$Ct..dR[ 12*(0:4)+12])
                                      )
  
  #Calucate Average of technical reps.
  input_mean[[time_point]] <- colMeans(input_values[[time_point]],na.rm=T)
  
  if(time_point == 2)
  { #Plate layout is reversed wrt control and target primers. See lab book.
    d_ct_values_means[[time_point]]<- ct_values_means[[time_point]][,2]-ct_values_means[[time_point]][,1]
    d_input_mean[[time_point]]<-input_mean[[time_point]][2]-input_mean[[time_point]][1]
  
    dd_ct_values_means[[time_point]]<- d_ct_values_means[[time_point]]-d_input_mean[[time_point]]
  } else {
    d_ct_values_means[[time_point]]<- ct_values_means[[time_point]][,1]-ct_values_means[[time_point]][,2]
    d_input_mean[[time_point]]<-input_mean[[time_point]][1]-input_mean[[time_point]][2]
    
    dd_ct_values_means[[time_point]]<- d_ct_values_means[[time_point]]-d_input_mean[[time_point]]
  }
}

df_ct <-data.frame(
          cbind(2^dd_ct_values_means[[time0]],
                2^dd_ct_values_means[[time45]], 
                2^dd_ct_values_means[[time90]]))
colnames(df_ct) <- c("0 min", "45 mins", "90 mins")

set.seed(1)
par(mar=c(5,5,3,1))
boxplot(df_ct,main="ERα occupany of the TFF1 promoter", range=1, ylab= expression(paste("Enrichment over Input, 2"^"ΔΔCt")),xlab="Time point", pch=NA)
stripchart(df_ct, add=T, vertical=T, pch=20,method="jitter", jitter=0.05)          
