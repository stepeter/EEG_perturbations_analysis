library(R.matlab)

pathname = "/Users/stepeter/Documents/Neuroimage_revisions/event_timing_data/Connectivity_median_plots/"
for(cond in c("pullsStn","pullsWalk")){
  for(plotBand_cort in c("theta","alpha")){
    netData_sbjs = readMat(paste(pathname,"netVals_",cond,"_aveTime_sbjs.mat",sep=""))
    
    head(netData_sbjs$netVals[1])
    df = as.data.frame(netData_sbjs$netVals)
    
    #Determine significant connections above 0
    pvals_theta = matrix(0L, nrow = 1, ncol = 16*16-16)
    pvals_alpha = matrix(0L, nrow = 1, ncol = 16*16-16)
    counter = 0
    counter_all = 0
    for(i in 1:16){
      for(j in 1:16){
        counter_all = counter_all+1
        if(i!=j){
          counter=counter+1
          pvals_theta[counter] = t.test(netData_sbjs$netVals[[1]][[counter_all]][[1]])$p.value
          pvals_alpha[counter] = t.test(netData_sbjs$netVals[[2]][[counter_all]][[1]])$p.value
        }
      }
    }
    
    #Reshape into 16 x 16
    pvals_theta_box = matrix(1L, nrow = 16, ncol = 16)
    pvals_alpha_box = matrix(1L, nrow = 16, ncol = 16)
    counter = 0
    for(i in 1:16){
      for(j in 1:16){
        if(i!=j){
          counter=counter+1
          pvals_theta_box[j,i] = pvals_theta[counter]
          pvals_alpha_box[j,i] = pvals_alpha[counter]
        }
      }
    }
    
    #FDR correction
    # for(i in 1:16){
    #   pvals_theta_box[,i]=p.adjust(pvals_theta_box[,i],method='fdr')
    #   pvals_alpha_box[,i]=p.adjust(pvals_alpha_box[,i],method='fdr')
    # }
    
    #Check for significance
    isSig_theta_box = matrix(0L, nrow = 16, ncol = 16)
    isSig_alpha_box = matrix(0L, nrow = 16, ncol = 16)
    counter = 0
    for(i in 1:16){
      for(j in 1:16){
        if(i!=j){
          counter=counter+1
          if(pvals_theta_box[i,j]<0.05){
            isSig_theta_box[i,j]=1
          }
          if(pvals_alpha_box[i,j]<0.05){
            isSig_alpha_box[i,j]=1
          }
        }
      }
    }
    
    #Load average data and mask it by significance
    netData_ave = readMat(paste(pathname,"netVals_",cond,"_aveTime.mat",sep=""))
    thetaAve = netData_ave$netVals.ave[[1]]
    alphaAve = netData_ave$netVals.ave[[2]]
    thetaAve[isSig_theta_box==0]=0
    alphaAve[isSig_alpha_box==0]=0
    
    #Chord plot (for cortico-cortical data)
    library(circlize)
    circos.clear()
    
    #Rearrange order in matrix
    exchangeInds = c(5,6,7,8,4,3,2,1)+8
    cortTheta = thetaAve[9:16,9:16]
    cortAlpha = alphaAve[9:16,9:16]
    for(i in 1:8){
      for(j in 1:8){
        cortTheta[i,j] = thetaAve[exchangeInds[i],exchangeInds[j]]
        cortAlpha[i,j] = alphaAve[exchangeInds[i],exchangeInds[j]]
      }
    }
    
    if(plotBand_cort=='theta'){
      dat_in = cortTheta
    }else{
      dat_in = cortAlpha
    }
    png(paste(pathname,cond,"_",plotBand_cort,"_emgEMG.png",sep=""),bg="transparent",width = 720, height = 720)
    color_vals = c('gold','purple','green','cyan','red','blue','orange','magenta')#'#FF00FF','#FF6600','#0000FF','#FF0000','#D4AA00','#800080','#00FF00','#00FFFF')
    color_vals_links = adjustcolor(color_vals, alpha.f = 0.6)
    factors = c("RTA","RSO","RMG","RPL","LPL","LMG","LSO","LTA")
    circos.par(cell.padding = c(0, 0, 0, 0),start.degree = 85,gap.degree=c(0,0,0,10,0,0,0,10)) 
    sign_vals = sign(dat_in)
    A = abs(dat_in)*10^4
    circos.initialize(factors, xlim = cbind(c(0,0,0,0,0,0,0,0), 20*c(1,1,1,1,1,1,1,1))) 
    
    
    circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.05,
                           #panel.fun for each sector
                           panel.fun = function(x, y) {
                             offset = 90
                             #select details of current sector
                             name = get.cell.meta.data("sector.index")
                             i = get.cell.meta.data("sector.numeric.index")
                             xlim = get.cell.meta.data("xlim")
                             ylim = get.cell.meta.data("ylim")
                             
                             #text direction (dd) and adjusmtents (aa)
                             theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
                             dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
                             aa = c(1, 0.5)
                             if(theta < 90 || theta > 270)  aa = c(0, 0.5)
                             
                             #plot cortical labels
                             #circos.text(x=mean(xlim), y=1.7, labels=name, facing = dd, cex=0.6,  adj = aa)
                             
                             #plot main sector
                             circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2], ytop=ylim[2], 
                                         col = color_vals[i], border=color_vals[i])
                             
                             #plot axis
                             circos.axis(labels.cex=2.8, direction = "outside", major.at=c(0,5,10,15,21),#seq(from=0,to=floor(12),by=4), 
                                         minor.ticks=1, labels.away.percentage = 0.15, lwd=2,labels.niceFacing=FALSE)
                           }) 
    
    #Add nonsignificant links first
    curr_starts_gray=matrix(0L, nrow = 1, ncol = 8)
    for(i in 1:8) {
      for(j in 1:8){
        if((A[i,j]!=0 | A[j,i]!=0) & (A[j,i]>A[i,j])){
          if(sign_vals[i,j]==0){
            if(plotBand_cort=='theta'){
              data_temp = abs(netData_ave$netVals.ave[[1]])*10^4
              circos.link(factors[i], c(curr_starts_gray[i],curr_starts_gray[i]+data_temp[i,j]), factors[j], c(curr_starts_gray[j],curr_starts_gray[j]+data_temp[j,i]), col = adjustcolor('gray', alpha.f = 0.4)) 
            }else{
              data_temp = abs(netData_ave$netVals.ave[[2]])*10^4
              circos.link(factors[i], c(curr_starts_gray[i],curr_starts_gray[i]+data_temp[i,j]), factors[j], c(curr_starts_gray[j],curr_starts_gray[j]+data_temp[j,i]), col = adjustcolor('gray', alpha.f = 0.4)) 
            }
            curr_starts_gray[i] = curr_starts_gray[i]+data_temp[i,j]
            curr_starts_gray[j] = curr_starts_gray[j]+data_temp[j,i]
          }
        }
      }
    }
    
    curr_starts=matrix(0L, nrow = 1, ncol = 8)
    for(i in 1:8) {
      for(j in 1:8){
        if((A[i,j]!=0 | A[j,i]!=0) & (A[j,i]>A[i,j])){
          if(sign_vals[i,j]>0 & sign_vals[j,i]>0){
            circos.link(factors[i], c(curr_starts[i],curr_starts[i]+A[i,j]), factors[j], c(curr_starts[j],curr_starts[j]+A[j,i]), col = adjustcolor('red', alpha.f = 0.6))#color_vals_links[8]) 
            curr_starts[i] = curr_starts[i]+A[i,j]
            curr_starts[j] = curr_starts[j]+A[j,i]
          }else if(sign_vals[i,j]<0 & sign_vals[j,i]<0){
            circos.link(factors[i], c(curr_starts[i],curr_starts[i]+A[i,j]), factors[j], c(curr_starts[j],curr_starts[j]+A[j,i]), col = adjustcolor('blue', alpha.f = 0.6)) 
            curr_starts[i] = curr_starts[i]+A[i,j]
            curr_starts[j] = curr_starts[j]+A[j,i]
          }else if(sign_vals[i,j]!=0 & sign_vals[j,i]!=0){
            circos.link(factors[i], c(curr_starts[i],curr_starts[i]+A[i,j]), factors[j], c(curr_starts[j],curr_starts[j]+A[j,i]), col = adjustcolor('magenta', alpha.f = 0.6)) 
            curr_starts[i] = curr_starts[i]+A[i,j]
            curr_starts[j] = curr_starts[j]+A[j,i]
          }
        }
      }
    }
    dev.off()
  }
}
